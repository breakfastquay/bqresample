
#include <iostream>
#include <vector>

#include <cmath>

using namespace std;

class BQResampler
{
public:
    //!!! todo: quality settings
    //!!! todo: avoid reallocations in RatioOftenChanging mode (& document)
    //!!! copy/assign etc
    //!!! channels
    
    enum Dynamism { RatioOftenChanging, RatioMostlyFixed };
    
    BQResampler(Dynamism dynamism, double rate) :
        m_fade_count(0),
        m_dynamism(dynamism),
        m_initialised(false),
        m_initial_rate(rate),
        m_p_multiple(82)
    {
        int proto_p = 160;
        m_proto_length = proto_p * m_p_multiple + 1;
        double snr = 130.0;
//        double bandwidth = 80.0; // Hz
//        double transition = (bandwidth * 2.0 * M_PI) / m_initial_rate;
        double transition = 0.005;
        m_kaiser_prototype = kaiser_for
            (snr, transition, m_proto_length, m_proto_length);
        sinc_multiply(proto_p, m_kaiser_prototype);
        m_kaiser_prototype.push_back(0.0); // interpolate without fear
    }

    int resample(float *const out,
                 int outspace,
                 const float *const in,
                 int incount,
                 double ratio,
                 bool final) {

        int fade_length = round(m_initial_rate / 1000.0);
        if (fade_length < 6) fade_length = 6;
        int max_fade = min(outspace, int(floor(incount * ratio))) / 2;
        if (fade_length > max_fade) fade_length = max_fade;
        
        if (!m_initialised) {
            m_s = state_for_ratio(ratio);
            m_initialised = true;
        } else if (ratio != m_s.parameters.ratio) {
            m_fading = m_s;
            m_s = state_for_ratio(ratio);
            m_fade_count = fade_length;
        }

        int i = 0, o = 0;
        int bufsize = m_s.buffer.size();

        while (o < outspace) {
            while (i < incount && m_s.fill < bufsize) {
                m_s.buffer[m_s.fill++] = in[i++];
            }
            if (m_s.fill == bufsize) {
                out[o++] = reconstruct_one(m_s);
            } else if (final && m_s.fill > m_s.centre) {
                out[o++] = reconstruct_one(m_s);
            } else {
                break;
            }
        }

        int fbufsize = m_fading.buffer.size();
        int fi = 0, fo = 0;
        while (fo < o && m_fade_count > 0) {
            while (fi < incount && m_fading.fill < fbufsize) {
                m_fading.buffer[m_fading.fill++] = in[fi++];
            }
            if (m_fading.fill == bufsize) {
                double r = reconstruct_one(m_fading);
                double fadeWith = out[fo];
                double extent = double(m_fade_count-1) / double(fade_length);
                double mixture = 0.5 * (1.0 - cos(M_PI * extent));
                double mixed = r * mixture + fadeWith * (1.0 - mixture);
                out[fo] = mixed;
                ++fo;
                --m_fade_count;
            } else {
                break;
            }
        }
        
        return o;
    }

private:
    struct params {
        double ratio;
        int numerator;
        int denominator;
        double effective;
        int peak_to_zero;
        double scale;
        params() : ratio(1.0), numerator(1), denominator(1),
                   effective(1.0), peak_to_zero(0), scale(1.0) { }
    };

    struct phase_rec {
        int next_phase;
        int drop;
        int zip_length;
        vector<double> filter;
        phase_rec() : next_phase(0), drop(0), zip_length(0) { }
    };
    struct state {
        params parameters;
        int current_phase;
        int filter_length;
        vector<double> filter;
        vector<phase_rec> phases;
        vector<double> buffer;
        int left;
        int centre;
        int fill;
        state() : current_phase(0), filter_length(0),
                  left(0), centre(0), fill(0) { }
    };
    state m_s;
    state m_fading;
    int m_fade_count;

    Dynamism m_dynamism;
    bool m_initialised;
    double m_initial_rate;
    int m_p_multiple;
    vector<double> m_kaiser_prototype;
    int m_proto_length;
    
    static int gcd(int a, int b) {
        int c = a % b;
        if (c == 0) return b;
        else return gcd(b, c);
    }
    
    static params calc_rational(double ratio, int suggested_denom) {
        params p;
        int raw_num = round(ratio * suggested_denom);
        int g = gcd (raw_num, suggested_denom);
        p.ratio = ratio;
        p.numerator = raw_num / g;
        p.denominator = suggested_denom / g;
        p.effective = double(p.numerator) / double(p.denominator);
        p.peak_to_zero = max(p.denominator, p.numerator);
        if (ratio < 1.0) {
//            p.peak_to_zero /= std::max (0.99, ratio);
        } else { 
//            p.peak_to_zero /= 0.99;
        }
        p.scale = double(p.numerator) / double(p.peak_to_zero);
        return p;
    }

    params pick_params(double ratio) const {
        params best;
        bool first = true;
        static int candidates[] = { 44100, 96000 };
        static int n_cand = int(sizeof(candidates)/sizeof(candidates[0]));
        for (int i = 0; i < n_cand; ++i) {
            int c = candidates[i];
            params p = calc_rational(ratio, c);
            if (first) {
                best = p;
                first = false;
            } else if (fabs(p.effective - ratio) < fabs(best.effective - ratio)
                       && p.peak_to_zero <= best.peak_to_zero) {
                best = p;
            } else if (p.peak_to_zero < best.peak_to_zero
                       && fabs(p.effective - ratio) < 1e-5) {
                best = p;
            }
        }
        return best;
    }

    vector<phase_rec> phase_data_for(int filterlen, const vector<double> &filter,
                                     int input_spacing, int output_spacing) const {
        int length = filterlen;
        vector<phase_rec> phases;
        for (int p = 0; p < input_spacing; ++p) {
            int next_phase = p - output_spacing;
            while (next_phase < 0) next_phase += input_spacing;
            next_phase %= input_spacing;
            double dspace = double(input_spacing);
            int zip_length = ceil(double(length - p) / dspace);
            int drop = ceil(double(max(0, output_spacing - p)) / dspace);
            phase_rec phase;
            phase.next_phase = next_phase;
            phase.drop = drop;
            phase.zip_length = zip_length;
            if (m_dynamism == RatioMostlyFixed) {
                vector<double> pfilt(zip_length, 0.0);
                for (int i = 0; i < zip_length; ++i) {
                    pfilt[i] = filter[i * input_spacing + p];
                }
                phase.filter = pfilt;
            }
            phases.push_back(phase);
        }
        return phases;
    }
    
    state state_for_ratio(double ratio) const {

        params parameters = pick_params(ratio);

        state s;
        s.parameters = parameters;
        s.filter_length = int(parameters.peak_to_zero * m_p_multiple + 1);

        if (m_dynamism == RatioMostlyFixed) {
            //!!! can we share the prototype among instances even?
            vector<double> filter(s.filter_length, 0.0);
            double ratio =
                double(m_proto_length - 1) / double(s.filter_length - 1);
            for (int i = 0; i < s.filter_length; ++i) {
                double ix = i * ratio;
                int iix = floor(ix);
                double remainder = ix - iix;
                double value = m_kaiser_prototype[iix] * (1.0 - remainder);
                value += m_kaiser_prototype[iix+1] * remainder;
                filter[i] = value;
            }
            s.filter = filter;
        }

        int half_length = s.filter_length / 2; // nb length is actually odd
        int input_spacing = parameters.numerator;
        int initial_phase = half_length % input_spacing;
        int buffer_left = half_length / input_spacing;
        int buffer_right = buffer_left + 1;

        int buffer_length = buffer_left + buffer_right;
        buffer_length = max(buffer_length, int(m_s.buffer.size()));
        
        s.current_phase = initial_phase;
        s.phases = phase_data_for
            (s.filter_length, s.filter, input_spacing, parameters.denominator);
        s.buffer = vector<double>(buffer_length, 0.0);
        s.centre = buffer_length / 2;
        s.left = s.centre - buffer_left;
        s.fill = s.centre;

        int n_phases = int(s.phases.size());
        
        if (m_s.buffer.size() > 0) {
            if (m_s.current_phase < n_phases) {
                double distance_through =
                    double(m_s.current_phase) / double(n_phases);
                s.current_phase = round(n_phases * distance_through);
                if (s.current_phase >= n_phases) {
                    cerr << "!!! -> Need to drop an input sample!?" << endl;
                    s.current_phase = 0;
                }
            }
            for (int i = 0; i < m_s.fill; ++i) {
                int offset = i - m_s.centre;
                int new_ix = offset + s.centre;
                if (new_ix >= 0 && new_ix < int(s.buffer.size())) {
                    s.buffer[new_ix] = m_s.buffer[i];
                    s.fill = new_ix + 1;
                }
            }
        }
        
        return s;
    }

    double reconstruct_one(state &s) const {
        const phase_rec &pr = s.phases[s.current_phase];
        int phase_length = pr.zip_length;
        double result = 0.0;

        if (m_dynamism == RatioMostlyFixed) {
            for (int i = 0; i < phase_length; ++i) {
                result += pr.filter[i] * s.buffer[s.left + i];
            }
        } else {
            double ratio =
                double(m_proto_length - 1) / double(s.filter_length - 1);
            for (int i = 0; i < phase_length; ++i) {
                double sample = s.buffer[s.left + i];
                int filter_index = i * s.parameters.numerator + s.current_phase;
                double proto_index = ratio * filter_index;
                int iix = floor(proto_index);
                double remainder = proto_index - iix;
                double filter_value = m_kaiser_prototype[iix] * (1.0 - remainder);
                filter_value += m_kaiser_prototype[iix+1] * remainder;
                result += filter_value * sample;
            }
        }
        
        for (int i = pr.drop; i < int(s.buffer.size()); ++i) {
            s.buffer[i - pr.drop] = s.buffer[i];
        }
        for (int i = 0; i < pr.drop; ++i) {
            s.buffer[s.buffer.size() - pr.drop + i] = 0.0;
        }
        s.fill -= pr.drop;
        s.current_phase = pr.next_phase;
        return result * s.parameters.scale;
    }

    static double bessel0(double x) {
        static double facsquared[] = {
            0.0, 1.0, 4.0, 36.0,
            576.0, 14400.0, 518400.0, 25401600.0,
            1625702400.0, 131681894400.0, 1.316818944E13, 1.59335092224E15,
            2.29442532803E17, 3.87757880436E19, 7.60005445655E21,
            1.71001225272E24, 4.37763136697E26, 1.26513546506E29,
            4.09903890678E31, 1.47975304535E34
        };
        static int nterms = sizeof(facsquared) / sizeof(facsquared[0]);
        double b = 1.0;
        for (int n = 1; n < nterms; ++n) {
            double ff = facsquared[n];
            double term = pow(x / 2.0, n * 2.0) / ff;
            b += term;
        }
        return b;
    }

    static vector<double> kaiser(double beta, int len) {
        double denominator = bessel0(beta);
        int half = (len % 2 == 0 ? len/2 : (len+1)/2);
        vector<double> v(len, 0.0);
        for (int n = 0; n < half; ++n) {
            double k = (2.0 * n) / (len-1) - 1.0;
            v[n] = bessel0 (beta * sqrt(1.0 - k*k)) / denominator;
        }
        for (int n = half; n < len; ++n) {
            v[n] = v[len-1 - n];
        }
        return v;
    }

    static vector<double> kaiser_for(double attenuation,
                                     double transition,
                                     int minlen,
                                     int maxlen) {
        int m;
        if (attenuation > 21.0) {
            m = 1 + ceil((attenuation - 7.95) / (2.285 * transition));
        } else {
            m = 1 + ceil(5.79 / transition);
        }
        int mb = m;
        if (maxlen > 0 && mb > maxlen - 1) {
            mb = maxlen - 1;
        } else if (minlen > 0 && mb < minlen) {
            mb = minlen;
        }
        if (mb % 2 == 0) ++mb;
        double beta = 0.0;
        if (attenuation > 50.0) {
            beta = 0.1102 * (attenuation - 8.7);
        } else if (attenuation > 21.0) {
            beta = 0.5842 * (pow (attenuation - 21.0, 0.4)) +
                0.07886 * (attenuation - 21.0);
        }
        return kaiser(beta, mb);
    }
    
    static void sinc_multiply(double peak_to_zero, vector<double> &buf)
    {
        int len = int(buf.size());
        if (len < 2) return;

        int left = len / 2;
        int right = (len + 1) / 2;
        double m = M_PI / peak_to_zero;

        for (int i = 1; i <= right; ++i) {
            double x = i * m;
            double sinc = sin(x) / x;
            if (i <= left) {
                buf[left - i] *= sinc;
            }
            if (i < right) {
                buf[i + left] *= sinc;
            }
        }
    }
};

