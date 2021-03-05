/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    bqresample

    A small library wrapping various audio sample rate conversion
    implementations in C++.

    Copyright 2007-2021 Particular Programs Ltd.

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use, copy,
    modify, merge, publish, distribute, sublicense, and/or sell copies
    of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
    ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
    CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

    Except as contained in this notice, the names of Chris Cannam and
    Particular Programs Ltd shall not be used in advertising or
    otherwise to promote the sale, use or other dealings in this
    Software without prior written authorization.
*/

#include "BQResampler.h"

#include <cmath>

#include <iostream>

#include <bqvec/Allocators.h>
#include <bqvec/VectorOps.h>

//!!! copy/assign etc
//!!! channels

using std::vector;
using std::cerr;
using std::endl;
using std::min;
using std::max;

namespace breakfastquay {

BQResampler::BQResampler(Parameters parameters) :
    m_qparams(parameters.quality),
    m_dynamism(parameters.dynamism),
    m_ratio_change(parameters.ratioChange),
    m_debug_level(parameters.debugLevel),
    m_initial_rate(parameters.referenceSampleRate),
    m_fade_count(0),
    m_initialised(false)
{
    if (m_debug_level > 0) {
        cerr << "BQResampler::BQResampler: "
             << (m_dynamism == RatioOftenChanging ? "often-changing" : "mostly-fixed")
             << ", "
             << (m_ratio_change == SmoothRatioChange ? "smooth" : "sudden")
             << " ratio changes, ref " << m_initial_rate << " Hz" << endl;
    }
    
    if (m_dynamism == RatioOftenChanging) {
        m_proto_length = m_qparams.proto_p * m_qparams.p_multiple + 1;
        if (m_debug_level > 0) {
            cerr << "BQResampler: creating prototype filter of length "
                 << m_proto_length << endl;
        }
        m_prototype = make_filter(m_proto_length, m_qparams.proto_p);
        m_prototype.push_back(0.0); // interpolate without fear
    }

    int phase_reserve = 2 * int(round(m_initial_rate));
    int buffer_reserve = 1000;
    m_state_a.phase_info.reserve(phase_reserve);
    m_state_a.buffer.reserve(buffer_reserve);

    if (m_dynamism == RatioOftenChanging) {
        m_state_b.phase_info.reserve(phase_reserve);
        m_state_b.buffer.reserve(buffer_reserve);
    }

    m_s = &m_state_a;
    m_fade = &m_state_b;
}

BQResampler::QualityParams::QualityParams(Quality q)
{
    switch (q) {
    case Fastest:
        p_multiple = 12;
        proto_p = 160;
        k_snr = 70.0;
        k_transition = 0.2;
        cut = 0.9;
        break;
    case FastestTolerable:
        p_multiple = 62;
        proto_p = 160;
        k_snr = 90.0;
        k_transition = 0.05;
        cut = 0.975;
        break;
    case Best:
        p_multiple = 122;
        proto_p = 800;
        k_snr = 100.0;
        k_transition = 0.01;
        cut = 0.995;
        break;
    }
}
    
int
BQResampler::resample(float *const out,
                      int outspace,
                      const float *const in,
                      int incount,
                      double ratio,
                      bool final) {

    int fade_length = round(m_initial_rate / 1000.0);
    if (fade_length < 6) {
        fade_length = 6;
    }
    int max_fade = min(outspace, int(floor(incount * ratio))) / 2;
    if (fade_length > max_fade) {
        fade_length = max_fade;
    }
        
    if (!m_initialised) {
        state_for_ratio(*m_s, ratio, *m_fade);
        m_initialised = true;
    } else if (ratio != m_s->parameters.ratio) {
        state *tmp = m_fade;
        m_fade = m_s;
        m_s = tmp;
        state_for_ratio(*m_s, ratio, *m_fade);
        if (m_ratio_change == SmoothRatioChange) {
            if (m_debug_level > 0) {
                cerr << "BQResampler: ratio changed, beginning fade of length "
                     << fade_length << endl;
            }
            m_fade_count = fade_length;
        }
    }

    int i = 0, o = 0;
    int bufsize = m_s->buffer.size();

    while (o < outspace) {
        while (i < incount && m_s->fill < bufsize) {
            m_s->buffer[m_s->fill++] = in[i++];
        }
        if (m_s->fill == bufsize) {
            out[o++] = reconstruct_one(m_s);
        } else if (final && m_s->fill > m_s->centre) {
            out[o++] = reconstruct_one(m_s);
        } else if (final && m_s->fill == m_s->centre &&
                   m_s->current_phase != m_s->initial_phase) {
            out[o++] = reconstruct_one(m_s);
        } else {
            break;
        }
    }

    int fbufsize = m_fade->buffer.size();
    int fi = 0, fo = 0;
    while (fo < o && m_fade_count > 0) {
        while (fi < incount && m_fade->fill < fbufsize) {
            m_fade->buffer[m_fade->fill++] = in[fi++];
        }
        if (m_fade->fill == bufsize) {
            double r = reconstruct_one(m_fade);
            double fadeWith = out[fo];
            double extent = double(m_fade_count - 1) / double(fade_length);
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

int
BQResampler::gcd(int a, int b) const
{
    int c = a % b;
    if (c == 0) return b;
    else return gcd(b, c);
}

double
BQResampler::bessel0(double x) const
{
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

vector<double>
BQResampler::kaiser(double beta, int len) const
{
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

void
BQResampler::kaiser_params(double attenuation,
                           double transition,
                           double &beta,
                           int &len) const
{
    if (attenuation > 21.0) {
        len = 1 + ceil((attenuation - 7.95) / (2.285 * transition));
    } else {
        len = 1 + ceil(5.79 / transition);
    }
    beta = 0.0;
    if (attenuation > 50.0) {
        beta = 0.1102 * (attenuation - 8.7);
    } else if (attenuation > 21.0) {
        beta = 0.5842 * (pow (attenuation - 21.0, 0.4)) +
            0.07886 * (attenuation - 21.0);
    }
}

vector<double>
BQResampler::kaiser_for(double attenuation,
                        double transition,
                        int minlen,
                        int maxlen) const
{
    double beta;
    int m;
    kaiser_params(attenuation, transition, beta, m);
    int mb = m;
    if (maxlen > 0 && mb > maxlen - 1) {
        mb = maxlen - 1;
    } else if (minlen > 0 && mb < minlen) {
        mb = minlen;
    }
    if (mb % 2 == 0) ++mb;
    if (m_debug_level > 0) {
        cerr << "BQResampler: window attenuation " << attenuation
             << ", transition " << transition
             << " -> length " << m << " adjusted to " << mb
             << ", beta " << beta << endl;
    }
    return kaiser(beta, mb);
}
    
void
BQResampler::sinc_multiply(double peak_to_zero, vector<double> &buf) const
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

BQResampler::params
BQResampler::fill_params(double ratio, int num, int denom) const
{
    params p;
    int g = gcd (num, denom);
    p.ratio = ratio;
    p.numerator = num / g;
    p.denominator = denom / g;
    p.effective = double(p.numerator) / double(p.denominator);
    p.peak_to_zero = max(p.denominator, p.numerator);
    p.peak_to_zero /= m_qparams.cut;
    p.scale = double(p.numerator) / double(p.peak_to_zero);

    if (m_debug_level > 0) {
        cerr << "BQResampler: ratio " << p.ratio
             << " -> fraction " << p.numerator << "/" << p.denominator
             << " with error " << p.effective - p.ratio
             << endl;
        cerr << "BQResampler: peak-to-zero " << p.peak_to_zero
             << ", scale " << p.scale
             << endl;
    }
    
    return p;
}
    
BQResampler::params
BQResampler::pick_params(double ratio) const
{
    // Farey algorithm, see
    // https://www.johndcook.com/blog/2010/10/20/best-rational-approximation/
    int max_denom = 192000;
    double a = 0.0, b = 1.0, c = 1.0, d = 0.0;
    double pa = a, pb = b, pc = c, pd = d;
    double eps = 1e-9;
    while (b <= max_denom && d <= max_denom) {
        double mediant = (a + c) / (b + d);
        if (fabs(ratio - mediant) < eps) {
            if (b + d <= max_denom) {
                return fill_params(ratio, a + c, b + d);
            } else if (d > b) {
                return fill_params(ratio, c, d);
            } else {
                return fill_params(ratio, a, b);
            }
        }
        if (ratio > mediant) {
            pa = a; pb = b;
            a += c; b += d;
        } else {
            pc = c; pd = d;
            c += a; d += b;
        }
    }
    if (fabs(ratio - (pc / pd)) < fabs(ratio - (pa / pb))) {
        return fill_params(ratio, pc, pd);
    } else {
        return fill_params(ratio, pa, pb);
    }
}

void
BQResampler::phase_data_for(vector<BQResampler::phase_rec> &target_phase_data,
                            floatbuf &target_phase_sorted_filter,
                            int filter_length,
                            const vector<double> *filter,
                            int initial_phase,
                            int input_spacing,
                            int output_spacing) const
{
    target_phase_data.clear();
    target_phase_data.reserve(input_spacing);
        
    for (int p = 0; p < input_spacing; ++p) {
        int next_phase = p - output_spacing;
        while (next_phase < 0) next_phase += input_spacing;
        next_phase %= input_spacing;
        double dspace = double(input_spacing);
        int zip_length = ceil(double(filter_length - p) / dspace);
        int drop = ceil(double(max(0, output_spacing - p)) / dspace);
        phase_rec phase;
        phase.next_phase = next_phase;
        phase.drop = drop;
        phase.length = zip_length;
        phase.start_index = 0; // we fill this in below if needed
        target_phase_data.push_back(phase);
    }

    if (m_dynamism == RatioMostlyFixed) {
        if (!filter) throw std::logic_error("filter required at phase_data_for in RatioMostlyFixed mode");
        target_phase_sorted_filter.clear();
        target_phase_sorted_filter.reserve(filter_length);
        for (int p = initial_phase; ; ) {
            phase_rec &phase = target_phase_data[p];
            phase.start_index = target_phase_sorted_filter.size();
            for (int i = 0; i < phase.length; ++i) {
                target_phase_sorted_filter.push_back
                    ((*filter)[i * input_spacing + p]);
            }
            p = phase.next_phase;
            if (p == initial_phase) {
                break;
            }
        }
    }
}

vector<double>
BQResampler::make_filter(int filter_length, double peak_to_zero) const
{
    vector<double> filter;
    filter.reserve(filter_length);

    vector<double> kaiser = kaiser_for(m_qparams.k_snr, m_qparams.k_transition,
                                       1, filter_length);
    int k_length = kaiser.size();

    if (k_length == filter_length) {
        sinc_multiply(peak_to_zero, kaiser);
        return kaiser;
    } else {
        kaiser.push_back(0.0);
        double m = double(k_length - 1) / double(filter_length - 1);
        for (int i = 0; i < filter_length; ++i) {
            double ix = i * m;
            int iix = floor(ix);
            double remainder = ix - iix;
            double value = 0.0;
            value += kaiser[iix] * (1.0 - remainder);
            value += kaiser[iix+1] * remainder;
            filter.push_back(value);
        }
        sinc_multiply(peak_to_zero, filter);
        return filter;
    }
}

void
BQResampler::state_for_ratio(BQResampler::state &target_state,
                             double ratio,
                             const BQResampler::state &BQ_R__ prev_state) const
{
    params parameters = pick_params(ratio);
    target_state.parameters = parameters;
    
    target_state.filter_length =
        int(parameters.peak_to_zero * m_qparams.p_multiple + 1);

    if (target_state.filter_length % 2 == 0) {
        ++target_state.filter_length;
    }

    int half_length = target_state.filter_length / 2; // nb length is odd
    int input_spacing = parameters.numerator;
    int initial_phase = half_length % input_spacing;

    target_state.initial_phase = initial_phase;
    target_state.current_phase = initial_phase;

    if (m_dynamism == RatioMostlyFixed) {
        
        if (m_debug_level > 0) {
            cerr << "BQResampler: creating filter of length "
                 << target_state.filter_length << endl;
        }

        vector<double> filter =
            make_filter(target_state.filter_length, parameters.peak_to_zero);

        phase_data_for(target_state.phase_info,
                       target_state.phase_sorted_filter,
                       target_state.filter_length, &filter,
                       target_state.initial_phase,
                       input_spacing,
                       parameters.denominator);
    } else {
        phase_data_for(target_state.phase_info,
                       target_state.phase_sorted_filter,
                       target_state.filter_length, 0,
                       target_state.initial_phase,
                       input_spacing,
                       parameters.denominator);
    }

    int buffer_left = half_length / input_spacing;
    int buffer_right = buffer_left + 1;

    int buffer_length = buffer_left + buffer_right;
    buffer_length = max(buffer_length, int(prev_state.buffer.size()));
    
    target_state.centre = buffer_length / 2;
    target_state.left = target_state.centre - buffer_left;
    target_state.fill = target_state.centre;

    int n_phases = int(target_state.phase_info.size());

    if (m_debug_level > 0) {
        cerr << "BQResampler: buffer left " << buffer_left
             << ", right " << buffer_right
             << ", total " << buffer_length << endl;
    
        cerr << "BQResampler: input spacing " << input_spacing
             << ", output spacing " << parameters.denominator
             << ", initial phase " << initial_phase
             << " of " << n_phases << endl;
    }

    if (prev_state.buffer.size() > 0) {
        if (int(prev_state.buffer.size()) == buffer_length) {
            target_state.buffer = prev_state.buffer;
            target_state.fill = prev_state.fill;
        } else {
            target_state.buffer = floatbuf(buffer_length, 0.0);
            for (int i = 0; i < prev_state.fill; ++i) {
                int offset = i - prev_state.centre;
                int new_ix = offset + target_state.centre;
                if (new_ix >= 0 && new_ix < buffer_length) {
                    target_state.buffer[new_ix] = prev_state.buffer[i];
                    target_state.fill = new_ix + 1;
                }
            }
        }

        int phases_then = int(prev_state.phase_info.size());
        double distance_through =
            double(prev_state.current_phase) / double(phases_then);
        target_state.current_phase = round(n_phases * distance_through);
        if (target_state.current_phase >= n_phases) {
            target_state.current_phase = n_phases - 1;
        }
    } else {
        target_state.buffer = floatbuf(buffer_length, 0.0);
    }
}

double
BQResampler::reconstruct_one(state *s) const
{
    const phase_rec &pr = s->phase_info[s->current_phase];
    int phase_length = pr.length;
    double result = 0.0;

    if (m_dynamism == RatioMostlyFixed) {
        int phase_start = pr.start_index;
        result = v_multiply_and_sum
            (s->phase_sorted_filter.data() + phase_start,
             s->buffer.data() + s->left,
             phase_length);
    } else {
        double m = double(m_proto_length - 1) / double(s->filter_length - 1);
        for (int i = 0; i < phase_length; ++i) {
            double sample = s->buffer[s->left + i];
            int filter_index = i * s->parameters.numerator + s->current_phase;
            double proto_index = m * filter_index;
            int iix = floor(proto_index);
            double remainder = proto_index - iix;
            double filter_value = m_prototype[iix] * (1.0 - remainder);
            filter_value += m_prototype[iix+1] * remainder;
            result += filter_value * sample;
        }
    }

    if (pr.drop > 0) {
        v_move(s->buffer.data(), s->buffer.data() + pr.drop,
               int(s->buffer.size()) - pr.drop);
        for (int i = 0; i < pr.drop; ++i) {
            s->buffer[s->buffer.size() - pr.drop + i] = 0.0;
        }
        s->fill -= pr.drop;
    }

    s->current_phase = pr.next_phase;
    return result * s->parameters.scale;
}

}
