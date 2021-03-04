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

#ifndef BQ_BQRESAMPLER_H
#define BQ_BQRESAMPLER_H

#include <vector>

#include <bqvec/Allocators.h>

namespace breakfastquay {

class BQResampler
{
public:
    enum Quality { Best, FastestTolerable, Fastest };
    enum Dynamism { RatioOftenChanging, RatioMostlyFixed };
    enum RatioChange { SmoothRatioChange, SuddenRatioChange };
    
    struct Parameters {
        Quality quality;
        Dynamism dynamism; 
        RatioChange ratioChange;
        double referenceSampleRate;
        int debugLevel;

        Parameters() :
            quality(FastestTolerable),
            dynamism(RatioMostlyFixed),
            ratioChange(SmoothRatioChange),
            referenceSampleRate(44100),
            debugLevel(0) { }
    };

    BQResampler(Parameters parameters);

    int resample(float *const out, int outspace,
                 const float *const in, int incount,
                 double ratio, bool final);

private:
    const Dynamism m_dynamism;
    const RatioChange m_ratio_change;
    const int m_debug_level;
    const double m_initial_rate;
    const int m_p_multiple;

    struct params {
        double ratio;
        int numerator;
        int denominator;
        double effective;
        double peak_to_zero;
        double scale;
        params() : ratio(1.0), numerator(1), denominator(1),
                   effective(1.0), peak_to_zero(0), scale(1.0) { }
    };

    struct phase_rec {
        int next_phase;
        int length;
        int start_index;
        int drop;
        phase_rec() : next_phase(0), length(0), start_index(0), drop(0) { }
    };

    typedef std::vector<float, breakfastquay::StlAllocator<float> > floatbuf;
    
    struct state {
        params parameters;
        int initial_phase;
        int current_phase;
        int filter_length;
        std::vector<phase_rec> phase_info;
        floatbuf phase_sorted_filter;
        floatbuf buffer;
        int left;
        int centre;
        int fill;
        state() : initial_phase(0), current_phase(0), filter_length(0),
                  left(0), centre(0), fill(0) { }
    };
    state m_s;
    state m_fading;
    int m_fade_count;
    
    std::vector<double> m_prototype;
    int m_proto_length;
    bool m_initialised;

    int gcd(int a, int b) const;
    double bessel0(double x) const;
    std::vector<double> kaiser(double beta, int len) const;
    std::vector<double> kaiser_for(double attenuation, double transition,
                                   int minlen, int maxlen) const;
    void sinc_multiply(double peak_to_zero, std::vector<double> &buf) const;

    params fill_params(double ratio, int num, int denom) const;
    params pick_params(double ratio) const;

    std::vector<double> make_filter(int filter_length,
                                    double peak_to_zero) const;
    
    std::vector<phase_rec> phase_data_for(int filter_length,
                                          const std::vector<double> &filter,
                                          floatbuf &phase_sorted_filter,
                                          int initial_phase,
                                          int input_spacing,
                                          int output_spacing) const;
    
    state state_for_ratio(double ratio) const;
    double reconstruct_one(state &s) const;
};

}

#endif
