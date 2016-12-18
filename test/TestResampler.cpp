/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "bqresample/Resampler.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;
using namespace breakfastquay;

BOOST_AUTO_TEST_SUITE(TestResampler)

#define LEN(a) (int(sizeof(a)/sizeof(a[0])))

static vector<float>
sine(double samplerate, double frequency, int nsamples)
{
    vector<float> v(nsamples, 0.f);
    for (int i = 0; i < nsamples; ++i) {
        v[i] = sin ((i * 2.0 * M_PI * frequency) / samplerate);
    }
    return v;
}

#define COMPARE_N(a, b, n)                    \
    for (int cmp_i = 0; cmp_i < n; ++cmp_i) { \
        BOOST_CHECK_SMALL((a)[cmp_i] - (b)[cmp_i], 1e-4f);      \
    }

static const float guard_value = -999.f;

BOOST_AUTO_TEST_CASE(interpolated_sine)
{
    // Interpolating a sinusoid should give us a sinusoid, once we've
    // dropped the first few samples
    vector<float> in = sine(8, 2, 1000); // 2Hz wave at 8Hz: [ 0, 1, 0, -1 ] etc
    vector<float> expected = sine(16, 2, 2000);
    vector<float> out(in.size() * 2 + 1, guard_value);
    Resampler r(Resampler::Parameters(), 1);
    r.resampleInterleaved(out.data(), out.size(), in.data(), in.size(), 2, true);
    const float *outf = out.data() + 200, *expectedf = expected.data() + 200;
    COMPARE_N(outf, expectedf, 600);
    // should have an exact number of output samples
    BOOST_CHECK_NE(out[out.size()-2], guard_value);
    BOOST_CHECK_EQUAL(out[out.size()-1], guard_value);
}

BOOST_AUTO_TEST_CASE(overrun_interleaved)
{
    // Check that the outcount argument is correctly used: any samples
    // already in the out buffer beyond outcount*channels must be left
    // untouched. We test this with short buffers (likely to be
    // shorter than the resampler filter length) and longer ones, with
    // resampler ratios that reduce, leave unchanged, and raise the
    // sample rate, and with all quality settings.

    // Options are ordered from most normal/likely to least.
    
    int channels = 2;

    int lengths[] = { 6000, 6 };
    int constructionBufferSizes[] = { 0, 1000 };
    double ratios[] = { 1.0, 10.0, 0.1 };
    Resampler::Quality qualities[] = {
        Resampler::FastestTolerable, Resampler::Best, Resampler::Fastest
    };

    bool failed = false;
    
    for (int li = 0; li < LEN(lengths); ++li) {
        for (int cbi = 0; cbi < LEN(constructionBufferSizes); ++cbi) {
            for (int ri = 0; ri < LEN(ratios); ++ri) {
                for (int qi = 0; qi < LEN(qualities); ++qi) {
                    
                    int length = lengths[li];
                    double ratio = ratios[ri];
                    Resampler::Parameters parameters;
                    parameters.quality = qualities[qi];
                    parameters.maxBufferSize = constructionBufferSizes[cbi];
                    Resampler r(parameters, channels);

                    float *inbuf = new float[length * channels];
                    for (int i = 0; i < length; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            inbuf[i*channels+c] =
                                sinf((i * 2.0 * M_PI * 440.0) / 44100.0);
                        }
                    }

                    int outcount = int(round(length * ratio));
                    outcount -= 10;
                    if (outcount < 1) outcount = 1;
                    int outbuflen = outcount + 10;
                    float *outbuf = new float[outbuflen * channels];
                    for (int i = outcount; i < outbuflen; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            outbuf[i*channels+c] = guard_value;
                        }
                    }
/*
                    cerr << "\nTesting with length = " << length << ", ratio = "
                         << ratio << ", outcount = " << outcount << ", final = false"
                         << endl;
*/                    
                    r.resampleInterleaved
                        (outbuf, outcount, inbuf, length, ratio, false);

                    for (int i = outcount; i < outbuflen; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            BOOST_CHECK_EQUAL(outbuf[i*channels+c], guard_value);
                            if (outbuf[i*channels+c] != guard_value) {
                                failed = true;
                            }
                        }
                    }

                    if (failed) {
                        cerr << "Test failed, abandoning remaining loop cycles"
                             << endl;
                        break;
                    }
/*
                    cerr << "\nContinuing with length = " << length << ", ratio = "
                         << ratio << ", outcount = " << outcount << ", final = true"
                         << endl;
*/
                    r.resampleInterleaved
                        (outbuf, outcount, inbuf, length, ratio, true);

                    for (int i = outcount; i < outbuflen; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            BOOST_CHECK_EQUAL(outbuf[i*channels+c], guard_value);
                            if (outbuf[i*channels+c] != guard_value) {
                                failed = true;
                            }
                        }
                    }

                    delete[] outbuf;
                    delete[] inbuf;

                    if (failed) {
                        cerr << "Test failed, abandoning remaining loop cycles"
                             << endl;
                        break;
                    }
                }

                if (failed) break;
            }
            if (failed) break;
        }
        if (failed) break;
    }
}

BOOST_AUTO_TEST_SUITE_END()

