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

BOOST_AUTO_TEST_CASE(overrun_interleaved)
{
    // Check that the outcount argument is correctly used: any samples
    // already in the out buffer beyond outcount*channels must be left
    // untouched. We test this with short buffers (likely to be
    // shorter than the resampler filter length) and longer ones, with
    // resampler ratios that reduce, leave unchanged, and raise the
    // sample rate, and with all quality settings.

    int channels = 2;

    int lengths[] = { 6, 6000 };
    int constructionBufferSizes[] = { 0, 1000 };
    double ratios[] = { 0.1, 1.0, 10.0 };
    Resampler::Quality qualities[] = {
        Resampler::Best, Resampler::FastestTolerable, Resampler::Fastest
    };

    for (int li = 0; li < LEN(lengths); ++li) {
        for (int cbi = 0; cbi < LEN(constructionBufferSizes); ++cbi) {
            for (int ri = 0; ri < LEN(ratios); ++ri) {
                for (int qi = 0; qi < LEN(qualities); ++qi) {
                    
                    int length = lengths[li];
                    int constructionBufferSize = constructionBufferSizes[cbi];
                    double ratio = ratios[ri];
                    Resampler::Quality quality = qualities[qi];
                    Resampler r(quality, channels, constructionBufferSize);

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
                            outbuf[i*channels+c] = -999.f;
                        }
                    }

                    r.resampleInterleaved
                        (outbuf, outcount, inbuf, length, ratio, false);

                    for (int i = outcount; i < outbuflen; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            BOOST_CHECK_EQUAL(outbuf[i*channels+c], -999.f);
                        }
                    }

                    r.resampleInterleaved
                        (outbuf, outcount, inbuf, length, ratio, true);

                    for (int i = outcount; i < outbuflen; ++i) {
                        for (int c = 0; c < channels; ++c) {
                            BOOST_CHECK_EQUAL(outbuf[i*channels+c], -999.f);
                        }
                    }

                    delete[] outbuf;
                    delete[] inbuf;
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

