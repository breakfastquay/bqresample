/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*- vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#ifndef BQ_RESAMPLER_H
#define BQ_RESAMPLER_H

#include "bqvec/Restrict.h"

namespace breakfastquay {

class Resampler
{
public:
    enum Quality { Best, FastestTolerable, Fastest };
    enum Exception { ImplementationError };

    /**
     * Construct a resampler with the given quality level and channel
     * count.  maxBufferSize gives a bound on the maximum incount size
     * that may be passed to the resample function before the
     * resampler needs to reallocate its internal buffers.
     */
    Resampler(Quality quality, int channels, int maxBufferSize = 0,
              int debugLevel = 0);
    ~Resampler();

    /**
     * Resample the given multi-channel buffers, where incount is the
     * number of frames in the input buffers.  Returns the number of
     * frames written to the output buffers.
     */
    int resample(const float *const BQ_R__ *const BQ_R__ in,
                 float *const BQ_R__ *const BQ_R__ out,
                 int incount,
                 float ratio,
                 bool final = false);

    /**
     * Resample the given interleaved buffer, where incount is the
     * number of frames in the input buffer (i.e. it has incount *
     * getChannelCount() samples).  Returns the number of frames
     * written to the output buffer.
     */
    int resampleInterleaved(const float *const BQ_R__ in,
                            float *const BQ_R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int getChannelCount() const;

    void reset();

    class Impl;

protected:
    Impl *d;
    int m_method;
};

}

#endif
