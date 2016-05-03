/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    bqresample

    A small library wrapping various audio sample rate conversion
    implementations in C++.

    Copyright 2007-2015 Particular Programs Ltd.

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

#include "bqresample/Resampler.h"

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <algorithm>

#include "bqvec/Allocators.h"
#include "bqvec/VectorOps.h"

#ifdef HAVE_IPP
#include <ipps.h>
#include <ippsr.h>
#include <ippac.h>
#endif

#ifdef HAVE_SAMPLERATE
#define HAVE_LIBSAMPLERATE 1
#endif

#ifdef HAVE_LIBSAMPLERATE
#include <samplerate.h>
#endif

#ifdef HAVE_LIBRESAMPLE
#include <libresample.h>
#endif

#ifdef USE_SPEEX
#include "../speex/speex_resampler.h"
#endif

#ifndef HAVE_IPP
#ifndef HAVE_LIBSAMPLERATE
#ifndef HAVE_LIBRESAMPLE
#ifndef USE_SPEEX
#error No resampler implementation selected!
#endif
#endif
#endif
#endif

using namespace std;

namespace breakfastquay {

class Resampler::Impl
{
public:
    virtual ~Impl() { }
    
    virtual int resample(const float *const BQ_R__ *const BQ_R__ in, 
                         float *const BQ_R__ *const BQ_R__ out,
                         int incount,
                         double ratio,
                         bool final) = 0;
    
    virtual int resampleInterleaved(const float *const BQ_R__ in, 
                                    float *const BQ_R__ out,
                                    int incount,
                                    double ratio,
                                    bool final) = 0;

    virtual int getChannelCount() const = 0;

    virtual void reset() = 0;
};

namespace Resamplers {

#ifdef HAVE_IPP

class D_IPP : public Resampler::Impl
{
public:
    D_IPP(Resampler::Quality quality, int channels, int maxBufferSize,
          int debugLevel);
    ~D_IPP();

    int resample(const float *const BQ_R__ *const BQ_R__ in,
                 float *const BQ_R__ *const BQ_R__ out,
                 int incount,
                 double ratio,
                 bool final);

    int resampleInterleaved(const float *const BQ_R__ in,
                            float *const BQ_R__ out,
                            int incount,
                            double ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    IppsResamplingPolyphase_32f **m_state;
    float **m_inbuf;
    size_t m_inbufsz;
    float **m_outbuf;
    size_t m_outbufsz;
    int m_bufsize;
    int m_channels;
    int m_window;
    float m_factor;
    int m_history;
    int *m_lastread;
    double *m_time;
    int m_debugLevel;
    
    void setBufSize(int);
};

D_IPP::D_IPP(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_state(0),
    m_channels(channels),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        cerr << "Resampler::Resampler: using IPP implementation"
                  << endl;
    }

    int nStep;
    IppHintAlgorithm hint;

    switch (quality) {

    case Resampler::Best:
        m_window = 64;
        nStep = 80;
        hint = ippAlgHintAccurate;
        break;

    case Resampler::FastestTolerable:
//        m_window = 48;
        nStep = 16;
        m_window = 16;
//        nStep = 8;
        hint = ippAlgHintFast;
        break;

    case Resampler::Fastest:
        m_window = 24;
        nStep = 64;
        hint = ippAlgHintFast;
        break;
    }

    m_factor = 8; // initial upper bound on m_ratio, may be amended later
    m_history = int(m_window * 0.5 * max(1.0, 1.0 / m_factor)) + 1;

    m_state = new IppsResamplingPolyphase_32f *[m_channels];

    m_lastread = new int[m_channels];
    m_time = new double[m_channels];

    m_bufsize = maxBufferSize + m_history;

    if (m_debugLevel > 1) {
        cerr << "bufsize = " << m_bufsize << ", window = " << m_window << ", nStep = " << nStep << ", history = " << m_history << endl;
    }

    for (int c = 0; c < m_channels; ++c) {
        ippsResamplePolyphaseInitAlloc_32f(&m_state[c],
                                           float(m_window),
                                           nStep,
                                           0.95f,
                                           9.0f,
                                           hint);
        m_lastread[c] = m_history;
        m_time[c] = m_history;
    }

    m_inbufsz = m_bufsize + m_history + 2;
    if (m_debugLevel > 1) {
        cerr << "inbuf allocating " << m_bufsize << " + " << m_history << " + 2 = " << m_inbufsz << endl;
    }

    m_outbufsz = lrintf(ceil((m_bufsize - m_history) * m_factor + 2));
    if (m_debugLevel > 1) {
        cerr << "outbuf allocating (" << m_bufsize << " - " << m_history << ") * " << m_factor << " + 2 = " << m_outbufsz << endl;
    }

    m_inbuf  = allocate_and_zero_channels<float>(m_channels, m_inbufsz);
    m_outbuf = allocate_and_zero_channels<float>(m_channels, m_outbufsz);

    if (m_debugLevel > 1) {
        cerr << "Resampler init done" << endl;
    }
}

D_IPP::~D_IPP()
{
    for (int c = 0; c < m_channels; ++c) {
        ippsResamplePolyphaseFree_32f(m_state[c]);
    }

    deallocate_channels(m_inbuf, m_channels);
    deallocate_channels(m_outbuf, m_channels);

    delete[] m_lastread;
    delete[] m_time;
    delete[] m_state;
}

void
D_IPP::setBufSize(int sz)
{
    if (m_debugLevel > 1) {
        cerr << "resize bufsize " << m_bufsize << " -> ";
    }

    m_bufsize = sz;

    cerr << m_bufsize << endl;

    int n1 = m_bufsize + m_history + 2;
    int n2 = (int)lrintf(ceil((m_bufsize - m_history) * m_factor + 2));

    if (m_debugLevel > 1) {
        cerr << "(outbufsize = " << n2 << ")" << endl;
    }

    m_inbuf = reallocate_and_zero_extend_channels
        (m_inbuf, m_channels, m_inbufsz, m_channels, n1);

    m_outbuf = reallocate_and_zero_extend_channels
        (m_outbuf, m_channels, m_outbufsz, m_channels, n2);
            
    m_inbufsz = n1;
    m_outbufsz = n2;
}

int
D_IPP::resample(const float *const BQ_R__ *const BQ_R__ in,
                float *const BQ_R__ *const BQ_R__ out,
                int incount,
                double ratio,
                bool final)
{
    int outcount = 0;

    if (ratio > m_factor) {
        m_factor = ratio;
        m_history = int(m_window * 0.5 * max(1.0, 1.0 / m_factor)) + 1;
    }

    for (int c = 0; c < m_channels; ++c) {
        if (m_lastread[c] + incount + m_history > m_bufsize) {
            setBufSize(m_lastread[c] + incount + m_history);
        }
    }

    for (int c = 0; c < m_channels; ++c) {

        for (int i = 0; i < incount; ++i) {
            m_inbuf[c][m_lastread[c] + i] = in[c][i];
        }
        m_lastread[c] += incount;
        
        ippsResamplePolyphase_32f(m_state[c],
                                  m_inbuf[c],
                                  m_lastread[c] - m_history - int(m_time[c]),
                                  m_outbuf[c],
                                  ratio,
                                  0.97f,
                                  &m_time[c],
                                  &outcount);

        v_copy(out[c], m_outbuf[c], outcount);

        ippsMove_32f(m_inbuf[c] + int(m_time[c]) - m_history,
                     m_inbuf[c],
                     m_lastread[c] + m_history - int(m_time[c]));

        m_lastread[c] -= int(m_time[c]) - m_history;
        m_time[c] -= int(m_time[c]) - m_history;

        if (final) {

            // Looks like this actually produces too many samples
            // (additionalcount is a few samples too large).

            // Also, we aren't likely to have enough space in the
            // output buffer as the caller won't have allowed for
            // all the samples we're retrieving here.

            // What to do?

            int additionalcount = 0;

            for (int i = 0; i < m_history; ++i) {
                m_inbuf[c][m_lastread[c] + i] = 0.f;
            }
            
            ippsResamplePolyphase_32f(m_state[c],
                                      m_inbuf[c],
                                      m_lastread[c] - int(m_time[c]),
                                      m_outbuf[c],
                                      ratio,
                                      0.97f,
                                      &m_time[c],
                                      &additionalcount);

            if (m_debugLevel > 2) {
                cerr << "incount = " << incount << ", outcount = " << outcount << ", additionalcount = " << additionalcount << ", sum " << outcount + additionalcount << ", est space = " << lrintf(ceil(incount * ratio)) <<endl;
            }

            v_copy(out[c] + outcount, m_outbuf[c], additionalcount);

            outcount += additionalcount;
        }
    }

    for (int c = 0; c < m_channels; ++c) {
        ippsThreshold_32f_I(out[c], outcount, 1.f, ippCmpGreater);
        ippsThreshold_32f_I(out[c], outcount, -1.f, ippCmpLess);
    }

    return outcount;
}

int
D_IPP::resampleInterleaved(const float *const BQ_R__ in,
                           float *const BQ_R__ out,
                           int incount,
                           double ratio,
                           bool final)
{
    int outcount = 0;

    if (ratio > m_factor) {
        m_factor = ratio;
        m_history = int(m_window * 0.5 * max(1.0, 1.0 / m_factor)) + 1;
    }

    for (int c = 0; c < m_channels; ++c) {
        if (m_lastread[c] + incount + m_history > m_bufsize) {
            setBufSize(m_lastread[c] + incount + m_history);
        }
    }

    for (int c = 0; c < m_channels; ++c) {

        for (int i = 0; i < incount; ++i) {
            m_inbuf[c][m_lastread[c] + i] = in[i * m_channels + c];
        }
        m_lastread[c] += incount;
        
        ippsResamplePolyphase_32f(m_state[c],
                                  m_inbuf[c],
                                  m_lastread[c] - m_history - int(m_time[c]),
                                  m_outbuf[c],
                                  ratio,
                                  0.97f,
                                  &m_time[c],
                                  &outcount);

        ippsMove_32f(m_inbuf[c] + int(m_time[c]) - m_history,
                     m_inbuf[c],
                     m_lastread[c] + m_history - int(m_time[c]));

        m_lastread[c] -= int(m_time[c]) - m_history;
        m_time[c] -= int(m_time[c]) - m_history;
    }

    v_interleave(out, m_outbuf, m_channels, outcount);

    if (final) {

        // Looks like this actually produces too many samples
        // (additionalcount is a few samples too large).

        // Also, we aren't likely to have enough space in the
        // output buffer as the caller won't have allowed for
        // all the samples we're retrieving here.

        // What to do?

        int additionalcount = 0;
        
        for (int c = 0; c < m_channels; ++c) {

            for (int i = 0; i < m_history; ++i) {
                m_inbuf[c][m_lastread[c] + i] = 0.f;
            }
            
            ippsResamplePolyphase_32f(m_state[c],
                                      m_inbuf[c],
                                      m_lastread[c] - int(m_time[c]),
                                      m_outbuf[c],
                                      ratio,
                                      0.97f,
                                      &m_time[c],
                                      &additionalcount);

            if (m_debugLevel > 2) {
                cerr << "incount = " << incount << ", outcount = " << outcount << ", additionalcount = " << additionalcount << ", sum " << outcount + additionalcount << ", est space = " << lrintf(ceil(incount * ratio)) <<endl;
            }
        }

        v_interleave(out + (outcount * m_channels),
                     m_outbuf,
                     m_channels,
                     additionalcount);

        outcount += additionalcount;
    }

    ippsThreshold_32f_I(out, outcount * m_channels, 1.f, ippCmpGreater);
    ippsThreshold_32f_I(out, outcount * m_channels, -1.f, ippCmpLess);

    return outcount;
}

void
D_IPP::reset()
{
    //!!!
}

#endif /* HAVE_IPP */

#ifdef HAVE_LIBSAMPLERATE

class D_SRC : public Resampler::Impl
{
public:
    D_SRC(Resampler::Quality quality, int channels, int maxBufferSize,
          int m_debugLevel);
    ~D_SRC();

    int resample(const float *const BQ_R__ *const BQ_R__ in,
                 float *const BQ_R__ *const BQ_R__ out,
                 int incount,
                 double ratio,
                 bool final);

    int resampleInterleaved(const float *const BQ_R__ in,
                            float *const BQ_R__ out,
                            int incount,
                            double ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    SRC_STATE *m_src;
    float *m_iin;
    float *m_iout;
    double m_lastRatio;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    int m_debugLevel;
};

D_SRC::D_SRC(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_src(0),
    m_iin(0),
    m_iout(0),
    m_lastRatio(1.0),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        cerr << "Resampler::Resampler: using libsamplerate implementation"
                  << endl;
    }

    int err = 0;
    m_src = src_new(quality == Resampler::Best ? SRC_SINC_BEST_QUALITY :
                    quality == Resampler::Fastest ? SRC_LINEAR :
                    SRC_SINC_FASTEST,
                    channels, &err);

    if (err) {
        cerr << "Resampler::Resampler: failed to create libsamplerate resampler: " 
                  << src_strerror(err) << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }

    reset();
}

D_SRC::~D_SRC()
{
    src_delete(m_src);
    deallocate(m_iin);
    deallocate(m_iout);
}

int
D_SRC::resample(const float *const BQ_R__ *const BQ_R__ in,
                float *const BQ_R__ *const BQ_R__ out,
                int incount,
                double ratio,
                bool final)
{
    SRC_DATA data;

    int outcount = (int)lrint(ceil(incount * ratio));

    if (m_channels == 1) {
        data.data_in = const_cast<float *>(*in); //!!!???
        data.data_out = *out;
    } else {
        if (incount * m_channels > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (outcount * m_channels > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data.data_in = m_iin;
        data.data_out = m_iout;
    }

    data.input_frames = incount;
    data.output_frames = outcount;
    data.src_ratio = ratio;
    data.end_of_input = (final ? 1 : 0);

    int err = src_process(m_src, &data);

    if (err) {
        cerr << "Resampler::process: libsamplerate error: "
                  << src_strerror(err) << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, (int)data.output_frames_gen);
    }

    m_lastRatio = ratio;

    return (int)data.output_frames_gen;
}

int
D_SRC::resampleInterleaved(const float *const BQ_R__ in,
                           float *const BQ_R__ out,
                           int incount,
                           double ratio,
                           bool final)
{
    SRC_DATA data;

    int outcount = (int)lrint(ceil(incount * ratio));

    data.data_in = const_cast<float *>(in);
    data.data_out = out;

    data.input_frames = incount;
    data.output_frames = outcount;
    data.src_ratio = ratio;
    data.end_of_input = (final ? 1 : 0);

    int err = src_process(m_src, &data);

    if (err) {
        cerr << "Resampler::process: libsamplerate error: "
                  << src_strerror(err) << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    m_lastRatio = ratio;

    return (int)data.output_frames_gen;
}

void
D_SRC::reset()
{
    src_reset(m_src);
}

#endif /* HAVE_LIBSAMPLERATE */

#ifdef HAVE_LIBRESAMPLE

class D_Resample : public Resampler::Impl
{
public:
    D_Resample(Resampler::Quality quality, int channels, int maxBufferSize,
          int m_debugLevel);
    ~D_Resample();

    int resample(const float *const BQ_R__ *const BQ_R__ in,
                 float *const BQ_R__ *const BQ_R__ out,
                 int incount,
                 double ratio,
                 bool final);

    int resampleInterleaved(const float *const BQ_R__ in,
                            float *const BQ_R__ out,
                            int incount,
                            double ratio,
                            bool final);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    void *m_src;
    float *m_iin;
    float *m_iout;
    float m_lastRatio;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    int m_debugLevel;
};

D_Resample::D_Resample(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_src(0),
    m_iin(0),
    m_iout(0),
    m_lastRatio(1.f),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        cerr << "Resampler::Resampler: using libresample implementation"
                  << endl;
    }

    float min_factor = 0.125f;
    float max_factor = 8.0f;

    m_src = resample_open(quality == Resampler::Best ? 1 : 0, min_factor, max_factor);

    if (!m_src) {
        cerr << "Resampler::Resampler: failed to create libresample resampler: " 
                  << endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }

    reset();
}

D_Resample::~D_Resample()
{
    resample_close(m_src);
    if (m_iinsize > 0) {
        deallocate(m_iin);
    }
    if (m_ioutsize > 0) {
        deallocate(m_iout);
    }
}

int
D_Resample::resample(const float *const BQ_R__ *const BQ_R__ in,
                     float *const BQ_R__ *const BQ_R__ out,
                     int incount,
                     double ratio,
                     bool final)
{
    float *data_in;
    float *data_out;
    int input_frames, output_frames, end_of_input, source_used;
    float src_ratio;

    int outcount = (int)lrint(ceil(incount * ratio));

    if (m_channels == 1) {
        data_in = const_cast<float *>(*in); //!!!???
        data_out = *out;
    } else {
        if (incount * m_channels > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (outcount * m_channels > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data_in = m_iin;
        data_out = m_iout;
    }

    input_frames = incount;
    output_frames = outcount;
    src_ratio = ratio;
    end_of_input = (final ? 1 : 0);

    int output_frames_gen = resample_process(m_src,
                                             src_ratio,
                                             data_in,
                                             input_frames,
                                             end_of_input,
                                             &source_used,
                                             data_out,
                                             output_frames);

    if (output_frames_gen < 0) {
        cerr << "Resampler::process: libresample error: "
                  << endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, output_frames_gen);
    }

    m_lastRatio = ratio;

    return output_frames_gen;
}

int
D_Resample::resampleInterleaved(const float *const BQ_R__ in,
                                float *const BQ_R__ out,
                                int incount,
                                double ratio,
                                bool final)
{
    int input_frames, output_frames, end_of_input, source_used;
    float src_ratio;

    int outcount = (int)lrint(ceil(incount * ratio));

    input_frames = incount;
    output_frames = outcount;
    src_ratio = ratio;
    end_of_input = (final ? 1 : 0);

    int output_frames_gen = resample_process(m_src,
                                             src_ratio,
                                             const_cast<float *>(in),
                                             input_frames,
                                             end_of_input,
                                             &source_used,
                                             out,
                                             output_frames);

    if (output_frames_gen < 0) {
        cerr << "Resampler::process: libresample error: "
                  << endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    m_lastRatio = ratio;

    return output_frames_gen;
}

void
D_Resample::reset()
{
}

#endif /* HAVE_LIBRESAMPLE */

#ifdef USE_SPEEX
    
class D_Speex : public Resampler::Impl
{
public:
    D_Speex(Resampler::Quality quality, int channels, int maxBufferSize,
            int debugLevel);
    ~D_Speex();

    int resample(const float *const BQ_R__ *const BQ_R__ in,
                 float *const BQ_R__ *const BQ_R__ out,
                 int incount,
                 double ratio,
                 bool final);

    int resampleInterleaved(const float *const BQ_R__ in,
                            float *const BQ_R__ out,
                            int incount,
                            double ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    SpeexResamplerState *m_resampler;
    float *m_iin;
    float *m_iout;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    float m_lastratio;
    bool m_initial;
    int m_debugLevel;

    void setRatio(double);
};

D_Speex::D_Speex(Resampler::Quality quality, int channels, int maxBufferSize,
                 int debugLevel) :
    m_resampler(0),
    m_iin(0),
    m_iout(0),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_lastratio(1),
    m_initial(true),
    m_debugLevel(debugLevel)
{
    int q = (quality == Resampler::Best ? 10 :
             quality == Resampler::Fastest ? 0 : 4);

    if (m_debugLevel > 0) {
        cerr << "Resampler::Resampler: using Speex implementation with q = "
                  << q 
                  << endl;
    }

    int err = 0;
    m_resampler = speex_resampler_init_frac(m_channels,
                                            1, 1,
                                            48000, 48000, // irrelevant
                                            q,
                                            &err);
    

    if (err) {
        cerr << "Resampler::Resampler: failed to create Speex resampler" 
                  << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }
}

D_Speex::~D_Speex()
{
    speex_resampler_destroy(m_resampler);
    deallocate<float>(m_iin);
    deallocate<float>(m_iout);
}

void
D_Speex::setRatio(double ratio)
{
    // Speex wants a ratio of two unsigned integers, not a single
    // float.  Let's do that.

    unsigned int big = 272408136U; 
    unsigned int denom = 1, num = 1;

    if (ratio < 1.f) {
        denom = big;
        double dnum = double(big) * double(ratio);
        num = (unsigned int)dnum;
    } else if (ratio > 1.f) {
        num = big;
        double ddenom = double(big) / double(ratio);
        denom = (unsigned int)ddenom;
    }
    
    if (m_debugLevel > 1) {
        cerr << "D_Speex: Desired ratio " << ratio << ", requesting ratio "
                  << num << "/" << denom << " = " << float(double(num)/double(denom))
                  << endl;
    }
    
    int err = speex_resampler_set_rate_frac
        (m_resampler, denom, num, 48000, 48000);

    if (err) {
        cerr << "Resampler::Resampler: failed to set rate on Speex resampler" 
                  << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }
    
    speex_resampler_get_ratio(m_resampler, &denom, &num);
    
    if (m_debugLevel > 1) {
        cerr << "D_Speex: Desired ratio " << ratio << ", got ratio "
                  << num << "/" << denom << " = " << float(double(num)/double(denom))
                  << endl;
    }
    
    m_lastratio = ratio;

    if (m_initial) {
        speex_resampler_skip_zeros(m_resampler);
        m_initial = false;
    }
}

int
D_Speex::resample(const float *const BQ_R__ *const BQ_R__ in,
                  float *const BQ_R__ *const BQ_R__ out,
                  int incount,
                  double ratio,
                  bool final)
{
    if (ratio != m_lastratio) {
        setRatio(ratio);
    }

    unsigned int uincount = incount;
    unsigned int outcount = (unsigned int)lrint(ceil(incount * ratio)); //!!! inexact now

    float *data_in, *data_out;

    if (m_channels == 1) {
        data_in = const_cast<float *>(*in);
        data_out = *out;
    } else {
        if (int(incount * m_channels) > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (int(outcount * m_channels) > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data_in = m_iin;
        data_out = m_iout;
    }

    int err = speex_resampler_process_interleaved_float(m_resampler,
                                                        data_in,
                                                        &uincount,
                                                        data_out,
                                                        &outcount);

    if (err) {
        cerr << "Resampler::Resampler: Speex resampler returned error "
                  << err
                  << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, outcount);
    }

    return outcount;
}

int
D_Speex::resampleInterleaved(const float *const BQ_R__ in,
                             float *const BQ_R__ out,
                             int incount,
                             double ratio,
                             bool final)
{
    if (ratio != m_lastratio) {
        setRatio(ratio);
    }

    unsigned int uincount = incount;
    unsigned int outcount = (unsigned int)lrint(ceil(incount * ratio)); //!!! inexact now

    float *data_in = const_cast<float *>(in);
    float *data_out = out;

    int err = speex_resampler_process_interleaved_float(m_resampler,
                                                        data_in,
                                                        &uincount,
                                                        data_out,
                                                        &outcount);

    if (err) {
        cerr << "Resampler::Resampler: Speex resampler returned error "
                  << err
                  << endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    return outcount;
}

void
D_Speex::reset()
{
    speex_resampler_reset_mem(m_resampler);
}

#endif

} /* end namespace Resamplers */

Resampler::Resampler(Resampler::Quality quality, int channels,
                     int maxBufferSize, int debugLevel)
{
    m_method = -1;
    
    switch (quality) {

    case Resampler::Best:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
        break;

    case Resampler::FastestTolerable:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
        break;

    case Resampler::Fastest:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
        break;
    }

    if (m_method == -1) {
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << endl;
        abort();
    }

    switch (m_method) {
    case 0:
#ifdef HAVE_IPP
        d = new Resamplers::D_IPP(quality, channels, maxBufferSize, debugLevel);
#else
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << endl;
        abort();
#endif
        break;

    case 1:
#ifdef HAVE_LIBSAMPLERATE
        d = new Resamplers::D_SRC(quality, channels, maxBufferSize, debugLevel);
#else
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << endl;
        abort();
#endif
        break;

    case 2:
#ifdef USE_SPEEX
        d = new Resamplers::D_Speex(quality, channels, maxBufferSize, debugLevel);
#else
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << endl;
        abort();
#endif
        break;

    case 3:
#ifdef HAVE_LIBRESAMPLE
        d = new Resamplers::D_Resample(quality, channels, maxBufferSize, debugLevel);
#else
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << endl;
        abort();
#endif
        break;
    }

    if (!d) {
        cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize
                  << "): Internal error: No implementation selected"
                  << endl;
        abort();
    }
}

Resampler::~Resampler()
{
    delete d;
}

int 
Resampler::resample(const float *const BQ_R__ *const BQ_R__ in,
                    float *const BQ_R__ *const BQ_R__ out,
                    int incount, double ratio, bool final)
{
    return d->resample(in, out, incount, ratio, final);
}

int 
Resampler::resampleInterleaved(const float *const BQ_R__ in,
                               float *const BQ_R__ out,
                               int incount, double ratio, bool final)
{
    return d->resampleInterleaved(in, out, incount, ratio, final);
}

int
Resampler::getChannelCount() const
{
    return d->getChannelCount();
}

void
Resampler::reset()
{
    d->reset();
}

}
