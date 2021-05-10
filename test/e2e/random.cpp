
#include "bqresample/Resampler.h"

#include <sndfile.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

void usage()
{
    cerr << "This is a test program for bqresample. Please do not try to use it in earnest." << endl;
    cerr << "Usage: random <infile> <outfile>" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    if (argc != 3) {
        usage();
    }

    string infilename = argv[1];
    string outfilename = argv[2];
    
    SF_INFO info_in;
    SNDFILE *file_in = sf_open(infilename.c_str(), SFM_READ, &info_in);
    if (!file_in) {
        cerr << "failed to open " << infilename << endl;
        return 1;
    }

    int channels = info_in.channels;
    
//    cerr << "input rate = " << info_in.samplerate << endl;

    double ratio = 0.5; // initially
        
    SF_INFO info_out;
    memset(&info_out, 0, sizeof(SF_INFO));
    info_out.channels = channels;
    info_out.format = info_in.format;
    info_out.samplerate = info_in.samplerate;
    SNDFILE *file_out = sf_open(outfilename.c_str(), SFM_WRITE, &info_out);
    if (!file_out) {
        cerr << "failed to open " << outfilename << endl;
        return 1;
    }

    int ibs = 1024;
    int obs = ceil(ibs * ratio * 10);
    float *ibuf = new float[ibs];
    float *obuf = new float[obs];

    breakfastquay::Resampler::Parameters parameters;
    parameters.quality = breakfastquay::Resampler::Best;
    parameters.dynamism = breakfastquay::Resampler::RatioOftenChanging;
    parameters.ratioChange = breakfastquay::Resampler::SmoothRatioChange;
    parameters.initialSampleRate = info_in.samplerate;
//    parameters.debugLevel = 1;
    breakfastquay::Resampler resampler(parameters, channels);

    int n = 0;
    while (true) {
        int count = sf_readf_float(file_in, ibuf, ibs);
        cerr << ".";
        if (count < 0) {
            cerr << "error: count = " << count << endl;
            break;
        }
        
        bool final = (count < ibs);
        int got = resampler.resampleInterleaved
            (obuf, obs, ibuf, count, ratio, final);
        
        if (got == 0 && final) {
            break;
        } else {
            for (int i = 0; i < got; ++i) {
                if (obuf[i] < -1.0) obuf[i] = -1.0;
                if (obuf[i] > 1.0) obuf[i] = 1.0;
            }
            sf_writef_float(file_out, obuf, got);
        }
        
        ratio = drand48() * 2.5 + 0.1;
        ++n;
    }

    cerr << endl;
    
    sf_close(file_in);
    sf_close(file_out);
    
    delete[] ibuf;
    delete[] obuf;
}
