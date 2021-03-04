
#include "bqresample/Resampler.h"

#include <sndfile.h>

#include <iostream>

#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

void usage()
{
    cerr << "Usage: resample -to <rate> <infile> <outfile>" << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    if (argc != 5 || strcmp(argv[1], "-to")) {
        usage();
    }

    double target = strtod(argv[2], 0);
    if (!target) {
        usage();
    }

    string infilename = argv[3];
    string outfilename = argv[4];
    
    SF_INFO info_in;
    SNDFILE *file_in = sf_open(infilename.c_str(), SFM_READ, &info_in);
    if (!file_in) {
        cerr << "failed to open " << infilename << endl;
        return 1;
    }

    if (info_in.channels > 1) { //!!! todo
        cerr << "mono required" << endl;
        return 2;
    }

    int output_rate = round(target);
    
    cerr << "input rate = " << info_in.samplerate << endl;
    cerr << "output rate = " << output_rate << endl;

    double ratio = target / info_in.samplerate;
    cerr << "ratio = " << ratio << endl;
        
    SF_INFO info_out;
    memset(&info_out, 0, sizeof(SF_INFO));
    info_out.channels = 1;
    info_out.format = info_in.format;
    info_out.samplerate = output_rate;
    SNDFILE *file_out = sf_open(outfilename.c_str(), SFM_WRITE, &info_out);
    if (!file_out) {
        cerr << "failed to open " << outfilename << endl;
        return 1;
    }

    int ibs = 1024;
    int obs = ceil(ibs * ratio * 10.0);
    float *ibuf = new float[ibs];
    float *obuf = new float[obs];

    breakfastquay::Resampler::Parameters parameters;
    parameters.quality = breakfastquay::Resampler::Best;
    parameters.dynamism = breakfastquay::Resampler::RatioOftenChanging;
    parameters.initialSampleRate = info_in.samplerate;
    parameters.debugLevel = 1;
    breakfastquay::Resampler resampler(parameters, 1);

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
        ++n;
    }

    cerr << endl;
    
    sf_close(file_in);
    sf_close(file_out);
}
