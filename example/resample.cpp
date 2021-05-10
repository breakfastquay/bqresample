
#include "bqresample/Resampler.h"

#include <sndfile.h>

#include <iostream>

#include <cstring>
#include <cstdlib>
#include <cmath>

using namespace std;

void usage()
{
    cerr << "Usage: resample [-v] [-c <converter>] -to <rate> <infile> <outfile>" << endl;
    cerr << "where <converter> may be 0, 1, or 2, for best, medium, or fastest respectively." << endl;
    cerr << "Supply -v for verbose output." << endl;
    exit(2);
}

int main(int argc, char **argv)
{
    double target = 0.0;
    int quality = 0;
    int arg;
    bool verbose = false;

    for (arg = 1; arg + 2 < argc; ++arg) {
        if (!strcmp(argv[arg], "-c")) {
            char *e = argv[arg+1];
            quality = strtol(argv[arg+1], &e, 10);
            if (*e || (quality < 0) || (quality > 2)) {
                cerr << "error: invalid converter \""
                     << argv[arg+1] << "\" (must be 0, 1, 2)" << endl;
                usage();
            }
            ++arg;
            continue;
        } else if (!strcmp(argv[arg], "-to")) {
            target = strtod(argv[arg+1], 0);
            if (!target) {
                cerr << "error: invalid target \"" << argv[arg+1]
                     << "\" (must be numeric)" << endl;
                usage();
            }
            ++arg;
            continue;
        } else if (!strcmp(argv[arg], "-v")) {
            verbose = true;
        } else {
            cerr << "error: unexpected option \"" << argv[arg] << "\"" << endl;
            usage();
        }
    }

    if (!target || arg + 2 != argc) {
        usage();
    }

    string infilename = argv[arg];
    string outfilename = argv[arg+1];
    
    SF_INFO info_in;
    SNDFILE *file_in = sf_open(infilename.c_str(), SFM_READ, &info_in);
    if (!file_in) {
        cerr << "failed to open " << infilename << endl;
        return 1;
    }

    int output_rate = round(target);
    int channels = info_in.channels;
    
    cerr << "input rate = " << info_in.samplerate << endl;
    cerr << "output rate = " << output_rate << endl;

    double ratio = target / info_in.samplerate;
    cerr << "ratio = " << ratio << endl;

    breakfastquay::Resampler::Parameters parameters;
    switch (quality)  {
    case 0:
        parameters.quality = breakfastquay::Resampler::Best;
        cerr << "quality = best" << endl;
        break;
    case 1:
        parameters.quality = breakfastquay::Resampler::FastestTolerable;
        cerr << "quality = middling" << endl;
        break;
    case 2:
        parameters.quality = breakfastquay::Resampler::Fastest;
        cerr << "quality = worst" << endl;
        break;
    }        

    SF_INFO info_out;
    memset(&info_out, 0, sizeof(SF_INFO));
    info_out.channels = channels;
    info_out.format = info_in.format;
    info_out.samplerate = output_rate;
    SNDFILE *file_out = sf_open(outfilename.c_str(), SFM_WRITE, &info_out);
    if (!file_out) {
        cerr << "failed to open " << outfilename << endl;
        return 1;
    }

    int ibs = 1024;
    int obs = ceil(ibs * ratio);
    float *ibuf = new float[ibs * channels];
    float *obuf = new float[obs * channels];

    parameters.dynamism = breakfastquay::Resampler::RatioMostlyFixed;
    parameters.ratioChange = breakfastquay::Resampler::SuddenRatioChange;
    parameters.initialSampleRate = info_in.samplerate;
    parameters.debugLevel = (verbose ? 1 : 0);
    breakfastquay::Resampler resampler(parameters, info_in.channels);

    int n = 0;
    while (true) {
        int count = sf_readf_float(file_in, ibuf, ibs);
        if (verbose) {
            cerr << ".";
        }
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
            for (int i = 0; i < got * channels; ++i) {
                if (obuf[i] < -1.0) obuf[i] = -1.0;
                if (obuf[i] > 1.0) obuf[i] = 1.0;
            }
            sf_writef_float(file_out, obuf, got);
        }
        ++n;
    }

    if (verbose) {
        cerr << endl;
    }
    
    sf_close(file_in);
    sf_close(file_out);
    
    delete[] ibuf;
    delete[] obuf;
}
