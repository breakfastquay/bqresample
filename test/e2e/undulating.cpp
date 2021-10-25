
#include "bqresample/Resampler.h"

#include <sndfile.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

static const string programName = "undulating";
static const double initialRatio = 0.5;
static double nextRatio(double ratio) {
    static double factor = 1.01;
    ratio *= factor;
    if (ratio > 3.0) factor = 0.99;
    else if (ratio < 0.3) factor = 1.01;
    return ratio;
}
static bool isRatioChanging() {
    return true;
}

#include "e2e.cpp"
