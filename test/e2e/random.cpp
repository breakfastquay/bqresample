
#include "bqresample/Resampler.h"

#include <sndfile.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

static const string programName = "random";
static const double initialRatio = 0.5;
static double nextRatio(double) {
    return drand48() * 2.5 + 0.1;
}
static bool isRatioChanging() {
    return true;
}

#include "e2e.cpp"
