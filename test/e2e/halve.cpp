
#include "bqresample/Resampler.h"

#include <sndfile.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

static const string programName = "halve";
static const double initialRatio = 0.5;
static double nextRatio(double ratio) {
    return ratio;
}
static bool isRatioChanging() {
    return false;
}

#include "e2e.cpp"
