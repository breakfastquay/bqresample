/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    bqresample

    A small library wrapping various audio sample rate conversion
    implementations in C++.

    Copyright 2007-2023 Particular Programs Ltd.

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

#include <cmath>

namespace breakfastquay {

void pickNearestRational(double ratio, int max_denom, int &num, int &denom)
{
    // Farey algorithm, see
    // https://www.johndcook.com/blog/2010/10/20/best-rational-approximation/
    double a = 0.0, b = 1.0, c = 1.0, d = 0.0;
    double pa = a, pb = b, pc = c, pd = d;
    double eps = 1e-9;
    while (b <= max_denom && d <= max_denom) {
        double mediant = (a + c) / (b + d);
        if (fabs(ratio - mediant) < eps) {
            if (b + d <= max_denom) {
                num = a + c;
                denom = b + d;
                return;
            } else if (d > b) {
                num = c;
                denom = d;
                return;
            } else {
                num = a;
                denom = b;
                return;
            }
        }
        if (ratio > mediant) {
            pa = a; pb = b;
            a += c; b += d;
        } else {
            pc = c; pd = d;
            c += a; d += b;
        }
    }
    if (fabs(ratio - (pc / pd)) < fabs(ratio - (pa / pb))) {
        num = pc;
        denom = pd;
    } else {
        num = pa;
        denom = pb;
    }
}

}

