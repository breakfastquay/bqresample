
bqresample
==========

A small C++ library wrapping various audio sample rate conversion
libraries.

Requires the bqvec library.

This code originated as part of the Rubber Band Library written by the
same authors (see https://hg.sr.ht/~breakfastquay/rubberband/).
It has been pulled out into a separate library and relicensed under a
more permissive licence.

C++ standard required: C++98 (does not use C++11 or newer features)

 * To compile: read and follow the notes in Makefile, edit the Makefile,
   then make test. Or else use one of the pre-edited Makefiles in the
   build directory.

 * Depends on: [bqvec](https://hg.sr.ht/~breakfastquay/bqvec)

 * See also: [bqfft](https://hg.sr.ht/~breakfastquay/bqfft) [bqaudioio](https://hg.sr.ht/~breakfastquay/bqaudioio) [bqthingfactory](https://hg.sr.ht/~breakfastquay/bqthingfactory) [bqaudiostream](https://hg.sr.ht/~breakfastquay/bqaudiostream)

[![Build status](https://builds.sr.ht/~breakfastquay/bqresample.svg)](https://builds.sr.ht/~breakfastquay/bqresample?)

Copyright 2007-2020 Particular Programs Ltd.
Uses Speex code, see speex/COPYING for copyright and licence information.

