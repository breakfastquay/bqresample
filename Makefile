
# Add to RESAMPLE_DEFINES the relevant options for your desired
# implementation and/or third-party library support.
#
# Available library options are
#
#  -DHAVE_LIBSAMPLERATE  The libsamplerate library is available (recommended)
#  -DHAVE_LIBSPEEXDSP    The speexdsp library is available (recommended)
#  -DHAVE_IPP            Intel's Integrated Performance Primitives are available
#  -DUSE_BQRESAMPLER     Compile the built-in BQ resampler (pretty good)
#  -DUSE_SPEEX           Compile the bundled Speex-derived resampler
#
# You may define more than one of these, and the implementation used
# will depend on the quality setting you request - but it is usually
# better to stick with a single known library. If no flags are
# supplied, the code will refuse to compile.
#
RESAMPLE_DEFINES	:= -DUSE_BQRESAMPLER


# Add to VECTOR_DEFINES and ALLOCATOR_DEFINES any options desired for
# the bqvec library (that are not already defined in RESAMPLE_DEFINES).
# See the bqvec build documentation for more details.
#
VECTOR_DEFINES 		:= 
ALLOCATOR_DEFINES 	:= 


# Add any related includes and libraries here (e.g. -lsamplerate in
# THIRD_PARTY_LIBS if using libsamplerate)
#
THIRD_PARTY_INCLUDES	:=
THIRD_PARTY_LIBS	:=


# If you are including a set of bq libraries into a project, you can
# override variables for all of them (including all of the above) in
# the following file, which all bq* Makefiles will include if found

-include ../Makefile.inc-bq


# This project-local Makefile describes the source files and contains
# no routinely user-modifiable parts

include build/Makefile.inc
