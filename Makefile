
# Add to RESAMPLE_DEFINES the relevant options for your desired
# third-party library support.
#
# Available options are
#
#  -DHAVE_IPP            Intel's Integrated Performance Primitives are available
#  -DHAVE_LIBRESAMPLE    The libresample library is available
#  -DHAVE_LIBSAMPLERATE  The libsamplerate library is available
#  -DUSE_SPEEX           Compile the built-in Speex-derived resampler
#
# You may define more than one of these. If you define USE_SPEEX, the
# code will be compiled in and will be used when it is judged to be
# the best available option for a given quality setting. The default,
# if no flags are supplied, is for the code to refuse to compile.
# 
# Add any relevant -I flags for include paths as well.
#
# Note that you must supply the same flags when including bqresample
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqresample source files to your
# application's build system and not build a bqresample library at all.)

RESAMPLE_DEFINES	:= -DUSE_SPEEX


# Add to VECTOR_DEFINES the relevant options for your desired
# third-party library support.
#
# Available options are
#
#  -DHAVE_IPP    Intel's Integrated Performance Primitives are available
#  -DHAVE_VDSP   Apple's Accelerate framework is available
#
# These are optional (they affect performance, not function) and you
# may define more than one of them.
# 
# Add any relevant -I flags for include paths as well.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

VECTOR_DEFINES		:= -I../bqvec


# Add to ALLOCATOR_DEFINES options relating to aligned malloc.
#
# Available options are
#
#  -DHAVE_POSIX_MEMALIGN     The posix_memalign call is available in sys/mman.h
#  -DLACK_POSIX_MEMALIGN     The posix_memalign call is not available
#
#  -DMALLOC_IS_ALIGNED       The malloc call already returns aligned memory
#  -DMALLOC_IS_NOT_ALIGNED   The malloc call does not return aligned memory
#
#  -DUSE_OWN_ALIGNED_MALLOC  No aligned malloc is available, roll your own
#
#  -DLACK_BAD_ALLOC          The C++ library lacks the std::bad_alloc exception
#
# Here "aligned" is assumed to mean "aligned enough for whatever
# vector stuff the space will be used for" which most likely means
# 16-byte alignment.
#
# The default is to use _aligned_malloc when building with Visual C++,
# system malloc when building on OS/X, and posix_memalign otherwise.
#
# Note that you must supply the same flags when including bqvec
# headers later as you are using now when compiling the library. (You
# may find it simplest to just add the bqvec source files to your
# application's build system and not build a bqvec library at all.)

ALLOCATOR_DEFINES 	:= 


SRC_DIR	:= src
SPEEX_DIR := speex
HEADER_DIR := bqresample

SOURCES	:= $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SPEEX_DIR)/*.c)
HEADERS	:= $(wildcard $(HEADER_DIR)/*.h) $(wildcard $(SRC_DIR)/*.h)

OBJECTS	:= $(SOURCES:.cpp=.o)
OBJECTS	:= $(OBJECTS:.c=.o)

CXXFLAGS := $(RESAMPLE_DEFINES) $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) -I. -Wall -Werror
CFLAGS	:= $(RESAMPLE_DEFINES) $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) -I. -Wall -Werror

LIBRARY	:= libbqresample.a

all:	$(LIBRARY)

$(LIBRARY):	$(OBJECTS)
	$(AR) rc $@ $^

clean:		
	rm -f $(OBJECTS)

distclean:	clean
	rm -f $(LIBRARY)

depend:
	makedepend -Y -fMakefile $(SOURCES) $(HEADERS)


# DO NOT DELETE
