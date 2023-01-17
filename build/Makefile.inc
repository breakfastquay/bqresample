
SRC_DIR	:= src
TEST_DIR := test
EXAMPLE_DIR := example
SPEEX_DIR := ext/speex
HEADER_DIR := bqresample

SOURCES	:= $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(SPEEX_DIR)/*.c)
HEADERS	:= $(wildcard $(HEADER_DIR)/*.h) $(wildcard $(SRC_DIR)/*.h)

OBJECTS	:= $(SOURCES:.cpp=.o)
OBJECTS	:= $(OBJECTS:.c=.o)

TEST_SOURCES	:= $(wildcard $(TEST_DIR)/*.cpp) $(wildcard $(TEST_DIR)/e2e/*.cpp)
TEST_OBJECTS	:= $(TEST_SOURCES:.cpp=.o)

EXAMPLE_SOURCES	:= $(wildcard $(EXAMPLE_DIR)/*.cpp)
EXAMPLE_OBJECTS	:= $(EXAMPLE_SOURCES:.cpp=.o)

OPTFLAGS ?= -O3 -ffast-math

CXXFLAGS ?= -std=c++98 -fpic -Wall -Wextra -Werror $(RESAMPLE_DEFINES) $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) $(OPTFLAGS) -I. $(THIRD_PARTY_INCLUDES) -I../bqvec

CFLAGS	:= -fpic -Wall -Wextra $(RESAMPLE_DEFINES) $(VECTOR_DEFINES) $(ALLOCATOR_DEFINES) $(OPTFLAGS) -I. $(THIRD_PARTY_INCLUDES)

LIBRARY	:= libbqresample.a

all:	$(LIBRARY)

test:	$(LIBRARY) test-resampler test/e2e/halve test/e2e/oversample test/e2e/random test/e2e/undulating
	./test-resampler
	mkdir -p ./test/e2e/out
	./test/e2e/halve ./test/e2e/sweep-log.wav ./test/e2e/out/halve.wav
	./test/e2e/oversample ./test/e2e/sweep-log.wav ./test/e2e/out/oversample.wav
	./test/e2e/random ./test/e2e/sweep-log.wav ./test/e2e/out/random.wav
	./test/e2e/undulating ./test/e2e/sweep-log.wav ./test/e2e/out/undulating.wav

example:	$(LIBRARY) example/resample

valgrind:	$(LIBRARY) test-resampler
	valgrind ./test-resampler

$(LIBRARY):	$(OBJECTS)
	$(AR) rc $@ $^

test-resampler:	test/TestResampler.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBRARY) -lboost_unit_test_framework -L../bqvec -lbqvec $(THIRD_PARTY_LIBS)

example/resample: example/resample.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBRARY) -lsndfile -L../bqvec -lbqvec $(THIRD_PARTY_LIBS)

test/e2e/%: test/e2e/%.o $(LIBRARY)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBRARY) -lsndfile -L../bqvec -lbqvec $(THIRD_PARTY_LIBS)

clean:		
	rm -f $(OBJECTS) $(TEST_OBJECTS) $(EXAMPLE_OBJECTS)

distclean:	clean
	rm -f $(LIBRARY) test-resampler test/e2e/halve test/e2e/oversample test/e2e/oversample test/e2e/random test/e2e/undulating

depend:
	makedepend $(RESAMPLE_DEFINES) -Y -fbuild/Makefile.inc $(SOURCES) $(HEADERS) $(TEST_SOURCES) $(EXAMPLE_SOURCES)


# DO NOT DELETE

src/BQResampler.o: src/BQResampler.h
src/Resampler.o: bqresample/Resampler.h src/BQResampler.h
test/TestResampler.o: bqresample/Resampler.h
test/e2e/halve.o: bqresample/Resampler.h test/e2e/e2e.cpp
test/e2e/oversample.o: bqresample/Resampler.h test/e2e/e2e.cpp
test/e2e/random.o: bqresample/Resampler.h test/e2e/e2e.cpp
test/e2e/undulating.o: bqresample/Resampler.h test/e2e/e2e.cpp
example/resample.o: bqresample/Resampler.h
