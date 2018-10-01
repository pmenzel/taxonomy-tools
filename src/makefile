CC = gcc
CXX = g++
CFLAGS = -O3 -DNDEBUG
CXXFLAGS = -O3 -pthread -std=c++11 -DNDEBUG
LDLIBS = -lpthread -lz

ifeq ($(uname -s), "Darwin")
LD_LIBS_STATIC = -Wl,-all_load -lpthread -lz -Wl,-noall_load
else
LD_LIBS_STATIC = -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive
endif

all: makefile blasthits2krona
	mkdir -p ../bin && cp blasthits2krona ../bin/

blasthits2krona: makefile blasthits2krona.o util.o
	$(CXX) $(LDFLAGS) -o blasthits2krona blasthits2krona.o util.o

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f -v blasthits2krona ../bin/*
	find . -name "*.o" -delete

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer
debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized
debug: all

.PHONY: clean debug static
