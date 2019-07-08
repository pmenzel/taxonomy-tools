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

all: makefile blasthits2krona lca acc2name subtree
	mkdir -p ../bin && cp blasthits2krona lca acc2name subtree ../bin/

acc2name: makefile acc2name.o util.o
	$(CXX) $(LDFLAGS) -o acc2name acc2name.o util.o

blasthits2krona: makefile blasthits2krona.o util.o
	$(CXX) $(LDFLAGS) -o blasthits2krona blasthits2krona.o util.o

lca: makefile lca.o util.o
	$(CXX) $(LDFLAGS) -o lca lca.o util.o

subtree: makefile subtree.o util.o
	$(CXX) $(LDFLAGS) -o subtree subtree.o util.o

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	rm -f -v blasthits2krona lca acc2name subtree ../bin/*
	find . -name "*.o" -delete

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer
debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized
debug: all

.PHONY: clean debug static test
