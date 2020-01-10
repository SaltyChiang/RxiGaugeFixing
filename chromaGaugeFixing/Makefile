# CXX = mpicxx
# CXXFLAGS = -fopenmp -O3 -std=c++11 -Drestrict=__restrict__ -march=native -I$(CHROMA)/include -I/usr/include/libxml2
# LDFLAGS = -L$(CHROMA)/lib
# LIBS = -lqdp -lXPathReader -lxmlWriter -lqio -llime -lxml2 -lqmp -lqmp -lfiledb -lfilehash
QDPXX = $(CHROMA)

QDP_CONFIG = $(QDPXX)/bin/qdp++-config
CXX = $(shell $(QDP_CONFIG) --cxx)
CXXFLAGS = $(shell $(QDP_CONFIG) --cxxflags)
LDFLAGS = $(shell $(QDP_CONFIG) --ldflags)
LIBS = $(shell $(QDP_CONFIG) --libs)

all : bin
	+ make bin/chromagfix.exe

SRC = $(wildcard src/*.cc)
OBJ = $(patsubst %.cc, bin/%.o, $(notdir $(SRC)))
INC = include/chromabase.h include/helpfunc.h

bin :
	mkdir $@

$(OBJ) : bin/%.o : src/%.cc include/%.h $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

bin/chromagfix.exe : chromaRxi.cc $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

clean :
	rm -rf bin/*
