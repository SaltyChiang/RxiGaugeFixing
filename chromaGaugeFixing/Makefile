# CC = clang
# CXX = clang++
# CFLAGS = -std=c99 -O3 -Wall -g -fopenmp
# CXXFLAGS = -std=c++11 -O3 -Wall -g -fopenmp
# LDFLAGS =
# LIBS =
CHROMA = /home/irisk/usqcd/chroma-install

QDP_CONFIG = $(CHROMA)/bin/qdp++-config
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
