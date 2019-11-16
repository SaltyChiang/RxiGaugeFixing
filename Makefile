CC = clang
CXX = clang++
CFLAGS = -std=c99 -O3 -Wall -g -fopenmp
CXXFLAGS = -std=c++11 -O3 -Wall -g -fopenmp
LDFLAGS =
LIBS =

all : bin/zetaGF.exe

SRC = $(wildcard src/*.cc)
OBJ = $(patsubst %.cc, bin/%.o, $(notdir $(SRC)))
INC = include/zetaGF_lattice.h include/zetaGF_time.h

$(OBJ) : bin/%.o : src/%.cc include/%.h $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

bin/zetaGF.exe : zetaGF_main.cc $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@
