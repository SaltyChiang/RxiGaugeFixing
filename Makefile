CC = clang
CXX = clang++
CFLAGS = -std=c99 -O3 -Wall -g -fopenmp
CXXFLAGS = -std=c++14 -O3 -Wall -g -fopenmp
LDFLAGS =
LIBS =

all : bin/RxiGF.exe

SRC = $(wildcard src/*.cc)
OBJ = $(patsubst %.cc, bin/%.o, $(notdir $(SRC)))
INC = include/RxiGF_lattice.h include/RxiGF_macro.h

$(OBJ) : bin/%.o : src/%.cc include/%.h $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

bin/RxiGF.exe : RxiGF_main.cc $(OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	rm -rf bin/*