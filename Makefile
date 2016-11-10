CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -Wno-unused-function -fPIC -pedantic #-fpermissive
CXXFLAGS += -I./
BIN=main

SRC=$(wildcard *.cpp)
OBJ=$(SRC:%.cpp=%.o)

%.o: %.c
	$(CXX) $@ -c $(CXXFLAGS) $<

all: $(OBJ)
	$(CXX) -o $(BIN) $^

clean:
	rm -f *.o
	rm -f $(BIN)
