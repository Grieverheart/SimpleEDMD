SRC=$(wildcard src/*.cpp)
SHAPE_SRC=$(wildcard src/shape/*.cpp)
OBJ=$(patsubst src/%.cpp, bin/%.o, $(SRC))
OBJ+=$(patsubst src/shape/%.cpp, bin/%.o, $(SHAPE_SRC))
EXE=main

CC=g++
CFLAGS=-Wall -g -O3 -std=c++0x -march=native -I./include
LDFLAGS= -lm
RM=rm

vpath %.o bin/

bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: src/shape/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: all
all: $(EXE)
	@echo Done

$(EXE): $(OBJ)
	$(CC) $(OBJ) $(LDFLAGS) -o $@
	
.PHONY: clean
clean:
	-$(RM) $(OBJ)
	@echo Clean Done!
