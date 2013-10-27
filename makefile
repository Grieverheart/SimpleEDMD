SRC=$(wildcard src/*.cpp)
OBJ=$(patsubst src/%.cpp, bin/%.o, $(SRC))
EXE=main

CC=g++
CFLAGS=-Wall -g -O3 -std=c++0x -I.
LDFLAGS= -lm
RM=rm

vpath %.o bin/

bin/%.o: src/%.cpp
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
