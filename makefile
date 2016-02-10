SRC=$(wildcard src/*.cpp)
SHAPE_SRC=$(wildcard src/shape/*.cpp)
OVERLAP_SRC=$(wildcard src/overlap/*.cpp)
IO_SRC=$(wildcard src/io/*.cpp)
SERIALIZATION_SRC=$(wildcard src/serialization/*.cpp)
TINYXML_SRC=$(wildcard external/tinyxml2/*.cpp)
OBJ=$(patsubst src/%.cpp, bin/%.o, $(SRC))
OBJ+=$(patsubst src/shape/%.cpp, bin/%.o, $(SHAPE_SRC))
OBJ+=$(patsubst src/overlap/%.cpp, bin/%.o, $(OVERLAP_SRC))
OBJ+=$(patsubst src/io/%.cpp, bin/%.o, $(IO_SRC))
OBJ+=$(patsubst src/serialization/%.cpp, bin/%.o, $(SERIALIZATION_SRC))
OBJ+=$(patsubst external/tinyxml2/%.cpp, bin/%.o, $(TINYXML_SRC))
EXE=main

CC=g++
CFLAGS=-Wall -Wno-unused-function -g -O3 -std=c++11 -march=native -DNDEBUG -I./include -I./external
LDFLAGS= -lm
RM=rm

vpath %.o bin/

bin/%.o: src/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: src/shape/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: src/overlap/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: src/io/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: src/serialization/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

bin/%.o: external/tinyxml2/%.cpp
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
