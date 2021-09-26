CC=g++
OPENMPI=mpic++
CFLAGS=-O3
BUILD_DIR=build
SRC_DIR=src
INCLUDE_DIR=./include
SOURCES := $(shell find $(SRC_DIR) -name '*.cpp')

$(info $(shell mkdir -p $(BUILD_DIR)))

default:
	$(CC) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS)
	./build/main
	@printf "\n"

openmp:
	$(CC) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS) -fopenmp
	./build/main
	@printf "\n"

mpi:
	$(OPENMPI) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS) 
	mpirun -np 4 ./build/main 

hybrid:
	$(OPENMPI) -o $(BUILD_DIR)/main -I$(INCLUDE_DIR) $(SOURCES) $(CFLAGS) -fopenmp
	mpirun -np 4 ./build/main

clean:
	rm test