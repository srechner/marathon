# Makefile for building marathon library
#
# Just run 'make' to compile cpp only library (without cuda support).
#
# Run 'make make GPU=true' to build library with cuda support.
# You may change $(COMP_CAP) to your own cuda compate capability
#

COMP_CAP = \
-gencode arch=compute_20,code=sm_20 \
-gencode arch=compute_30,code=sm_30 \
-gencode arch=compute_35,code=sm_35 \
-gencode arch=compute_50,code=sm_50 \
-gencode arch=compute_52,code=sm_52

# All of the sources participating in the build are defined here
CPP_SRCS = \
./src/marathon/marathon_nocuda.cpp \
./src/marathon/common/rational.cpp \
./src/marathon/common/state_graph.cpp \
./src/marathon/common/transition.cpp \
./src/marathon/cpu/canonical_path.cpp \
./src/marathon/cpu/eigenvalues.cpp \
./src/marathon/cpu/mixing_time.cpp \
./src/marathon/cpu/transition_matrix.cpp \
./src/marathon/cpu/variation_distance.cpp \
./src/marathon/cpu/shortest_paths.cpp \
./src/marathon/chains/matching/bipartite_matching.cpp \
./src/marathon/chains/matching/matching_chain_JS89.cpp \
./src/marathon/chains/matching/matching_chain_JSV04.cpp \
./src/marathon/chains/matching/sparse_bipartite_graph.cpp \
./src/marathon/chains/sequences/switch_chain_bipartite.cpp \
./src/marathon/chains/sequences/switch_chain_bipartite_berger.cpp \
./src/marathon/chains/sequences/dense_bipartite_graph.cpp \
./src/marathon/chains/sequences/havel_hakimi.cpp 

CU_SRCS = \
./src/marathon/marathon_cuda.cpp \
./src/marathon/common/rational.cpp \
./src/marathon/common/state_graph.cpp \
./src/marathon/common/transition.cpp \
./src/marathon/cpu/canonical_path.cpp \
./src/marathon/cpu/eigenvalues.cpp \
./src/marathon/cpu/mixing_time.cpp \
./src/marathon/cpu/transition_matrix.cpp \
./src/marathon/cpu/variation_distance.cpp \
./src/marathon/cpu/shortest_paths.cpp \
./src/marathon/chains/matching/bipartite_matching.cpp \
./src/marathon/chains/matching/matching_chain_JS89.cpp \
./src/marathon/chains/matching/matching_chain_JSV04.cpp \
./src/marathon/chains/matching/sparse_bipartite_graph.cpp \
./src/marathon/chains/sequences/switch_chain_bipartite.cpp \
./src/marathon/chains/sequences/switch_chain_bipartite_berger.cpp \
./src/marathon/chains/sequences/dense_bipartite_graph.cpp \
./src/marathon/chains/sequences/havel_hakimi.cpp \
./src/marathon/gpu/cuda_functions.cu \
./src/marathon/gpu/init_finalize.cpp \
./src/marathon/gpu/mixing_time.cpp \
./src/marathon/gpu/transition_matrix.cpp \
./src/marathon/gpu/variation_distance.cpp \
./src/marathon/hybrid/cuda_functions.cu \
./src/marathon/hybrid/init_finalize.cpp \
./src/marathon/hybrid/mixing_time.cpp \
./src/marathon/hybrid/transition_matrix.cpp

CPP_OBJECTS := $(CPP_SRCS:.cpp=.o) 
CU_OBJECTS := $(addsuffix .o,$(basename $(CU_SRCS)))

GPU = false

CC = g++
NVCC = /usr/local/cuda/bin/nvcc
SRCS = $(CPP_SRCS)
OBJECTS = $(CPP_OBJECTS)

CFLAGS = -c -std=c++11 -O3 -fPIC -fopenmp 
LDFLAGS = -L/usr/local/lib/
LIBS = -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu 

ifeq ($(GPU),true)
	CC = $(NVCC)
	CFLAGS = -O3 -Xcompiler -fPIC -Xcompiler -fopenmp \
	         -std=c++11 --compile --relocatable-device-code=true $(COMP_CAP)
	SRCS = $(CU_SRCS)
	OBJECTS = $(CU_OBJECTS) 
	LIBS += -lcublas
endif

RM := rm -rf

# Add inputs and outputs from these tool invocations to the build variables 

all: shared

static: $(OBJECTS)
	ar rcs libmarathon.a $(OBJECTS)

shared: $(OBJECTS)
	$(CC) -shared -o "libmarathon.so" $(OBJECTS) 

%.o: %.cu
	$(CC) $(CFLAGS) -o "$@" "$<"

%.o: %.cpp
	$(CC) $(CFLAGS) -o "$@" "$<"

# Other Targets
clean:
	-$(RM) $(CPP_OBJECTS) $(CU_OBJECTS) libmarathon.a libmarathon.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

