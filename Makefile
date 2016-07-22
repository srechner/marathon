# Makefile for building marathon library
#
# Just run 'make' to compile cpp only library (without cuda support).
#
# Run 'make make CUDA=true' to build library with CUDA support.
# You may change $(COMP_CAP) to your own CUDA compate capability
#

CC = g++
CUDA_PATH = /usr/
NVCC = $(CUDA_PATH)/bin/nvcc

CUDA = false

# All of the sources participating in the build are defined here
CPP_SRCS = \
./src/marathon/marathon.cpp \
./src/marathon/Rational.cpp \
./src/marathon/StateGraph.cpp \
./src/marathon/Transition.cpp \
./src/marathon/MarkovChain.cpp \
./src/marathon/PathCongestion.cpp \
./src/marathon/Eigenvalues.cpp \
./src/marathon/TotalMixingTime.cpp \
./src/marathon/ShortestPaths.cpp \
./src/marathon/InitFinalize.cpp \
./src/marathon/Random.cpp \
./src/marathon/Combinatorics.cpp \
./src/marathon/TransitionMatrixCBLAS.cpp \
./src/marathon/TransitionMatrixCuBLAS.cpp \
./src/marathon/TransitionMatrixCuBLASXt.cpp \
./src/marathon/chain/matching/BipartiteMatching.cpp \
./src/marathon/chain/matching/Broder86.cpp \
./src/marathon/chain/matching/JSV04.cpp \
./src/marathon/chain/matching/SparseBipartiteGraph.cpp \
./src/marathon/chain/matching/JS89CanPath.cpp \
./src/marathon/chain/bipgraph/Curveball.cpp \
./src/marathon/chain/bipgraph/CurveballForbiddenEntries.cpp \
./src/marathon/chain/bipgraph/SwitchChain.cpp \
./src/marathon/chain/bipgraph/LindaChain.cpp \
./src/marathon/chain/bipgraph/SwitchChainBerger.cpp \
./src/marathon/chain/bipgraph/BinaryMatrix.cpp \
./src/marathon/chain/bipgraph/HavelHakimi.cpp \
./src/marathon/chain/bipgraph/KannanCanPath.cpp

CU_SRCS = \
./src/marathon/CudaVariationDistance.cu \
./src/marathon/CudaWrapper.cu \
./src/marathon/CudaInitFinalize.cu 

CPP_OBJS = $(CPP_SRCS:.cpp=.o)
CU_OBJS = $(addsuffix .o,$(basename $(CU_SRCS))) 
OBJECTS =  $(CPP_OBJS)

CFLAGS = -c -std=c++11 -O3 -fPIC -fopenmp 
LDFLAGS = -L/usr/local/lib/
LIBS = -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu 

ifeq ($(CUDA),true)
	CC = $(NVCC)
	CFLAGS = -O3 -Xcompiler -fPIC -Xcompiler -fopenmp \
	         -std=c++11 --compile --relocatable-device-code=true \
	         -DCUDA -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED 
	OBJECTS += $(CU_OBJS)
	LIBS += -lcublas
endif

RM := rm -rf

# Add inputs and outputs from these tool invocations to the build variables 

all: shared

shared: $(OBJECTS)
	$(CC) -shared -o "libmarathon.so" $(OBJECTS) 

%.o: %.cu
	$(CC) $(CFLAGS) -o "$@" "$<"

%.o: %.cpp
	$(CC) $(CFLAGS) -o "$@" "$<"

# Other Targets
clean:
	-$(RM) $(CPP_OBJS) $(CU_OBJS) libmarathon.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

