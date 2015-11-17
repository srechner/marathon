# Makefile for building marathon library
#
# Just run 'make cpp' to compile cpp only library (without cuda support).
#
# Run 'make cuda' to build library with cuda support.
# You may change $(COMP_CAP) to your own cuda compate capability
#

NVCC := /usr/local/cuda/bin/nvcc

COMP_CAP = \
-gencode arch=compute_52,code=sm_52 \
-gencode arch=compute_30,code=sm_30 \
-gencode arch=compute_35,code=sm_35

RM := rm -rf

# All of the sources participating in the build are defined here
CPP_SRCS = \
./src/marathon/rational.cpp \
./src/marathon/state_graph.cpp \
./src/marathon/transition.cpp \
./src/marathon/cpu/canonical_path.cpp \
./src/marathon/cpu/eigenvalues.cpp \
./src/marathon/cpu/mixing_time.cpp \
./src/marathon/cpu/transition_matrix.cpp \
./src/marathon/cpu/variation_distance.cpp \
./src/marathon/cpu/shortest_paths.cpp \
./src/chains/matching/bipartite_matching.cpp \
./src/chains/matching/broder86.cpp \
./src/chains/matching/jerrum_sinclair_vigoda04.cpp \
./src/chains/matching/sparse_bipartite_graph.cpp \
./src/chains/sequences/switch_bipartite.cpp \
./src/chains/sequences/switch_bipartite_fast.cpp \
./src/chains/sequences/curveball.cpp \
./src/chains/sequences/curveball2.cpp \
./src/chains/sequences/dense_bipartite_graph.cpp \
./src/chains/sequences/havel_hakimi.cpp 

CUDA_SRCS = \
./src/marathon/gpu/cuda_functions.cu \
./src/marathon/gpu/mixing_time.cpp \
./src/marathon/gpu/transition_matrix.cpp \
./src/marathon/gpu/variation_distance.cpp \
./src/marathon/hybrid/cuda_functions.cu \
./src/marathon/hybrid/mixing_time.cpp \
./src/marathon/hybrid/transition_matrix.cpp


CPP_OBJS = \
./src/marathon/rational.o \
./src/marathon/state_graph.o \
./src/marathon/transition.o \
./src/marathon/cpu/canonical_path.o \
./src/marathon/cpu/eigenvalues.o \
./src/marathon/cpu/mixing_time.o \
./src/marathon/cpu/transition_matrix.o \
./src/marathon/cpu/variation_distance.o \
./src/marathon/cpu/shortest_paths.o \
./src/chains/matching/bipartite_matching.o \
./src/chains/matching/broder86.o \
./src/chains/matching/jerrum_sinclair_vigoda04.o \
./src/chains/matching/sparse_bipartite_graph.o \
./src/chains/sequences/switch_bipartite.o \
./src/chains/sequences/switch_bipartite_fast.o \
./src/chains/sequences/curveball.o \
./src/chains/sequences/curveball2.o \
./src/chains/sequences/dense_bipartite_graph.o \
./src/chains/sequences/havel_hakimi.o 

CUDA_OBJS = \
./src/marathon/gpu/cuda_functions.o \
./src/marathon/gpu/mixing_time.o \
./src/marathon/gpu/transition_matrix.o \
./src/marathon/gpu/variation_distance.o \
./src/marathon/hybrid/cuda_functions.o \
./src/marathon/hybrid/mixing_time.o \
./src/marathon/hybrid/transition_matrix.o


CPP_DEPS = \
./src/marathon/rational.d \
./src/marathon/state_graph.d \
./src/marathon/transition.d \
./src/marathon/cpu/canonical_path.d \
./src/marathon/cpu/eigenvalues.d \
./src/marathon/cpu/mixing_time.d \
./src/marathon/cpu/transition_matrix.d \
./src/marathon/cpu/variation_distance.d \
./src/marathon/cpu/shortest_paths.d \
./src/chains/matching/bipartite_matching.d \
./src/chains/matching/broder86.d \
./src/chains/matching/jerrum_sinclair_vigoda04.d \
./src/chains/matching/sparse_bipartite_graph.d \
./src/chains/sequences/switch_bipartite.d \
./src/chains/sequences/switch_bipartite_fast.d \
./src/chains/sequences/curveball.d \
./src/chains/sequences/curveball2.d \
./src/chains/sequences/dense_bipartite_graph.d \
./src/chains/sequences/havel_hakimi.d 

CUDA_DEPS = \
./src/marathon/gpu/cuda_functions.d \
./src/marathon/gpu/mixing_time.d \
./src/marathon/gpu/transition_matrix.d \
./src/marathon/gpu/variation_distance.d \
./src/marathon/hybrid/cuda_functions.d \
./src/marathon/hybrid/mixing_time.d \
./src/marathon/hybrid/transition_matrix.d

# Add inputs and outputs from these tool invocations to the build variables 

all: cpp

# Compile CPU Code only (No CUDA)
cpp: $(CPP_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -shared -o "libmarathon.so" $(CPP_OBJS) 
	@echo 'Finished building target: $@'
	@echo ' '
	
cuda: $(CUDA_OBJS) $(CPP_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC Linker'
	$(NVCC) --cudart static -shared -std=c++11 --relocatable-device-code=true $(COMP_CAP) -link -o "libmarathon.so" $(CPP_OBJS) $(CUDA_OBJS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o: %.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	$(NVCC) -O2 -Xcompiler -fPIC -Xcompiler -fopenmp -std=c++11 $(COMP_CAP) -M -o "$(@:%.o=%.d)" "$<"
	$(NVCC) -O2 -Xcompiler -fPIC -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=true $(COMP_CAP) -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Each subdirectory must supply rules for building sources it contributes
%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++11 -O2 -c -fmessage-length=0 -fPIC -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(LIBRARIES) $(CC_DEPS) $(CPP_DEPS) $(C_UPPER_DEPS) $(CXX_DEPS) $(CPP_OBJS) $(CPP_DEPS) $(C_DEPS) $(CUDA_OBJS) $(CUDA_DEPS) libmarathon.so
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

