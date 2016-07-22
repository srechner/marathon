# Makefile for building marathon library
#
# Just run 'make' to compile cpp only library (without cuda support).
#
# Run 'make make CUDA=true' to build library with CUDA support.
# You may change $(COMP_CAP) to your own CUDA compate capability
#

CC = g++
CUDA_PATH = /usr/local/cuda/

RM := rm -rf

CUDA = false

EXECUTABLE = "RandomSamples"
CFLAGS = -c -std=c++11 -O3 -fopenmp -I"../../include"
LDFLAGS = -L"../.." 
LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu

SRCS = RandomSamples.cpp 
OBJS = $(SRCS:.cpp=.o) 
DEPS = $(SRCS:.cpp=.d)

ifeq ($(CUDA),true)
	LDFLAGS += -L"$(CUDA_PATH)/lib64/"
	LIBS += -lcublas
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: $(OBJS)
	$(CC) $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -o "$@" "$<"

# Other Targets
clean:
	-$(RM) $(EXECUTABLE) $(OBJS) $(DEPS)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

