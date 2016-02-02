
RM := rm -rf

GPU = false

EXECUTABLE = "States"
CFLAGS = -c -std=c++11 -O3 -fopenmp -I"../../include"
LDFLAGS = -L"../.."
LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu

SRCS = states.cpp 
OBJS = $(SRCS:.cpp=.o) 
DEPS = $(SRCS:.cpp=.d)

ifeq ($(GPU),true)
	LDFLAGS += -L"/usr/local/cuda/lib64/"
	LIBS += -lcublas
endif

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: $(OBJS)
	g++ $(LDFLAGS) -o $(EXECUTABLE) $(OBJS) $(LIBS)

%.o: %.cpp
	g++ $(CFLAGS) -o "$@" "$<"

# Other Targets
clean:
	-$(RM) $(EXECUTABLE) $(OBJS) $(DEPS)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

