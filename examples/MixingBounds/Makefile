
RM := rm -rf

LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu

SRCS = mixing_bounds.cpp 

OBJS = mixing_bounds.o 

DEPS = mixing_bounds.d 

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: cpp

# Tool invocations
cuda: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L"../../" -L"/usr/local/cuda/lib64/" -o "mixing_bounds" $(OBJS) $(LIBS) -lcublas
	@echo 'Finished building target: $@'
	@echo ' '

cpp: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L"../../" -o "mixing_bounds" $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"../../include" -O2 -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS)$(OBJS) mixing_bounds
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

