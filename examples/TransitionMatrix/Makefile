
RM := rm -rf

LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu

SRCS = transition_matrix.cpp 

OBJS = transition_matrix.o 

DEPS = transition_matrix.d 

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: cpp

# Tool invocations
cpp: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L"../../" -o "transition_matrix" $(OBJS) $(LIBS)
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
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS)$(OBJS) transition_matrix
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

