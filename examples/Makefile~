
RM := rm -rf

LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu -lcublas

CPP_SRCS = totalMixingTime.cpp 

OBJS = totalMixingTime.o 

CPP_DEPS = totalMixingTime.d 

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: totalMixingTime

# Tool invocations
totalMixingTime: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L"../" -L"/usr/local/cuda/lib64/" -o "totalMixingTime" $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I"../include" -O2 -c -fmessage-length=0 -std=c++11 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) totalMixingTime
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

