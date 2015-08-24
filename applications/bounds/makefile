
RM := rm -rf

LIB_DIR := -L../../ -L/usr/local/cuda/lib64/

LIBS := -lmarathon -lgomp -lpthread -lopenblas -larpack++ -larpack -lsuperlu -lcublas

CPP_SRCS += \
./src/bounds.cpp 

OBJS += \
./src/bounds.o 

CPP_DEPS += \
./src/bounds.d 

# All Target
all: bounds

# Tool invocations
bounds: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ $(LIB_DIR) -o "bounds" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o: %.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) bounds
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
