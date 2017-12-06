######
# Makefile for l2knn by David C. Anastasiu
#	Required libraries: m, gcc_s
######

# Define libraries and library directories
LIBDIRS := -L"${HOME}/programs/lib" -L/usr/lib -L/usr/local/lib
LIBS := -lgcc_s -lm 
# Source include directories
INC += -I"${HOME}/programs/include" -I/usr/local/include/
# OMP options
OMPOPTIONS = -fopenmp -D__OPENMP__
# C flags
CFLAGS += -c -O3 -fno-strict-aliasing $(OMPOPTIONS) -DLINUX -D_FILE_OFFSET_BITS=64 -std=c99 -Wall -Wno-unknown-pragmas -Wno-unused-function -Wno-unused-label -Wno-unused-variable -Wno-parentheses -Wsequence-point
# Other compile choices
DEBUG := -DNDEBUG # change -DNDEBUG to nothing to enable internal debug messages.
RM := rm -rf
EXE := pl2ap 
CC := gcc

###
# Execute make
###

# Add inputs and outputs from these tool invocations to the build variables 
HEADERS := $(shell ls ../src/*.h)
C_SRCS  := $(shell cd ../src && ls *.c)
C_OBJS  := $(C_SRCS:%.c=%.o)
C_DEPS  := $(C_SRCS:%.c=%.d)

# All Targets
all: pl2ap

# Objects depend on its source and all headers
%.o: ../src/%.c $(HEADERS)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc $(DEBUG) $(INC) $(CFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Program depends on objects
pl2ap: $(C_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc $(LIBDIRS) $(OMPOPTIONS) -o $(EXE) $(C_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Clean Target
clean:
	$(RM) *.o *.d $(EXE)
	@echo ' '

# These targets do not produce files
.PHONY: all clean

