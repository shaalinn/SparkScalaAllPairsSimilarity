######
# Makefile for apss by David C. Anastasiu
# Dependencies: perl (for CDEFINES parsing below - feel free to remove dependency)
#	and the following libraries: GKlib, m
######

###
# Compile choices
###

# Define libraries and library directories
LIBDIRS := -L"${HOME}/programs/lib" -L/usr/lib -L/usr/local/lib
LIBS := -lGKlib -lm
# Source include directories
INC += -I"${HOME}/programs/include" -I/usr/local/include/
# C flags
CFLAGS += -O3 -c -fno-strict-aliasing -DLINUX -D_FILE_OFFSET_BITS=64 -std=c99 -Wall -Wno-unknown-pragmas -Wno-unused-function -Wno-unused-label -Wno-unused-variable -Wno-parentheses -Wsequence-point
# Other compile choices
DEBUG := -DNDEBUG # change to nothing to enable internal debug messages.
RM := rm -rf
EXE := apss 
CC := gcc

###
# Execute make
###

# Parse pruning parameters
CDEFS = L2PS RS4 L2CG PSCV DP5 L2CV TL1 TL2 MKPL
CDEFINES := $(shell perl -e '@a=split(/\s+/,"${CDEFS}"); for($$i=0; $$i<scalar(@a); $$i++){$$a[$$i]="-D".uc($$a[$$i]);} print join(" ",@a);')

# Add inputs and outputs from these tool invocations to the build variables 
HEADERS := $(shell ls ../*.h)
C_SRCS := $(shell cd ../ && ls *.c)
C_OBJS := $(C_SRCS:%.c=%.o)
C_DEPS := $(C_SRCS:%.c=%.d)

# All Targets
all: apss

# Objects depend on its source and all headers
%.o: ../%.c $(HEADERS)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc $(CDEFINES) $(DEBUG) $(INC) $(CFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# Program depends on objects
apss: $(C_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc $(LIBDIRS) -o $(EXE) $(C_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Clean Target
clean:
	$(RM) *.o *.d $(EXE)
	@echo ' '

# These targets do not produce files
.PHONY: all clean

