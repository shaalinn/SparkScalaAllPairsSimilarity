# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.9

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.9.4_1/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.9.4_1/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/samin/Downloads/l2ap/build/GKlib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64

# Include any dependencies generated for this target.
include test/CMakeFiles/grKx.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/grKx.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/grKx.dir/flags.make

test/CMakeFiles/grKx.dir/grKx.c.o: test/CMakeFiles/grKx.dir/flags.make
test/CMakeFiles/grKx.dir/grKx.c.o: ../../test/grKx.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object test/CMakeFiles/grKx.dir/grKx.c.o"
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/grKx.dir/grKx.c.o   -c /Users/samin/Downloads/l2ap/build/GKlib/test/grKx.c

test/CMakeFiles/grKx.dir/grKx.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/grKx.dir/grKx.c.i"
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/samin/Downloads/l2ap/build/GKlib/test/grKx.c > CMakeFiles/grKx.dir/grKx.c.i

test/CMakeFiles/grKx.dir/grKx.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/grKx.dir/grKx.c.s"
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test && /Library/Developer/CommandLineTools/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/samin/Downloads/l2ap/build/GKlib/test/grKx.c -o CMakeFiles/grKx.dir/grKx.c.s

test/CMakeFiles/grKx.dir/grKx.c.o.requires:

.PHONY : test/CMakeFiles/grKx.dir/grKx.c.o.requires

test/CMakeFiles/grKx.dir/grKx.c.o.provides: test/CMakeFiles/grKx.dir/grKx.c.o.requires
	$(MAKE) -f test/CMakeFiles/grKx.dir/build.make test/CMakeFiles/grKx.dir/grKx.c.o.provides.build
.PHONY : test/CMakeFiles/grKx.dir/grKx.c.o.provides

test/CMakeFiles/grKx.dir/grKx.c.o.provides.build: test/CMakeFiles/grKx.dir/grKx.c.o


# Object files for target grKx
grKx_OBJECTS = \
"CMakeFiles/grKx.dir/grKx.c.o"

# External object files for target grKx
grKx_EXTERNAL_OBJECTS =

test/grKx: test/CMakeFiles/grKx.dir/grKx.c.o
test/grKx: test/CMakeFiles/grKx.dir/build.make
test/grKx: libGKlib.a
test/grKx: test/CMakeFiles/grKx.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable grKx"
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/grKx.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/grKx.dir/build: test/grKx

.PHONY : test/CMakeFiles/grKx.dir/build

test/CMakeFiles/grKx.dir/requires: test/CMakeFiles/grKx.dir/grKx.c.o.requires

.PHONY : test/CMakeFiles/grKx.dir/requires

test/CMakeFiles/grKx.dir/clean:
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test && $(CMAKE_COMMAND) -P CMakeFiles/grKx.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/grKx.dir/clean

test/CMakeFiles/grKx.dir/depend:
	cd /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/samin/Downloads/l2ap/build/GKlib /Users/samin/Downloads/l2ap/build/GKlib/test /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64 /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test /Users/samin/Downloads/l2ap/build/GKlib/build/Darwin-x86_64/test/CMakeFiles/grKx.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/grKx.dir/depend

