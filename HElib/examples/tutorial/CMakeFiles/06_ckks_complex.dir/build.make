# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/c/Ubu20/HElib/examples

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Ubu20/HElib/examples

# Include any dependencies generated for this target.
include tutorial/CMakeFiles/06_ckks_complex.dir/depend.make

# Include the progress variables for this target.
include tutorial/CMakeFiles/06_ckks_complex.dir/progress.make

# Include the compile flags for this target's objects.
include tutorial/CMakeFiles/06_ckks_complex.dir/flags.make

tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o: tutorial/CMakeFiles/06_ckks_complex.dir/flags.make
tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o: tutorial/06_ckks_complex.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Ubu20/HElib/examples/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o"
	cd /mnt/c/Ubu20/HElib/examples/tutorial && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o -c /mnt/c/Ubu20/HElib/examples/tutorial/06_ckks_complex.cpp

tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.i"
	cd /mnt/c/Ubu20/HElib/examples/tutorial && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Ubu20/HElib/examples/tutorial/06_ckks_complex.cpp > CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.i

tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.s"
	cd /mnt/c/Ubu20/HElib/examples/tutorial && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Ubu20/HElib/examples/tutorial/06_ckks_complex.cpp -o CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.s

# Object files for target 06_ckks_complex
06_ckks_complex_OBJECTS = \
"CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o"

# External object files for target 06_ckks_complex
06_ckks_complex_EXTERNAL_OBJECTS =

bin/06_ckks_complex: tutorial/CMakeFiles/06_ckks_complex.dir/06_ckks_complex.cpp.o
bin/06_ckks_complex: tutorial/CMakeFiles/06_ckks_complex.dir/build.make
bin/06_ckks_complex: /usr/local/helib_pack/lib/libhelib.a
bin/06_ckks_complex: /usr/local/helib_pack/lib/libntl.so
bin/06_ckks_complex: /usr/local/helib_pack/lib/libgmp.so
bin/06_ckks_complex: tutorial/CMakeFiles/06_ckks_complex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Ubu20/HElib/examples/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/06_ckks_complex"
	cd /mnt/c/Ubu20/HElib/examples/tutorial && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/06_ckks_complex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tutorial/CMakeFiles/06_ckks_complex.dir/build: bin/06_ckks_complex

.PHONY : tutorial/CMakeFiles/06_ckks_complex.dir/build

tutorial/CMakeFiles/06_ckks_complex.dir/clean:
	cd /mnt/c/Ubu20/HElib/examples/tutorial && $(CMAKE_COMMAND) -P CMakeFiles/06_ckks_complex.dir/cmake_clean.cmake
.PHONY : tutorial/CMakeFiles/06_ckks_complex.dir/clean

tutorial/CMakeFiles/06_ckks_complex.dir/depend:
	cd /mnt/c/Ubu20/HElib/examples && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Ubu20/HElib/examples /mnt/c/Ubu20/HElib/examples/tutorial /mnt/c/Ubu20/HElib/examples /mnt/c/Ubu20/HElib/examples/tutorial /mnt/c/Ubu20/HElib/examples/tutorial/CMakeFiles/06_ckks_complex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tutorial/CMakeFiles/06_ckks_complex.dir/depend

