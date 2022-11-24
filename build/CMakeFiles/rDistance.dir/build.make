# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/ruohawang2/miniconda3/envs/pytorch_gpu/bin/cmake

# The command to remove a file.
RM = /home/ruohawang2/miniconda3/envs/pytorch_gpu/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ruohawang2/12.Phage_assem/pipeline/seqGraph

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build

# Include any dependencies generated for this target.
include CMakeFiles/rDistance.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/rDistance.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/rDistance.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/rDistance.dir/flags.make

CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o: CMakeFiles/rDistance.dir/flags.make
CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o: ../utils/reads_overlap.cpp
CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o: CMakeFiles/rDistance.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o -MF CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o.d -o CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o -c /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/utils/reads_overlap.cpp

CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/utils/reads_overlap.cpp > CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.i

CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/utils/reads_overlap.cpp -o CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.s

# Object files for target rDistance
rDistance_OBJECTS = \
"CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o"

# External object files for target rDistance
rDistance_EXTERNAL_OBJECTS =

rDistance: CMakeFiles/rDistance.dir/utils/reads_overlap.cpp.o
rDistance: CMakeFiles/rDistance.dir/build.make
rDistance: CMakeFiles/rDistance.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable rDistance"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/rDistance.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/rDistance.dir/build: rDistance
.PHONY : CMakeFiles/rDistance.dir/build

CMakeFiles/rDistance.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/rDistance.dir/cmake_clean.cmake
.PHONY : CMakeFiles/rDistance.dir/clean

CMakeFiles/rDistance.dir/depend:
	cd /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ruohawang2/12.Phage_assem/pipeline/seqGraph /home/ruohawang2/12.Phage_assem/pipeline/seqGraph /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build /home/ruohawang2/12.Phage_assem/pipeline/seqGraph/build/CMakeFiles/rDistance.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/rDistance.dir/depend
