# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_COMMAND = /opt/clion-2018.3.4/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion-2018.3.4/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/SSE.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SSE.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SSE.dir/flags.make

CMakeFiles/SSE.dir/main.c.o: CMakeFiles/SSE.dir/flags.make
CMakeFiles/SSE.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/SSE.dir/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/SSE.dir/main.c.o   -c /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/main.c

CMakeFiles/SSE.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/SSE.dir/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/main.c > CMakeFiles/SSE.dir/main.c.i

CMakeFiles/SSE.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/SSE.dir/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/main.c -o CMakeFiles/SSE.dir/main.c.s

# Object files for target SSE
SSE_OBJECTS = \
"CMakeFiles/SSE.dir/main.c.o"

# External object files for target SSE
SSE_EXTERNAL_OBJECTS =

SSE: CMakeFiles/SSE.dir/main.c.o
SSE: CMakeFiles/SSE.dir/build.make
SSE: CMakeFiles/SSE.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable SSE"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/SSE.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SSE.dir/build: SSE

.PHONY : CMakeFiles/SSE.dir/build

CMakeFiles/SSE.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/SSE.dir/cmake_clean.cmake
.PHONY : CMakeFiles/SSE.dir/clean

CMakeFiles/SSE.dir/depend:
	cd /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/SSE/cmake-build-debug/CMakeFiles/SSE.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SSE.dir/depend

