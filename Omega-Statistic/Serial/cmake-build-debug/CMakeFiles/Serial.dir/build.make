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
CMAKE_SOURCE_DIR = /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Serial.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Serial.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Serial.dir/flags.make

CMakeFiles/Serial.dir/main.c.o: CMakeFiles/Serial.dir/flags.make
CMakeFiles/Serial.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/Serial.dir/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/Serial.dir/main.c.o   -c /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/main.c

CMakeFiles/Serial.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/Serial.dir/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/main.c > CMakeFiles/Serial.dir/main.c.i

CMakeFiles/Serial.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/Serial.dir/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/main.c -o CMakeFiles/Serial.dir/main.c.s

# Object files for target Serial
Serial_OBJECTS = \
"CMakeFiles/Serial.dir/main.c.o"

# External object files for target Serial
Serial_EXTERNAL_OBJECTS =

Serial: CMakeFiles/Serial.dir/main.c.o
Serial: CMakeFiles/Serial.dir/build.make
Serial: CMakeFiles/Serial.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable Serial"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Serial.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Serial.dir/build: Serial

.PHONY : CMakeFiles/Serial.dir/build

CMakeFiles/Serial.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Serial.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Serial.dir/clean

CMakeFiles/Serial.dir/depend:
	cd /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug /home/stefanos/TUC_Projects/TUC_Parallel_Computer_Architecture/Omega-Statistic/Serial/cmake-build-debug/CMakeFiles/Serial.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Serial.dir/depend

