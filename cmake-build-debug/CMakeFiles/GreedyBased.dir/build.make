# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.17

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
CMAKE_COMMAND = /home/carlosrey/clion-2020.1.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/carlosrey/clion-2020.1.2/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/carlosrey/CLionProjects/GBQ

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/carlosrey/CLionProjects/GBQ/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/GreedyBased.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/GreedyBased.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GreedyBased.dir/flags.make

CMakeFiles/GreedyBased.dir/main.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/GreedyBased.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/main.cpp.o -c /home/carlosrey/CLionProjects/GBQ/main.cpp

CMakeFiles/GreedyBased.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/main.cpp > CMakeFiles/GreedyBased.dir/main.cpp.i

CMakeFiles/GreedyBased.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/main.cpp -o CMakeFiles/GreedyBased.dir/main.cpp.s

CMakeFiles/GreedyBased.dir/localsearch.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/localsearch.cpp.o: ../localsearch.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/GreedyBased.dir/localsearch.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/localsearch.cpp.o -c /home/carlosrey/CLionProjects/GBQ/localsearch.cpp

CMakeFiles/GreedyBased.dir/localsearch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/localsearch.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/localsearch.cpp > CMakeFiles/GreedyBased.dir/localsearch.cpp.i

CMakeFiles/GreedyBased.dir/localsearch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/localsearch.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/localsearch.cpp -o CMakeFiles/GreedyBased.dir/localsearch.cpp.s

CMakeFiles/GreedyBased.dir/utils.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/utils.cpp.o: ../utils.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/GreedyBased.dir/utils.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/utils.cpp.o -c /home/carlosrey/CLionProjects/GBQ/utils.cpp

CMakeFiles/GreedyBased.dir/utils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/utils.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/utils.cpp > CMakeFiles/GreedyBased.dir/utils.cpp.i

CMakeFiles/GreedyBased.dir/utils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/utils.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/utils.cpp -o CMakeFiles/GreedyBased.dir/utils.cpp.s

CMakeFiles/GreedyBased.dir/Solution.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/Solution.cpp.o: ../Solution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/GreedyBased.dir/Solution.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/Solution.cpp.o -c /home/carlosrey/CLionProjects/GBQ/Solution.cpp

CMakeFiles/GreedyBased.dir/Solution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/Solution.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/Solution.cpp > CMakeFiles/GreedyBased.dir/Solution.cpp.i

CMakeFiles/GreedyBased.dir/Solution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/Solution.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/Solution.cpp -o CMakeFiles/GreedyBased.dir/Solution.cpp.s

CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o: ../QMKP_01qpCPX2.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o -c /home/carlosrey/CLionProjects/GBQ/QMKP_01qpCPX2.cpp

CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/QMKP_01qpCPX2.cpp > CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.i

CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/QMKP_01qpCPX2.cpp -o CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.s

CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o: ../qkp_grdy.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o -c /home/carlosrey/CLionProjects/GBQ/qkp_grdy.cpp

CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/qkp_grdy.cpp > CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.i

CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/qkp_grdy.cpp -o CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.s

CMakeFiles/GreedyBased.dir/quadknap.c.o: CMakeFiles/GreedyBased.dir/flags.make
CMakeFiles/GreedyBased.dir/quadknap.c.o: ../quadknap.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object CMakeFiles/GreedyBased.dir/quadknap.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/GreedyBased.dir/quadknap.c.o   -c /home/carlosrey/CLionProjects/GBQ/quadknap.c

CMakeFiles/GreedyBased.dir/quadknap.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/GreedyBased.dir/quadknap.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/carlosrey/CLionProjects/GBQ/quadknap.c > CMakeFiles/GreedyBased.dir/quadknap.c.i

CMakeFiles/GreedyBased.dir/quadknap.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/GreedyBased.dir/quadknap.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/carlosrey/CLionProjects/GBQ/quadknap.c -o CMakeFiles/GreedyBased.dir/quadknap.c.s

# Object files for target GreedyBased
GreedyBased_OBJECTS = \
"CMakeFiles/GreedyBased.dir/main.cpp.o" \
"CMakeFiles/GreedyBased.dir/localsearch.cpp.o" \
"CMakeFiles/GreedyBased.dir/utils.cpp.o" \
"CMakeFiles/GreedyBased.dir/Solution.cpp.o" \
"CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o" \
"CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o" \
"CMakeFiles/GreedyBased.dir/quadknap.c.o"

# External object files for target GreedyBased
GreedyBased_EXTERNAL_OBJECTS =

GreedyBased: CMakeFiles/GreedyBased.dir/main.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/localsearch.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/utils.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/Solution.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/QMKP_01qpCPX2.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/qkp_grdy.cpp.o
GreedyBased: CMakeFiles/GreedyBased.dir/quadknap.c.o
GreedyBased: CMakeFiles/GreedyBased.dir/build.make
GreedyBased: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libilocplex.a
GreedyBased: /opt/ibm/ILOG/CPLEX_Studio1210/cplex/lib/x86-64_linux/static_pic/libcplex.a
GreedyBased: /opt/ibm/ILOG/CPLEX_Studio1210/concert/lib/x86-64_linux/static_pic/libconcert.a
GreedyBased: /home/carlosrey/gurobi903/linux64/lib/libgurobi_c++.a
GreedyBased: /home/carlosrey/gurobi903/linux64/lib/libgurobi90.so
GreedyBased: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
GreedyBased: /usr/lib/x86_64-linux-gnu/libpthread.so
GreedyBased: CMakeFiles/GreedyBased.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable GreedyBased"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GreedyBased.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GreedyBased.dir/build: GreedyBased

.PHONY : CMakeFiles/GreedyBased.dir/build

CMakeFiles/GreedyBased.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GreedyBased.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GreedyBased.dir/clean

CMakeFiles/GreedyBased.dir/depend:
	cd /home/carlosrey/CLionProjects/GBQ/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/carlosrey/CLionProjects/GBQ /home/carlosrey/CLionProjects/GBQ /home/carlosrey/CLionProjects/GBQ/cmake-build-debug /home/carlosrey/CLionProjects/GBQ/cmake-build-debug /home/carlosrey/CLionProjects/GBQ/cmake-build-debug/CMakeFiles/GreedyBased.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GreedyBased.dir/depend

