# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/kevinsun/Projects/Math_calculation/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/kevinsun/Projects/Math_calculation/build

# Include any dependencies generated for this target.
include cal_kin/CMakeFiles/cal_kin.dir/depend.make

# Include the progress variables for this target.
include cal_kin/CMakeFiles/cal_kin.dir/progress.make

# Include the compile flags for this target's objects.
include cal_kin/CMakeFiles/cal_kin.dir/flags.make

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o: cal_kin/CMakeFiles/cal_kin.dir/flags.make
cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o: /home/kevinsun/Projects/Math_calculation/src/cal_kin/src/cal_kin.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/kevinsun/Projects/Math_calculation/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o"
	cd /home/kevinsun/Projects/Math_calculation/build/cal_kin && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o -c /home/kevinsun/Projects/Math_calculation/src/cal_kin/src/cal_kin.cpp

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cal_kin.dir/src/cal_kin.cpp.i"
	cd /home/kevinsun/Projects/Math_calculation/build/cal_kin && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/kevinsun/Projects/Math_calculation/src/cal_kin/src/cal_kin.cpp > CMakeFiles/cal_kin.dir/src/cal_kin.cpp.i

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cal_kin.dir/src/cal_kin.cpp.s"
	cd /home/kevinsun/Projects/Math_calculation/build/cal_kin && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/kevinsun/Projects/Math_calculation/src/cal_kin/src/cal_kin.cpp -o CMakeFiles/cal_kin.dir/src/cal_kin.cpp.s

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.requires:
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.requires

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.provides: cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.requires
	$(MAKE) -f cal_kin/CMakeFiles/cal_kin.dir/build.make cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.provides.build
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.provides

cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.provides.build: cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o

# Object files for target cal_kin
cal_kin_OBJECTS = \
"CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o"

# External object files for target cal_kin
cal_kin_EXTERNAL_OBJECTS =

/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: cal_kin/CMakeFiles/cal_kin.dir/build.make
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/libroscpp.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_signals.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/librosconsole.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/librosconsole_log4cxx.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/librosconsole_backend_interface.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/liblog4cxx.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_regex.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/libxmlrpcpp.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/libroscpp_serialization.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/librostime.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_date_time.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /opt/ros/indigo/lib/libcpp_common.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_system.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libboost_thread.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libpthread.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: /usr/lib/x86_64-linux-gnu/libconsole_bridge.so
/home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin: cal_kin/CMakeFiles/cal_kin.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable /home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin"
	cd /home/kevinsun/Projects/Math_calculation/build/cal_kin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cal_kin.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
cal_kin/CMakeFiles/cal_kin.dir/build: /home/kevinsun/Projects/Math_calculation/devel/lib/cal_kin/cal_kin
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/build

cal_kin/CMakeFiles/cal_kin.dir/requires: cal_kin/CMakeFiles/cal_kin.dir/src/cal_kin.cpp.o.requires
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/requires

cal_kin/CMakeFiles/cal_kin.dir/clean:
	cd /home/kevinsun/Projects/Math_calculation/build/cal_kin && $(CMAKE_COMMAND) -P CMakeFiles/cal_kin.dir/cmake_clean.cmake
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/clean

cal_kin/CMakeFiles/cal_kin.dir/depend:
	cd /home/kevinsun/Projects/Math_calculation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/kevinsun/Projects/Math_calculation/src /home/kevinsun/Projects/Math_calculation/src/cal_kin /home/kevinsun/Projects/Math_calculation/build /home/kevinsun/Projects/Math_calculation/build/cal_kin /home/kevinsun/Projects/Math_calculation/build/cal_kin/CMakeFiles/cal_kin.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : cal_kin/CMakeFiles/cal_kin.dir/depend

