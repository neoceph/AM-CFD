# Northwestern AM-CFD Code (C++) 
BetaVersion
Developed by Abdullah Al Amin. email: abdullah.amin@northwestern.edu 

Initially started by Kevontrez Jones and Zhengtao Gan

## Installation 
---
1. Create a build directory {"mkdir build"}
2. Copy "configure.sh" into the build directory {"cp configure.sh build/"}
3. Go into the build directory {"cd build/"}
4. Execute the configure.sh file in the build directory {"./configure"}
5. Compile the code {"make -j4"}

## OpenMp Version 4.0 or greater is REQUIRED!

## Building the distribution package
1. Build the project by either clicking the 'Build' button on the VSCode bottom ribbon or the usual CMake build process.
2. go into the ./build directory and execute 'cpack' comman. It will generate the 'ProgramName-version.number-Linux.sh'. 