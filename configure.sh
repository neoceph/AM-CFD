# Cleanup old cache before we configure
# Note:  This does not remove files produced by make.  Use "make clean" for this.
# "-DCMAKE_BUILD_TYPE=Debug " for debug and "-DCMAKE_BUILD_TYPE=Release " for release
find . -name "CMakeFiles" -exec rm -rf {} \;
rm -f CMakeCache.txt

cmake \
    -DCMAKE_BUILD_TYPE=Release \
    ..

