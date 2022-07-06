find . -name "CMakeFiles" -exec rm -rf {} \;
rm -f CMakeCache.txt

cmake \
    -DCMAKE_BUILD_TYPE=Release \
    ..
