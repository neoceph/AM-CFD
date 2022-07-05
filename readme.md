# A comprehensive framework for C++ development
- This CMake structure includes googletest as an external library using git submodule and organizes unittests in a separate executable.
- Organizes codes on a separate src directory where headers and sources have corresponding headers and sources directory.
- Developed on Ubuntu 18.04

## Requirements
- CMake 3.21 or above
- Doxygen 1.9.4 and above
    -`sudo apt update` 
    - `sudo apt install doxygen`
- Sphinx 4.4.0 and above
    `sudo apt-get install python3-sphinx`
- Sphinx readthedocs theme
    `pip3 install sphinx_rtd_theme`
- Breathe
    - `pip3 install breathe`
    - `pip3 install sphinx-sitemap`
- Latexmk (required for pdf generation)
    - If not available can be installed with `sudo apt-get update -y` and next `sudo apt-get install -y latexmk`
    - For .sty file not found use `apt install texlive-latex-extra`
- Rpmbuild (required for .rpm package generation) 
    - Creating rpm will require an rpmbuild executable. To get rpm installed issue `sudo apt-get install rpm`
- NSIS 3.0 and up (required for windows installer package generation)
    - Creating NSIS will require an nsis executable. In Ubuntu-18.04 nsis installer in the form of .deb might be required as it doesn't come with the distribution. The installation file could be found [here](https://ubuntu.pkgs.org/18.04/ubuntu-universe-amd64/nsis_2.51-1_amd64.deb.html). It might also require a dependency on nsis-common which can be downloaded from [here](http://ftp.de.debian.org/debian/pool/main/n/nsis/nsis-common_3.08-2_all.deb). The debian packages are searched and found on the website [here](https://debian.pkgs.org/). Installation of .deb package is done using the command `sudo dpkg -i DEB_PACKAGE_NAME`
## Packing installer
1. Assigning a git tag is easy 
    - `git tag -a v1.0.3 -m 'Version 1.0.3 - Feature 1'`
    - `git tag` lists all the available tags.
    - Check other `git tag` details [here].(https://git-scm.com/book/en/v2/Git-Basics-Tagging)
2. Build the code 
    - create a directory named 'build' with `mkdir build`
    - change to build directory with `cd build`
    - copy the configure.sh to build directory with `cp ../configure.sh .`
    - execute configure .sh with `./configure.sh`. You may need to change executable permission with `chmod +x configure.sh`
    - compile the code `make -j $(nproc)`. The  `$(nproc)` part of the code identifies the number of cpus and use all of them to compile the code.
3. From within the build directory issue `cpack` command to create installer. All the .sh, .tar.gz, and .tar.Z files are the installer in different format.

## Changing the name of the project
- Project name needs to be changed in two places
    1. Line #3 of the root directory CMakeLists.txt file `set(Project Project-Name)`. The 'Project-Name' can be changed to any appropriate name.
    2. Line #35 `PROJECT_NAME = "AM-CFD"` of the [Doxyfile.in](./documentation/Doxyfile.in) located in the *./documentation* directory. 

## Generating documentation
- Under the build directory after a successful build `make docs` will generate api documentation under doc/api directory
    - in the doc/api/latex/ directory issuing `make pdf` will generate pdf from the generated latex files
- Under the doc/usermanual directory issuing `make html` will generate html documentation under doc/readthedocs/_build/html directory
- `make latexpdf` will generate pdf documentation under doc/readthedocs/_build/latex directory

    