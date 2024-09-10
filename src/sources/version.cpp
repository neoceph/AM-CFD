#include <../version.h>
#include <iostream>
#include <string.h>


void VersionInfo::getVersionInfo(char* arg)
{
    if (strcmp(arg, "--version") == 0)
    {
        std::cout << gGIT_VERSION << std::endl;
    }

    else if (strcmp(arg, "-v") == 0)
    {
        std::cout << gGIT_VERSION_SHORT << std::endl;
    }   
}