#ifndef VERSION_H
#define VERSION_H

#include <string>

extern const std::string gGIT_VERSION;
extern const std::string gGIT_VERSION_SHORT;


class VersionInfo
{
    public:
        void getVersionInfo(char* arg);
};

#endif //VERSION_H
