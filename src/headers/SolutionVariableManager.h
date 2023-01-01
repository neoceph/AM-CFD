#ifndef SolutionVariableManager_h
#define SolutionVariableManager_h
 
// standard headers
#include <vector>
#include <Mesh.h>

class SolutionVariableManager
{
    public:
        // Constructor
        SolutionVariableManager();

        // Destructor
        virtual ~SolutionVariableManager() {}

        // Local variables
        map < string , double*> mapSolutionVars_;
        map < string , double*> mapGeometryVars_;
};

#endif
