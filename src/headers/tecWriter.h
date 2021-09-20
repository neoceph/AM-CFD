#ifndef tecWriter_h
#define tecWriter_h

// standard headers
#include <vector>
#include <math.h>
#include <map>
#include <unordered_map>
#include <stdlib.h>
#include <string>

// user-defined headers
#include <PrintManager.h>

class tecWriter : public PrintManager
{
    public:
        // Constructor
        tecWriter(Mesh *meshObj,
                  DomainManager *domainMgr,
                  SolutionVariableManager *solVarMgr,
                  BoundCondManager *bcMgr,
                  CFDSolverManager *cfdSolvMgr,
                  FreeSurfaceManager *freesurMgr,
                  string &FileOut,
                  string &directoryName,
                  string &baseOutName);
    
        // Destructor
        virtual ~tecWriter() {}
    
        // Local variables
        ofstream outputFile_;
        unordered_map<string, double*> mappedScalarDataBase_, mappedVectorDataBase_;
         
        // Local functions/subroutine
        void visualOut();

        void writeTec_openFile();

        void writeTec_data();

        void writeTec_closeFile();
        
        void outputTecFile();

        void outputCustomFile();
        
        void outputOtherFile();
        
        void outputUserDefinedFile();
        
        void outputSurfaceFile();

    private: 
        // Private variables
        double outputtemp_;
          // term for determining if simulation is being initialized
        bool isinit_;
          // calculate how many grids should be output in different axis
        int gridx, gridy, gridz;  
};

#endif

