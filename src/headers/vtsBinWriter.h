#ifndef vtsBinWriter_h
#define vtsBinWriter_h

// standard headers
#include <vector>
#include <math.h>
#include <map>
#include <unordered_map>
#include <stdlib.h>
#include <string>

// user-defined headers
#include <PrintManager.h>

class vtsBinWriter : public PrintManager
{
    public:
        // Constructor
        vtsBinWriter(Mesh *meshObj,
                     DomainManager *domainMgr,
                     SolutionVariableManager *solVarMgr,
                     BoundCondManager *bcMgr,
                     CFDSolverManager *cfdSolvMgr,
                     FreeSurfaceManager *freesurMgr,
                     string &FileOut,
                     string &directoryName,
                     string &baseOutName);

        // Destructor
        virtual ~vtsBinWriter() {}
    
        // Local variables
        string pvdName_;
        ofstream outputFile_, outputPVTS_, outputPVD_;
        vector<string> pvtsOutList_;
        vector<double> timeStampList_;
        int offSetCtr_, nPoints_, nCells_, byteCtr_, nDim_;
        unordered_map<string, double*> mappedScalarDataBase_, mappedVectorDataBase_;
    
        // Local functions/subroutines
        virtual void visualOut();
    
    private: 
        // Private functions/subroutines
        void writeVTS_coordinates();
    
        void writeVTS_coordsHeader();
    
        void writeVTS_coordsEnd();
    
        void writeVTS_pointData(double *pointVec);
    
        void writeVTS_header();
    
        void writeVTS_end();
    
        void writeVTS_pointDataHeader();
    
        void writeVTS_pointDataEnd();
    
        void writeVTS_pointOutHeader(string &outScalarName);
    
        void writeVTS_pointOutEnd();
    
        void writeVTS_appendHeader();
    
        void writeVTS_appendEnd();
    
        void initializePVD();
    
        void appendPVD();
    
        void closePVD();
    
        void writeVTS_cellDataHeader();
    
        void writeVTS_cellOutHeader();
    
        void writeVTS_cellDataEnd();
    
        void writeVTS_cellData();
    
        void writeVTS_pointOutVecHeader(string &outVectorName);
    
        void writeVTS_pointVecData(double *pointVec);

        void transferVelNodes();

        void writeVTS_GhostPoints();
};

#endif
