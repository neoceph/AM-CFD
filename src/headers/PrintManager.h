#ifndef PrintManager_h
#define PrintManager_h

// standard headers
#include <vector>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h> 
#include <string>

// definitions for memory monitoring
#if defined(__APPLE__)
#include<mach/task.h>
#include<mach/mach_init.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#endif
#ifdef __linux__
#  ifdef BGQ_LWK
#    include <spi/include/kernel/memory.h>
#    include <spi/include/kernel/location.h>
#  else
#    define PROCFS
#  endif
#endif

// user-defined headers
#include <Mesh.h>
#include <DomainManager.h>
#include <BoundCondManager.h>
#include <SolutionVariableManager.h>
#include <CFDSolverManager.h>
#include <FreeSurfaceManager.h>

// define some classes to be used
class Mesh;
class DomainManager;
class BoundCondManager;
class SolutionVariableManager;
class CFDSolverManager;
class FreeSurfaceManager;

class PrintManager
{
    public:
        // Constructor
        PrintManager(Mesh *meshObj,
                     DomainManager *domainMgr,
                     SolutionVariableManager *solVarMgr,
                     BoundCondManager *bcMgr,
                     CFDSolverManager *cfdSolvMgr,
                     FreeSurfaceManager *freesurMgr,
                     string &FileOut,
                     string &directoryName,
                     string &baseOutName);

        // Destructor
        virtual ~PrintManager() {};
         
        // Input
        Mesh *meshObj_;
        DomainManager *domainMgr_;
        SolutionVariableManager *solVarMgr_;
        BoundCondManager *bcMgr_;
        CFDSolverManager *cfdSolvMgr_;
        FreeSurfaceManager *freesurMgr_;
        string &FileOut_, &directoryName_, &baseOutName_;
         
        // Local variables
          // solution variables
        double *temp_;
        double *uvel_;
        double *vvel_;
        double *wvel_;
        double *rate_;
        double *grad_;
        double *tempr_;
          // geometric variables
        double *x_; 
        double *y_; 
        double *z_; 
        double *zr_;
          // velocity field at central nodes
        double auvel_[_nz_*_ny_*_nx_];
        double avvel_[_nz_*_ny_*_nx_];
        double awvel_[_nz_*_ny_*_nx_];
        int gridz_;
          // outputs
        vector<string> outScalarNames_, outVectorNames_, outCellNames_;
          // timer used outputting section file between main outputs
        double secttime_ = 0;
         
        // Local functions/subroutines
        void initialization();

        void scanOutputData();
        
        void outputResiduals();

        void outputSection();
        
        void outputFinalSteady();

        void outputFinalGeometry();
        
        void endTime(double &comptime, double &dumptime);

        virtual void visualOut() = 0;

    private: 
        // Private variables
          // definition for consistent easy printing
        const std::string SecondEleWidth_  = "  ";
          // max values for velocity in each direction
        double umax_, vmax_, wmax_;
          // timer for outputting an extracted section
        double tracksect_ = 0.0;
        double dumpsect_ = 0.0;
          // residual, final geometry, and extraction section files
        ofstream residFile_, finalgeomFile_, sectionFile_;
          // location for the height of the substate for powder-bed
        int baseheightk_;

        // Private functions/subroutines
        void outputInitHeader();

        void outputAssignInputTimer();

        void outputGridTime();

        void outputTimestepStats();

        void outputSolverMode();

        void outputSectionInfo();

        void outputInitFooter();

        void initializePrinter();

        void outputFinalHeader();

        void outputInfoCPU(double &comptime, double &dumptime);

        void outputInfoMemory();

        double reportMemoryUsed();

        double reportMaxMemoryUsed();

        void getMemoryAvailable(double &avail);

        void outputFinalFooter();
};

#endif
