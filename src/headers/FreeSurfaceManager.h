#ifndef FreeSurfaceManager_h
#define FreeSurfaceManager_h

// standard headers
#include <vector>

// user-defined headers
#include <Mesh.h>
#include <DomainManager.h>
#include <SolutionVariableManager.h>
#include <BoundCondManager.h>

// define some classes to be used
class Mesh;
class DomainManager;
class SolutionVariableManager;
class BoundCondManager;

class FreeSurfaceManager
{
    public:
        // Constructor
        FreeSurfaceManager(Mesh *meshObj,
                           DomainManager *domainMgr,
                           SolutionVariableManager *solVarMgr,
                           BoundCondManager *bcMgr);
    
        // Destructor
        virtual ~FreeSurfaceManager() {};

        // Input
        Mesh *meshObj_;
        DomainManager *domainMgr_;
        SolutionVariableManager *solVarMgr_;
        BoundCondManager *bcMgr_;
         
        // Local variables
          // solution variables
        double *temp_;
          // geometric variables
        double *areaij_;
        double *xu_;
        double *yv_;
        double *x_; 
        double *y_; 
        double *z_; 

        // Local functions/subroutines
        void execute();

        void initializeFreeSurface();
    private:
        // Private variables
          // surface defomation (Z minus is positive direction), previous surf defomation
        double fi_[_ny_*_nx_]; 
        double zzm_[_ny_*_nx_]; 
        double apdis_[_ny_*_nx_]; 
        double tempr_[_nz_*_ny_*_nx_];
        double zr_[_nz_*_ny_*_nx_];
          // file for ouputting data
        ofstream freesurfFile_;

        // Private functions/subroutines
        double getDeltaVol(double &alambda);

        void ebal(double &alambda);
};

#endif
