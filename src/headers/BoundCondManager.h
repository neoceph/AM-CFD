#ifndef BoundCondManager_h
#define BoundCondManager_h

// standard headers

// user-defined headers
#include <Mesh.h>
#include <DomainManager.h>
#include <SolutionVariableManager.h>

// define some classes to be used
class Mesh;
class DomainManager;
class SolutionVariableManager;

class BoundCondManager
{
    public:
        // Constructor
        BoundCondManager(Mesh *meshObj,
                         DomainManager *domainMgr,
                         SolutionVariableManager *solVarManager);
    
        // Destructor
        virtual ~BoundCondManager() {};

        // Input
        Mesh *meshObj_;
        DomainManager *domainMgr_;
        SolutionVariableManager *solVarMgr_;

        // Local variables
          // solution variables
        double *pp_;
        double *vis_;
        double *rho_;
        double *uvel_;
        double *vvel_;
        double *wvel_;
        double *temp_;
        double *diff_;
        double *hnot_;
        double *fracl_;
        double *fraclnot_;
        double *enthalpy_;
        double *sourceinput_;      
        double *avgconcentration_;
          // geometric variables
        double *volume_;
        double *areaij_;
        double *areajk_;
        double *areaik_;
        double *dxpwinv_;
        double *dypsinv_;
        double *dzpbinv_;
        double *fracx_;
        double *fracy_;
        double *x_;
        double *y_;
        double *z_;
          // flux terms and peak temperature
        double ahtoploss_, heatout_, accul_, heatvol_, ratio_, tpeak_;
        double fluxwest_, fluxeast_, fluxtop_, fluxbottom_, fluxnorth_, fluxsouth_;
          // length, depth, width of melt pool, and steady dimension indicator
        double alen_ = 0.0;
        double depth_ = 0.0; 
        double width_ = 0.0; 
        bool steadylength_ = false;
          // indicies for meltpool size
        int istat_, jstat_, kstat_, iend_, jend_, istatp1_, iendm1_;
          // public laser parameters
        int jstart_;
        double heatinlaser_, beamposx_, beamposy_, beamposz_;
        double oldbeamposx_, oldbeamposy_,oldbeamposz_;
          // laser power index (assume it is initially on)
        int powerindicator_ = 1;
          // indiciator to determine if the simulation is newly intialized
        bool init_ = true;
          // index to track number of laser tracks
        int trackindx_ = 1;

        // Local functions/subroutines
        void initializeBoundaries();

        void updateLaser();

        void getMeltPoolShape(); 

        void applyBCs(int &ivar);

        double getSurfaceTensionGradient(double &t, double &a0);

        void calculateFluxes();
    private: 
        // Private variables
          // private laser parameters
        int istart_, imin_, imax_, jmin_, jmax_, kmin_ ,kmax_;
        double peakhin_;
        double heatin_[_ny_*_nx_];
         
        // Private functions/subroutines
          // function to determine transient laser coordinates
        void getLaserCoords();
          // function to alteratively calculate power law (faster)
        double fastPow(double base, int exponent);
};

#endif


