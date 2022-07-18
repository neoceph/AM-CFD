#ifndef DomainManager_h
#define DomainManager_h

// standard headers

// user-defined headers
#include <Mesh.h>

// define some classes to be used
class Mesh;

class DomainManager
{
    public:
        // Constructor
        DomainManager(Mesh *meshObj);
        
        // Destructor
        virtual ~DomainManager() {}; 
        
        // Input
        Mesh *meshObj_;

        // Local variables
          // input variables 
        int nzx_, nzy_, nzz_;
        double alaspow_, alasfact_, lasertime_, xstart_, ystart_, scanvel_, apowseta_, heatrb_, avolpow_, avolfact_;
        double alaseta_ = 0.0;
        double alasrb_ = 1.0e-6;
        double layerheight_ = 0.0;
        double heatthick_ = 0.0;
        double dens_, denl_, viscos_, tsolid_, tliquid_, hsmelt_, hlfriz_, acpa_, acpb_, acpl_, thconsa_, thconsb_, thconsc_, 
               thconla_, thconlb_, beta_, emiss_, dgdtp_, darcyKo_;
        double urfu_, urfv_, urfw_, urfp_, urfh_;
        double htci_, htcj_, htck1_, htckn_, tempwest_, tempeast_, tempnorth_, tempbottom_, temppreheat_, a_oxygen_;
          // relevant grid information
        double *xzone_;
        double *ncvx_;
        double *powrx_;
        double *yzone_;
        double *ncvy_;
        double *powry_;
        double *zzone_;
        double *ncvz_;
        double *powrz_;
          // number of nodes for solving (total and internal)
        int ni_, nim1_, nim2_, nj_, njm1_, njm2_, nk_, nkm1_, nkm2_, coordinateIndex[700*300*150][3];
          // current simulation time, latent heat, steady-state logical 
        double timet_, hlcal_, hlatnt_;
        bool steady_;
          // fixed enthalpy for bottom of substrate 
        double enthbottom_;
          // vaporization latent heat, temperature, factor and molar mass
        double evaplatent_, evaptemp_, evapfact_, molarmass_;
          // species initial concentration, diffusion coeff, partition coeff, under-relax factor
        double initconcentration_, massdiff_, kp_, urfc_;
          // energy-based free surface: bi-section interval, surface tension, 
          // max iterations, tolerances, under-relax factor, and solver output frequency
        double initlam1_, initlam2_, surftens_, maxiterenergy_, tolenergy_, maxiterbisect_, 
               tolbisect_, urfree_; 
          // material properties for powder-bed process
        double powderden_, powdercpa_, powdercpb_, powderthcona_, powderthconb_;
        int numlayers_, ncvzlayer_;
          // toolpath vectors for GAMMA formated file
        vector< vector<double>> tooltxyz_;
        vector<int> laserOn_, trackID;
        // vector for reading heat source parameter as an input
        vector< vector<double>> heatSourceParameterValues;
          // outputting a section of the domain: box dimensions and output time
        double gxmin_, gymin_, gzmin_, gxmax_, gymax_, gzmax_, toutsect_;
          // tracker to know when simulation is finished
        bool finish_ = false;
        
        // Local functions/subroutines
        void assignInputValues();

        void getToolpath();
        
        void getHeatSourceParameters();
};

        // global maximum node numbers (for performance)
        const int _nx_ = 700;
        const int _ny_ = 300;
        const int _nz_ = 300;

#endif

