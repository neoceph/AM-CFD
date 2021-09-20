#ifndef CFDSolverManager_h
#define CFDSolverManager_h

// standard headers

// user-defined headers
#include <Mesh.h>
#include <DomainManager.h>
#include <SolutionVariableManager.h>
#include <BoundCondManager.h>

// define some classes to be used
class DomainManager;
class Mesh;
class SolutionVariableManager;
class BoundCondManager;

class CFDSolverManager
{
    public:
        // Constructor
        CFDSolverManager(Mesh *meshObj,
                         DomainManager *domainMgr,
                         SolutionVariableManager *solVarMgr,
                         BoundCondManager *bcMgr);
        
        // Destructor
        virtual ~CFDSolverManager() {};

        // Input
        Mesh *meshObj_;
        DomainManager *domainMgr_;
        SolutionVariableManager *solVarMgr_;
        BoundCondManager *bcMgr_;

        // Local variables
          // public geometric variables
        double x_[_nx_];
        double y_[_ny_];
        double z_[_nz_];
        int meshXMovementDirection, meshYMovementDirection;
        double dxpwinv_[_nx_];
        double dypsinv_[_ny_];
        double dzpbinv_[_nz_];
        double volume_[_nz_*_ny_*_nx_];
        double areaij_[_ny_*_nx_];
        double areajk_[_nz_*_ny_];
        double areaik_[_nz_*_nx_];
        double fracx_[_nx_];
        double fracy_[_ny_];
          // viscosity matrix 
        double vis_[_nz_*_ny_*_nx_];
          // diffusion matrix
        double diff_[_nz_*_ny_*_nx_];
          // uvw velocities
        double uvel_[_nz_*_ny_*_nx_];
        double vvel_[_nz_*_ny_*_nx_];
        double wvel_[_nz_*_ny_*_nx_];
          // enthalpy, temperature, and previous pressure and enthalpy
        double enthalpy_[_nz_*_ny_*_nx_];
        double temp_[_nz_*_ny_*_nx_];
        double pp_[_nz_*_ny_*_nx_];
        double hnot_[_nz_*_ny_*_nx_];
          // current and previous volume fractions of liquid
        double fracl_[_nz_*_ny_*_nx_];
        double fraclnot_[_nz_*_ny_*_nx_];
          // temperature gradient and solidification growth rate at solidus
        double grad_[_nz_*_ny_*_nx_];
        double rate_[_nz_*_ny_*_nx_];
          // heat input from source
        double sourceinput_[_nz_*_ny_*_nx_];
          // defined activated elements
        //int ***active_; 
          // consolidated fraction for powder
        double csfrac_[_nz_*_ny_*_nx_];
        double rho_[_nz_*_ny_*_nx_];
        double solfrac_[_nz_*_ny_*_nx_];
          // residual error of velocities, mass and enthalpy
	double resoru_, resorv_, resorw_, resorm_, resorh_;
          // number of iterations
        int itertot_, niter_;
          // indicator for initialization
        bool init_ = true;
          // concentration (for species solve) and residual error
        double avgconcentration_[_nz_*_ny_*_nx_];
        double avgconcentrationnot_[_nz_*_ny_*_nx_];
        double concentration_[_nz_*_ny_*_nx_];
        double resorc_;

        // Local functions/subroutines
        void generateGrid();

        void initializeSolver();

        void iterativeLoopSolve();

        void getSolidificationParameters();

        void updateSolutions();
    private: 
          // private geometric variables
        double xu_[_nx_];
        double yv_[_ny_];
        double zw_[_nz_];
        double volume_u_[_nz_*_ny_*_nx_]; 
        double volume_v_[_nz_*_ny_*_nx_];
        double volume_w_[_nz_*_ny_*_nx_];
        double areauij_[_ny_*_nx_];
        double areauik_[_nz_*_nx_]; 
        double areavjk_[_nz_*_ny_]; 
        double areavij_[_ny_*_nx_]; 
        double areawik_[_nz_*_nx_]; 
        double areawjk_[_nz_*_ny_];
        double fracz_[_nz_];
          // previous timestep velocities and pressure
        double unot_[_nz_*_ny_*_nx_];
        double vnot_[_nz_*_ny_*_nx_];
        double wnot_[_nz_*_ny_*_nx_];
        double pressure_[_nz_*_ny_*_nx_];
          // coeffcients for discretized equations
        double ap_[_nz_*_ny_*_nx_];
        double an_[_nz_*_ny_*_nx_];
        double as_[_nz_*_ny_*_nx_];
        double ae_[_nz_*_ny_*_nx_];
        double aw_[_nz_*_ny_*_nx_];
        double at_[_nz_*_ny_*_nx_];
        double ab_[_nz_*_ny_*_nx_];
        double apnot_[_nz_*_ny_*_nx_]; 
          // velocity rise matrix, source term matrix su and sp
        double dux_[_nz_*_ny_*_nx_];
        double dvy_[_nz_*_ny_*_nx_];
        double dwz_[_nz_*_ny_*_nx_];
        double su_[_nz_*_ny_*_nx_];
        double sp_[_nz_*_ny_*_nx_];
          // momentum reference for residuals 
        double refmom_; 
          // some values for sanity checks
        const double small_ = 1.0e-6;
        const double large_ = 1.0e20;
        
        // Private functions/subroutines
          // discretization functions
        double (*discFunct)(double);
        static double discretPower(double pe);
        static double discretUpwind(double pe);
        static double discretExponential(double pe);

        void enhanceConvergenceSpeed(int &ivar);

        void convertEnthalpyToTemp();

        void calculateResiduals(int &ivar);

        void solveLiquidDomain(int &ivar);

        void solveEntireDomain(int &ivar);

        void cleanVelocities();

        void discretizeEquations(int &ivar);

        void calculateSource(int &ivar);

        void correctPressure(int &ivar);

          // material property functions
        void (CFDSolverManager::*updateThermalProps)();
        void bareplateMaterialProps();
        void powderMaterialProps();

        void remapSolutions(const vector <double> &oldx, const vector <double> &oldy,  const vector <double> &oldz);
        double interpolate(double x0, double x1, double y0, double y1, double x);
};

#endif

