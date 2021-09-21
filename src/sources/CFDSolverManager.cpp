// standard headers
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string.h>

// user-defined headers
#include <CFDSolverManager.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
CFDSolverManager::CFDSolverManager(Mesh *meshObj,
                                   DomainManager *domainMgr,
                                   SolutionVariableManager *solVarMgr,
                                   BoundCondManager *bcMgr)
                                   : domainMgr_(domainMgr),
                                     meshObj_(meshObj),
                                     solVarMgr_(solVarMgr),
                                     bcMgr_(bcMgr)
{
}

//////////////////////////////////////////////////////
//		generateGrid()	                    //
//////////////////////////////////////////////////////
void
CFDSolverManager::generateGrid()
{
    // grab some stuff from the domain
    int nzx = domainMgr_->nzx_; 
    double *xzone = domainMgr_->xzone_;
    double *ncvx = domainMgr_->ncvx_;
    double *powrx = domainMgr_->powrx_;
    int nzy = domainMgr_->nzy_;
    double *yzone = domainMgr_->yzone_;
    double *ncvy = domainMgr_->ncvy_;
    double *powry = domainMgr_->powry_;
    int nzz = domainMgr_->nzz_;
    double *zzone = domainMgr_->zzone_;
    double *ncvz = domainMgr_->ncvz_;
    double *powrz = domainMgr_->powrz_;
    
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;

    // declare some local variables
    int i, j, k, ist;
    double statloc, term;

    // copy the old x coordinates (for moving mesh)
    std::vector <double> oldx(x_, x_+ni);
    std::vector <double> oldy(y_, y_+nj);
    std::vector <double> oldz(z_, z_+nk);
    
    //declare strings for test
    // string timeFileOutput, filename;


    
        // update the x-zone
    if(!init_)
    {
        // update the mesh zone
        double dx = bcMgr_->beamposx_ - bcMgr_->oldbeamposx_;
        double dy = bcMgr_->beamposy_ - bcMgr_->oldbeamposy_;

        // Code block for determining the X-direction movements
        if (dx < 0)
        {
            meshXMovementDirection = -1;
        }
        else
        {
            meshXMovementDirection = 1;
        }


        // Code block for determining the Y-direction movement
        if (dy < 0)
        {
            meshYMovementDirection = -1;
        }
        else
        {
            meshYMovementDirection = 1;
        }
        
        int indx = meshObj_->xmoveindex_;
        int indy = meshObj_->ymoveindex_;
        if(bcMgr_->steadylength_)
        {
            xzone[indx-1] += dx;
            xzone[indx+1] -= dx;

            yzone[indy-1] += dy;
            yzone[indy+1] -= dy;
        }
        else
        {
            // xzone[indx] += dx;
            // xzone[indx+1] -= dx;
            xzone[indx-1] += dx;
            xzone[indx+1] -= dx;

            yzone[indx-1] += dy;
            yzone[indx+1] -= dy;
        }//end if

        // update the number of control volumes in "outer" zones
        // ncvx[0] = floor((ni - 2.0 - ncvx[1])/(1.0 + xzone[2]/xzone[0]));
        // ncvx[2] = ni - 2.0 - ncvx[1] - ncvx[0];
    }//end if
    
    // create the multi-dimensional arrays
    // x-grid-----------------------------------------
    fill(&xu_[0], &xu_[2], 0.0);
    ist = 2;
    statloc = 0.0;

    // loop through "x" sub-regions
    for(i=0; i<nzx; i++)
    {
        for(j=0; j<ncvx[i]; j++)
        {
            if(powrx[i] >= 0.0)  
            {
                // grids transit from fine to coarse
                term = pow((double)(j+1)/ncvx[i], powrx[i]);
            }
            else
            {
                // grids transit from coarse to fine
                term = 1.0 - pow(1.0-(double)(j+1)/ncvx[i], -powrx[i]);
            }//end if
            xu_[j+ist] = statloc + xzone[i]*term;
        }//end for(j)
        ist += (int)ncvx[i];
        statloc += xzone[i];
    }//end for(i)

    // central interpolation
    for (i=0; i<nim1; i++)
        x_[i] = (xu_[i+1] + xu_[i])*0.5;

    // last coordinate value of scalar node equals to that of uVel nodes
    x_[nim1] = xu_[nim1];
     
    // y-grid-----------------------------------------
    fill(&yv_[0], &yv_[2], 0.0);
    ist = 2;
    statloc = 0.0;

    // loop through "y" sub-regions
    for (i=0; i<nzy; i++)
    {
        for (j=0; j<ncvy[i]; j++)
        {
            if(powry[i] >= 0.0)
            {
                term = pow((double)(j+1) / ncvy[i], powry[i]);
            }
            else
            {
                term = 1.0 - pow(1.0-(double)(j+1) / ncvy[i], -powry[i]);
            }//end if
            yv_[j+ist] = statloc + yzone[i]*term;
        }//end for(j)
        ist += (int)ncvy[i];
        statloc += yzone[i];
    }//end for(i)
    
    // central interpolation
    for (i=0; i<njm1; i++)
        y_[i] = (yv_[i+1] + yv_[i])*0.5;

    // last coordinate value of scalar node equals to that of vVel nodes
    y_[njm1] = yv_[njm1];

    // z-grid-----------------------------------------
    fill(&zw_[0], &zw_[2], 0.0);
    ist = 2;
    statloc = 0.0;

    // loop through "z" sub-regions
    for (i=0; i<nzz; i++)
    {
        for (j=0; j < ncvz[i]; j++)
        {
            if(powrz[i] >= 0.0)
            {
                term = pow((double)(j+1) / ncvz[i], powrz[i]);
            }
            else
            {
                term = 1.0 - pow(1.0-(double)(j+1) / ncvz[i], -powrz[i]);
            }//end if
            zw_[j+ist] = statloc + zzone[i]*term;
        }//end for(j)
        ist += (int)ncvz[i];
        statloc += zzone[i];
    }//end for(i)

    // central interpolation
    for (i=0; i<nkm1; i++)
        z_[i] = (zw_[i+1] + zw_[i])*0.5;

    // last coordinate value of scalar node equals to that of wVel nodes
    z_[nkm1] = zw_[nkm1];

    // loop over velocity nodes and calculate reciprical of node distances
    for (i=1; i<ni; i++)    
        dxpwinv_[i] = 1.0/(x_[i] - x_[i-1]);

    for (j=1; j<nj; j++)
        dypsinv_[j] = 1.0/(y_[j] - y_[j-1]);

    for (k=1; k<nk; k++)
        dzpbinv_[k] = 1.0/(z_[k] - z_[k-1]);

    // interpolate and calculate distance between nodes and adjacent interfaces
    for (i=0; i<nim1; i++)
        fracx_[i] = (x_[i+1] - xu_[i+1])/(x_[i+1] - x_[i]);

    for (j=0; j<njm1; j++)
        fracy_[j] = (y_[j+1] - yv_[j+1])/(y_[j+1] - y_[j]);

    for (k=0; k<nkm1; k++)
        fracz_[k] = (z_[k+1] - zw_[k+1])/(z_[k+1] - z_[k]);
    
    // calculate volumes of CV
    for (k=1; k<nkm1; k++)
    {
        for (j=1; j<njm1; j++)
        {
            for (i=1; i<nim1; i++)
            {
                volume_[(k*nj*ni)+(j*ni)+i] = (xu_[i+1]-xu_[i]) * (yv_[j+1]-yv_[j]) * (zw_[k+1]-zw_[k]);
                volume_u_[(k*nj*ni)+(j*ni)+i] = (x_[i]-x_[i-1]) * (yv_[j+1]-yv_[j]) * (zw_[k+1]-zw_[k]);
                volume_v_[(k*nj*ni)+(j*ni)+i] = (xu_[i+1]-xu_[i]) * (y_[j]-y_[j-1]) * (zw_[k+1]-zw_[k]);
                volume_w_[(k*nj*ni)+(j*ni)+i] = (xu_[i+1]-xu_[i]) * (yv_[j+1]-yv_[j]) * (z_[k]-z_[k-1]);
            }//end for(i)
        }//end for(j)
    }//end for(k)

    // calculate areas in three directions of CV
    for (j=1; j<njm1; j++)
    {
        for (i=1; i<nim1; i++)
        {
            areaij_[(j*ni)+i] = (xu_[i+1]-xu_[i]) * (yv_[j+1]-yv_[j]);
            areauij_[(j*ni)+i] = (x_[i]-x_[i-1]) * (yv_[j+1]-yv_[j]);
            areavij_[(j*ni)+i] = (xu_[i+1]-xu_[i]) * (y_[j]-y_[j-1]);
        }//end for(i)
    }//end for(j)
    
    for (k=1; k<nkm1; k++)
    {
        for (i=1; i<nim1; i++)
        {
            areaik_[(k*ni)+i] = (xu_[i+1]-xu_[i]) * (zw_[k+1]-zw_[k]);
            areawik_[(k*ni)+i] = (xu_[i+1]-xu_[i]) * (z_[k]-z_[k-1]);
            areauik_[(k*ni)+i] = (x_[i]-x_[i-1]) * (zw_[k+1]-zw_[k]);
        }//end for(i)
    }//end for(k)

    for (k=1; k<nkm1; k++)
    {
        for (j=1; j<njm1; j++)
        {
            areajk_[(k*nj)+j] = (yv_[j+1]-yv_[j]) * (zw_[k+1]-zw_[k]);
            areavjk_[(k*nj)+j] = (y_[j]-y_[j-1]) * (zw_[k+1]-zw_[k]);
            areawjk_[(k*nj)+j] = (yv_[j+1]-yv_[j]) * (z_[k]-z_[k-1]);
        }//end for(j)
    }//end for(k)
   
    // map geometric variables to the viewing window
    solVarMgr_->mapGeometryVars_["volume"] = &volume_[0];       
    solVarMgr_->mapGeometryVars_["areaij"] = &areaij_[0];       
    solVarMgr_->mapGeometryVars_["areajk"] = &areajk_[0];       
    solVarMgr_->mapGeometryVars_["areaik"] = &areaik_[0];       
    solVarMgr_->mapGeometryVars_["dxpwinv"] = &dxpwinv_[0];       
    solVarMgr_->mapGeometryVars_["dypsinv"] = &dypsinv_[0];       
    solVarMgr_->mapGeometryVars_["dzpbinv"] = &dzpbinv_[0];       
    solVarMgr_->mapGeometryVars_["fracx"] = &fracx_[0];       
    solVarMgr_->mapGeometryVars_["fracy"] = &fracy_[0];       
    solVarMgr_->mapGeometryVars_["xu"] = &xu_[0];       
    solVarMgr_->mapGeometryVars_["yv"] = &yv_[0];       
    solVarMgr_->mapGeometryVars_["x"] = &x_[0];       
    solVarMgr_->mapGeometryVars_["y"] = &y_[0];       
    solVarMgr_->mapGeometryVars_["z"] = &z_[0];       





// // writing variable to disc for further analysis
//     timeFileOutput = "time="+to_string(this->domainMgr_->timet_);
//     filename = "oldTemp"+timeFileOutput+".txt";

//     ofstream myoldfile (filename);
//     if (myoldfile.is_open())
//     {
//         for (int loc = 0; loc<ni*nj*nk; loc++)
//             myoldfile << temp_[loc] <<"\n";      
//     }
    
//     myoldfile.close();
    
//     filename = "oldX"+timeFileOutput+".txt";

//     ofstream myoldXfile (filename);
//     if (myoldXfile.is_open())
//     {
//         for (int i = 0; i<ni; i++)
//             myoldXfile << oldx[i] <<"\n";      
//     }

//     myoldXfile.close();

//     filename = "oldY"+timeFileOutput+".txt";

//     ofstream myoldYfile (filename);
//     if (myoldYfile.is_open())
//     {
//         for (int j = 0; j<nj; j++)
//             myoldYfile << oldy[j] <<"\n";      
//     }

//     myoldYfile.close();

// end of first stage of file writing.

    if(!init_)
        remapSolutions(oldx, oldy, oldz);
    
 
// //second stage of file writing. 
//     filename = "newTemp"+timeFileOutput+".txt";
    
//     ofstream mynewfile (filename);
//     if (mynewfile.is_open())
//     {
//         for (int loc = 0; loc<ni*nj*nk; loc++)
//             mynewfile << temp_[loc] <<"\n";
//     }

//     mynewfile.close();

//     filename = "X_"+timeFileOutput+".txt";

//     ofstream myX_file (filename);
//     if (myX_file.is_open())
//     {
//         for (int i = 0; i<ni; i++)
//             myX_file << x_[i] <<"\n";      
//     }

//     myX_file.close();

//     filename = "Y_"+timeFileOutput+".txt";

//     ofstream myY_file (filename);
//     if (myY_file.is_open())
//     {
//         for (int j = 0; j<nj; j++)
//             myY_file << y_[j] <<"\n";      
//     }

//     myY_file.close();


//     filename = "Z_"+timeFileOutput+".txt";

//     ofstream myZ_file (filename);
//     if (myZ_file.is_open())
//     {
//         for (int k = 0; k<nk; k++)
//             myZ_file << z_[k] <<"\n";      
//     }

//     myZ_file.close();    
// //end of second stage of file writing.


}//end generateGrid

//////////////////////////////////////////////////////
//		initializeSolver                    //
//////////////////////////////////////////////////////
void 
CFDSolverManager::initializeSolver()
{
    // grab some stuff from domain Mgr
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    double viscos = domainMgr_->viscos_;
    double dens = domainMgr_->dens_;
    double acpa = domainMgr_->acpa_;
    double acpb = domainMgr_->acpb_;
    double thconsa = domainMgr_->thconsa_;
    double thconsb = domainMgr_->thconsb_;
    double thconsc = domainMgr_->thconsc_;
    double temppreheat = domainMgr_->temppreheat_;
    double tempwest = domainMgr_->tempwest_;
    double tempeast = domainMgr_->tempeast_;
    double tempnorth = domainMgr_->tempnorth_;
    double tempbottom = domainMgr_->tempbottom_;
    
    // define some local variables
    double enthalpyWest, enthalpyEast, enthalpyNorth, enthalpyBottom, enthalpyPreheat;
    int i, j, k;

    // translate preheat temperature to enthalpy at the boundaries
    enthalpyPreheat = 0.5*acpa*temppreheat*temppreheat + acpb*temppreheat; 
    enthalpyWest = 0.5*acpa*tempwest*tempwest + acpb*tempwest;            
    enthalpyEast = 0.5*acpa*tempeast*tempeast + acpb*tempeast;            
    enthalpyNorth = 0.5*acpa*tempnorth*tempnorth + acpb*tempnorth;        
    enthalpyBottom = 0.5*acpa*tempbottom*tempbottom + acpb*tempbottom;    

    // intialize values for all the arrays
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        vis_[loc] = viscos;
        rho_[loc] = dens;
        diff_[loc] = (thconsa*temppreheat*temppreheat + thconsb*temppreheat + thconsc)/
                        (acpa*temppreheat+acpb); 
        enthalpy_[loc] = enthalpyPreheat;
        hnot_[loc] = enthalpyPreheat;
        temp_[loc] = temppreheat;
        dux_[loc] = 0.0;
        dvy_[loc] = 0.0;
        dwz_[loc] = 0.0;
        fracl_[loc] = 0.0;
        fraclnot_[loc] = 0.0;
        sourceinput_[loc] = 0.0;
        su_[loc] = 0.0;
        sp_[loc] = 0.0;
        ap_[loc] = 0.0;
        an_[loc] = 0.0;
        as_[loc] = 0.0;
        ae_[loc] = 0.0;
        aw_[loc] = 0.0;
        at_[loc] = 0.0;
        ab_[loc] = 0.0;
        apnot_[loc] = 0.0; 
        uvel_[loc] = 0.0;
        vvel_[loc] = 0.0;
        wvel_[loc] = 0.0;
        unot_[loc] = 0.0;
        vnot_[loc] = 0.0;
        wnot_[loc] = 0.0;
        pressure_[loc] = 0.0;
        pp_[loc] = 0.0;
        if(meshObj_->solidificationparams_)
        {
            grad_[loc] = 0.0;
            rate_[loc] = 0.0;
        }
        if(meshObj_->species_)
        {
            avgconcentration_[loc] = domainMgr_->initconcentration_;
            avgconcentrationnot_[loc] = domainMgr_->initconcentration_;
            concentration_[loc] = domainMgr_->initconcentration_;
        }

        // added here for powder simulations
        if(meshObj_->powderbed_)
        {
            solfrac_[loc] = 0.0;
            int k = loc/(ni*nj);
            if(z_[k] < z_[nk-1] - ((double)domainMgr_->numlayers_*domainMgr_->layerheight_))
                csfrac_[loc] = 1.0;
            else
                csfrac_[loc] = 0.0;
        }//end if
        //active_[loc] = 1;
    }//end for(i)

    //----enthalpy BC: initial enthalpy at i=0/i=nim1 planes
    for (k=0; k<nk; k++)
    {
        for (j=0; j<nj; j++)
        {
            enthalpy_[(k*nj*ni)+(j*ni)+0] = enthalpyWest;
            enthalpy_[(k*nj*ni)+(j*ni)+nim1] = enthalpyEast;
        }//end for(j)
    }//end for(k)

    //----k=0 plane----------------
    for (j=0; j<nj; j++)
    {
        for (i=0; i<ni; i++)
        {
            enthalpy_[(0*nj*ni)+(j*ni)+i] = enthalpyBottom;
        }//end for(i)
    }//end for(j)

    //----j=njm1 plane---------------
    for (k=0; k<nk; k++)
    {
        for (i=0; i<ni; i++)
        {
            enthalpy_[(k*nj*ni)+(njm1*ni)+i] = enthalpyNorth;
        }//end for(i)
    }//end for(k)

    // set fixed enthalpy
    domainMgr_->enthbottom_ = enthalpyBottom;
     
    // map solution variables to the viewing window
    solVarMgr_->mapSolutionVars_["pp"] = &pp_[0];       
    solVarMgr_->mapSolutionVars_["vis"] = &vis_[0];      
    solVarMgr_->mapSolutionVars_["rho"] = &rho_[0];      
    solVarMgr_->mapSolutionVars_["uvel"] = &uvel_[0];     
    solVarMgr_->mapSolutionVars_["vvel"] = &vvel_[0];     
    solVarMgr_->mapSolutionVars_["wvel"] = &wvel_[0];     
    solVarMgr_->mapSolutionVars_["temp"] = &temp_[0];     
    solVarMgr_->mapSolutionVars_["diff"] = &diff_[0];     
    solVarMgr_->mapSolutionVars_["hnot"] = &hnot_[0];     
    solVarMgr_->mapSolutionVars_["fracl"] = &fracl_[0];   
    solVarMgr_->mapSolutionVars_["fraclnot"] = &fraclnot_[0];   
    solVarMgr_->mapSolutionVars_["pressure"] = &pressure_[0];
    solVarMgr_->mapSolutionVars_["enthalpy"] = &enthalpy_[0]; 
    solVarMgr_->mapSolutionVars_["sourceinput"] = &sourceinput_[0];
    if(meshObj_->solidificationparams_)
    {
        solVarMgr_->mapSolutionVars_["rate"] = &rate_[0];     
        solVarMgr_->mapSolutionVars_["grad"] = &grad_[0];     
    }//end if
    if(meshObj_->species_)
    {
        solVarMgr_->mapSolutionVars_["avgconcentration"] = &avgconcentration_[0];
        solVarMgr_->mapSolutionVars_["concentration"] = &concentration_[0];
    }//end if

    // determine if the material is a power or bulk material
    if(meshObj_->powderbed_)
    {
        updateThermalProps = &CFDSolverManager::powderMaterialProps;
        solVarMgr_->mapSolutionVars_["csfrac"] = &csfrac_[0];   
        solVarMgr_->mapSolutionVars_["solfrac"] = &solfrac_[0];   
    }
    else
    {
        updateThermalProps = &CFDSolverManager::bareplateMaterialProps;
    }//end if

    // determine the type of discretization that should be used
    if(meshObj_->upwind_) 
    {
        discFunct = discretUpwind;
    }
    else if(meshObj_->powerlaw_) 
    {
        discFunct = discretPower;
    }
    else if(meshObj_->exponential_) 
    {
        discFunct = discretExponential;
    }//end if

    // zero out some variables 
    resoru_ = 0.0;
    resorv_ = 0.0;
    resorw_ = 0.0;
    resorm_ = 0.0;
    resorh_ = 0.0;
    resorc_ = 0.0;
    itertot_ = 0; 
}//end initializeSolver

//////////////////////////////////////////////////////
//		iterativeLoopSolve                  //
//////////////////////////////////////////////////////
void 
CFDSolverManager::iterativeLoopSolve()
{
    // grab some stuff from the mesh
    double tol = meshObj_->nonlintol_;

    // define some local variables 
    int ivar; 
    double ratio;
    double amaxres = 0.0;
    double alenold = bcMgr_->alen_; 
    // string filename, timeFileOutput;

    // intialize the number of iterations to zero
    niter_ = 0;

    // itertatively solve for all desired solution variables
    while (niter_ < meshObj_->maxit_)
    { 
        // increment the iterators
        niter_++;
        itertot_++;
        
        // first solve energy equation
        ivar = 5;                        // ivar = 5 => solve the energy equation
        bcMgr_->applyBCs(ivar);          // thermal bondary condition
        discretizeEquations(ivar);       // calc coeffs of discretized enthalpy conservation equation
        calculateSource(ivar);           // calc source term coefficients to complete discretization
        calculateResiduals(ivar);        // calc residual error in the domain
        (this->*updateThermalProps)();   // update temperature-dependent thermo-properties

        enhanceConvergenceSpeed(ivar);   // calc residual error in each x direction slice to enhance converge speed (TDMA)

        // timeFileOutput = "time="+to_string(this->domainMgr_->timet_);
        // filename = "Temp_before_solve"+timeFileOutput+".txt";
        

        // ofstream myoldfile (filename);
        // if (myoldfile.is_open())
        // {
        //     for (int loc = 0; loc<this->domainMgr_->ni_*domainMgr_->nj_*domainMgr_->nk_; loc++)
        //         myoldfile << temp_[loc] <<"\n";      
        // }
    
        // myoldfile.close();


        solveEntireDomain(ivar);         // solve enthalpy conservation equation (Line-by-line TDMA)
        convertEnthalpyToTemp();         // translate from enthalpy to temperature

        // filename = "Temp_after_solve"+timeFileOutput+".txt";

        // ofstream mynewfile (filename);
        // if (mynewfile.is_open())
        // {
        //     for (int loc = 0; loc<this->domainMgr_->ni_*domainMgr_->nj_*domainMgr_->nk_; loc++)
        //         mynewfile << temp_[loc] <<"\n";      
        // }
    
        // mynewfile.close();


        bcMgr_->getMeltPoolShape();      // get melt pool dimension, start and end index of i,j,k to detemine fluid region

        // if peak temperature less than solidus, do not need to solve  momentum eqution
        if(bcMgr_->tpeak_ > domainMgr_->tsolid_ && meshObj_->navierstokes_)
        {
            // zero out velocites outside the melt pool
            cleanVelocities();

            // solve for for momentum and mass
            for (ivar = 1; ivar < 5; ivar++)
            {
                 bcMgr_->applyBCs(ivar);
                 discretizeEquations(ivar);
                 calculateSource(ivar);
                 calculateResiduals(ivar);
                 solveLiquidDomain(ivar);
                 correctPressure(ivar);
            }//end for(ivar)
        }//end if
        
        // solve for species conservation (if desired)
        if(meshObj_->species_)
        {
            ivar = 6;
            bcMgr_->applyBCs(ivar);
            discretizeEquations(ivar);
            calculateSource(ivar);
            calculateResiduals(ivar);
            enhanceConvergenceSpeed(ivar);   
            solveEntireDomain(ivar);         
        }//end if
        
        // calculate maximum residual
        amaxres = max({resorm_,resoru_,resorv_,resorw_,resorc_});
         
        // calculate fluxes
        bcMgr_->calculateFluxes();
        ratio = bcMgr_->ratio_;

        // check convergence of all equations
        if(domainMgr_->steady_)
        {
            if(ratio <= 1.001 && ratio >= 0.999 && amaxres < tol) 
                return;
        }      
        else
        {
            // transient-state criteria (heating stage) 
            if(bcMgr_->powerindicator_ == 1)
            {
               if(ratio <= 1.01 && ratio >= 0.99 && amaxres < tol) 
                   return;  
            }
            // transient-state criteria (cooling stage)
            else
            {
               if(resorh_ < 1e-6 && amaxres < tol) 
                   return;
            }//end if
        }//end if
    }//end while(niter)
    
    // check if the meltpool has steady state dimensions (for moving mesh)
    if(abs(bcMgr_->alen_ - alenold) < 1.0e-6)
        bcMgr_->steadylength_ = true;
}//end iterativeLoopSolve

//////////////////////////////////////////////////////
//		 discretizeEquations                //
//////////////////////////////////////////////////////
void 
CFDSolverManager::discretizeEquations(int &ivar)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    int nim2 = domainMgr_->nim2_;
    int njm2 = domainMgr_->njm2_;
    int nkm2 = domainMgr_->nkm2_;
    
    // grab some stuff from the mesh
    double delt = meshObj_->delt_;
    
    // grab some stuff from the bc
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;
     
    // define some local variables
    int i, j, k;
    double vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb;
    double visn, viss, visv, visvw, visve;
    double visu, visun, visus;
    double vise, visw, visut, visub;
    double vist, visb;
    
    double visvt, visvb;
    double ds, dn, de, dw, dt, db;
    double delf, cp0, cp1;
    double dudxp, dudxm, dvdxp, dvdxm, dwdxp, dwdxm;
    double dudyp, dudym, dvdyp, dvdym, dwdyp, dwdym;
    double vis_w, viswn, visws, viswe, visww;
    double dudzp, dudzm, dvdzp, dvdzm, dwdzp, dwdzm;
    double rhobydt, difn, dife, dift;

    // discretize the respective conservation equation using FVM 
    switch(ivar)
    {
        case 1:
        {
            #pragma omp parallel for \
                    private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, visu, visun, visus, visn, viss, vise, visw, visut, visub, vist, visb, dn, ds, dw, de, dt, db, delf, cp0, cp1, dudxp, dudxm, dvdxp, dvdxm, dwdxp, dwdxm) \
                    collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        //----velocities at cv faces----------------------------
                        vn = vvel_[(k*nj*ni)+((j+1)*ni)+i]*(1.0 - fracx_[i-1]) + vvel_[(k*nj*ni)+((j+1)*ni)+i-1]*fracx_[i-1];
                        vs = vvel_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracx_[i-1]) + vvel_[(k*nj*ni)+(j*ni)+i-1]*fracx_[i-1];
                        ue = (uvel_[(k*nj*ni)+(j*ni)+i+1] + uvel_[(k*nj*ni)+(j*ni)+i])*0.5;
                        uw = (uvel_[(k*nj*ni)+(j*ni)+i-1] + uvel_[(k*nj*ni)+(j*ni)+i])*0.5;
                        wt = wvel_[((k+1)*nj*ni)+(j*ni)+i]*(1.0-fracx_[i-1]) + wvel_[((k+1)*nj*ni)+(j*ni)+i-1]*fracx_[i-1];
                        wb = wvel_[(k*nj*ni)+(j*ni)+i]*(1.0-fracx_[i-1]) + wvel_[(k*nj*ni)+(j*ni)+i-1]*fracx_[i-1];
                        
                        //----convection coefficients---------------------------
                        fn = vn*rho_[(k*nj*ni)+(j*ni)+i]*areauik_[(k*ni)+i];
                        fs = vs*rho_[(k*nj*ni)+(j*ni)+i]*areauik_[(k*ni)+i];
                        fe = ue*rho_[(k*nj*ni)+(j*ni)+i]*areajk_[(k*nj)+j];
                        fw = uw*rho_[(k*nj*ni)+(j*ni)+i]*areajk_[(k*nj)+j];
                        ft = wt*rho_[(k*nj*ni)+(j*ni)+i]*areauij_[(j*ni)+i];
                        fb = wb*rho_[(k*nj*ni)+(j*ni)+i]*areauij_[(j*ni)+i];
                        
                        //----viscosity at cv faces-----------------------------
                        visu = vis_[(k*nj*ni)+(j*ni)+i]*vis_[(k*nj*ni)+(j*ni)+i-1]/(fracx_[i-1]*vis_[(k*nj*ni)+(j*ni)+i-1] + (1.0-fracx_[i-1])*vis_[(k*nj*ni)+(j*ni)+i]);
                        visun = vis_[(k*nj*ni)+((j+1)*ni)+i]*vis_[(k*nj*ni)+((j+1)*ni)+i-1]/(fracx_[i-1]*vis_[(k*nj*ni)+((j+1)*ni)+i-1] + (1.0-fracx_[i-1])*vis_[(k*nj*ni)+((j+1)*ni)+i]);
                        visus = vis_[(k*nj*ni)+((j-1)*ni)+i]*vis_[(k*nj*ni)+((j-1)*ni)+i-1]/(fracx_[i-1]*vis_[(k*nj*ni)+((j-1)*ni)+i-1] + (1.0-fracx_[i-1])*vis_[(k*nj*ni)+((j-1)*ni)+i]);
                        if(j == njm2)
                        {
                            visn = visu;
                        }
                        else 
                        {
                            visn = visu*visun/((1.0 - fracy_[j])*visun + fracy_[j]*visu);
                        }//end if

                        if(j == 1)  
                        {
                            viss = visu;
                        }
                        else
                        {
                            viss = visu*visus/(fracy_[j-1]*visus + (1.0 - fracy_[j-1])*visu);
                        }//end if
                        
                        vise = vis_[(k*nj*ni)+(j*ni)+i];
                        visw = vis_[(k*nj*ni)+(j*ni)+i-1];
                        visut = vis_[((k+1)*nj*ni)+(j*ni)+i]*vis_[((k+1)*nj*ni)+(j*ni)+i-1]/(fracx_[i-1]*vis_[((k+1)*nj*ni)+(j*ni)+i-1] +
                                    (1.0-fracx_[i-1])*vis_[((k+1)*nj*ni)+(j*ni)+i]);
                        visub = vis_[((k-1)*nj*ni)+(j*ni)+i]*vis_[((k-1)*nj*ni)+(j*ni)+i-1]/(fracx_[i-1]*vis_[((k-1)*nj*ni)+(j*ni)+i-1] + 
                                    (1.0-fracx_[i-1])*vis_[((k-1)*nj*ni)+(j*ni)+i]);
                        
                        if(k == nkm2)
                        {
                            vist = visu;
                        }
                        else
                        {
                            vist = visut*visu/(fracz_[k]*visu + (1.0 - fracz_[k])*visut);
                        }//end if

                        if(k == 1)
                        {
                            visb = visu;
                        }
                        else
                        {
                            visb = visub*visu/((1.0 - fracz_[k-1])*visu + fracz_[k-1]*visub);
                        }//end if
                        
                        //----diffusion coefficients----------------------------
                        dn = visn*areauik_[(k*ni)+i]*dypsinv_[j+1];
                        ds = viss*areauik_[(k*ni)+i]*dypsinv_[j];
                        de = vise*areajk_[(k*nj)+j]/(xu_[i+1]-xu_[i]);
                        dw = visw*areajk_[(k*nj)+j]/(xu_[i]-xu_[i-1]);
                        dt = vist*areauij_[(j*ni)+i]*dzpbinv_[k+1];
                        db = visb*areauij_[(j*ni)+i]*dzpbinv_[k];
                        
                        //----matrix coefficients-------------------------------
                        an_[(k*nj*ni)+(j*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                        as_[(k*nj*ni)+(j*ni)+i] = ds*discFunct(abs(fs)/ds) + max(0.0,fs);
                        ae_[(k*nj*ni)+(j*ni)+i] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                        aw_[(k*nj*ni)+(j*ni)+i] = dw*discFunct(abs(fw)/dw) + max(0.0,fw);
                        at_[(k*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                        ab_[(k*nj*ni)+(j*ni)+i] = db*discFunct(abs(fb)/db) + max(0.0,fb);
                        apnot_[(k*nj*ni)+(j*ni)+i] = rho_[(k*nj*ni)+(j*ni)+i]*volume_u_[(k*nj*ni)+(j*ni)+i]/delt;
                        
                        //----su and sp--------------------------------
                        delf = fn - fs + fe - fw + ft - fb;
                        cp0 = max(0.0,delf);
                        cp1 = min(0.0,delf);
                        su_[(k*nj*ni)+(j*ni)+i] = -cp1*uvel_[(k*nj*ni)+(j*ni)+i];
                        su_[(k*nj*ni)+(j*ni)+i] += areajk_[(k*nj)+j]*(pressure_[(k*nj*ni)+(j*ni)+i-1] - pressure_[(k*nj*ni)+(j*ni)+i]);
                        sp_[(k*nj*ni)+(j*ni)+i] = -cp0;
                        su_[(k*nj*ni)+(j*ni)+i] += apnot_[(k*nj*ni)+(j*ni)+i]*unot_[(k*nj*ni)+(j*ni)+i];
                        
                        dudxp = (uvel_[(k*nj*ni)+(j*ni)+i+1] - uvel_[(k*nj*ni)+(j*ni)+i])/(xu_[i+1] - xu_[i]);
                        dudxm = (uvel_[(k*nj*ni)+(j*ni)+i] - uvel_[(k*nj*ni)+(j*ni)+i-1])/(xu_[i] - xu_[i-1]);
                        su_[(k*nj*ni)+(j*ni)+i] += (vise*dudxp-visw*dudxm)*areajk_[(k*nj)+j];
                        dvdxp = (vvel_[(k*nj*ni)+((j+1)*ni)+i] - vvel_[(k*nj*ni)+((j+1)*ni)+i-1])*dxpwinv_[i];
                        dvdxm = (vvel_[(k*nj*ni)+(j*ni)+i] - vvel_[(k*nj*ni)+(j*ni)+i-1])*dxpwinv_[i];
                        su_[(k*nj*ni)+(j*ni)+i] += (visn*dvdxp - viss*dvdxm)*areauik_[(k*ni)+i];
                        dwdxp = (wvel_[((k+1)*nj*ni)+(j*ni)+i] - wvel_[((k+1)*nj*ni)+(j*ni)+i-1])*dxpwinv_[i];
                        dwdxm = (wvel_[(k*nj*ni)+(j*ni)+i] - wvel_[(k*nj*ni)+(j*ni)+i-1])*dxpwinv_[i];
                        su_[(k*nj*ni)+(j*ni)+i] += (vist*dwdxp - visb*dwdxm)*areauij_[(j*ni)+i];
                    }//end for(i)
                }//end (j)
            }// end for(k)
            break;
        }//end case(1)
        case 2:
        {
            #pragma omp parallel for \
                    private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, visn, viss, visv, visve, visvw, vise, visw, visvt, visvb, vist, visb, dn, ds, dw, de, dt, db, delf, cp0, cp1, dudyp, dudym, dvdyp, dvdym, dwdyp, dwdym) \
                    collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        vn = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j+1)*ni)+i])*0.5;
                        vs = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j-1)*ni)+i])*0.5;
                        ue = uvel_[(k*nj*ni)+(j*ni)+i+1]*(1.0 - fracy_[j-1]) + uvel_[(k*nj*ni)+((j-1)*ni)+i+1]*fracy_[j-1];
                        uw = uvel_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracy_[j-1]) + uvel_[(k*nj*ni)+((j-1)*ni)+i]*fracy_[j-1];
                        wt = wvel_[((k+1)*nj*ni)+(j*ni)+i]*(1.0 - fracy_[j-1]) + wvel_[((k+1)*nj*ni)+((j-1)*ni)+i]*fracy_[j-1];
                        wb = wvel_[(k*nj*ni)+(j*ni)+i]*(1.0- fracy_[j-1]) + wvel_[(k*nj*ni)+((j-1)*ni)+i]*fracy_[j-1];

                        //----convection coefficients---------------------------
                        fn = vn*rho_[(k*nj*ni)+(j*ni)+i]*areaik_[(k*ni)+i];
                        fs = vs*rho_[(k*nj*ni)+(j*ni)+i]*areaik_[(k*ni)+i];
                        fe = ue*rho_[(k*nj*ni)+(j*ni)+i]*areavjk_[(k*nj)+j];
                        fw = uw*rho_[(k*nj*ni)+(j*ni)+i]*areavjk_[(k*nj)+j];
                        ft = wt*rho_[(k*nj*ni)+(j*ni)+i]*areavij_[(j*ni)+i];
                        fb = wb*rho_[(k*nj*ni)+(j*ni)+i]*areavij_[(j*ni)+i];
                        
                        //----convection coefficients---------------------------
                        visn = vis_[(k*nj*ni)+(j*ni)+i];
                        viss = vis_[(k*nj*ni)+((j-1)*ni)+i];
                        visv = vis_[(k*nj*ni)+(j*ni)+i]*vis_[(k*nj*ni)+((j-1)*ni)+i]/(fracy_[j-1]*vis_[(k*nj*ni)+((j-1)*ni)+i]+(1.0 - fracy_[j-1])*vis_[(k*nj*ni)+(j*ni)+i]);
                        visve = vis_[(k*nj*ni)+(j*ni)+i+1]*vis_[(k*nj*ni)+((j-1)*ni)+i+1]/(fracy_[j-1]*vis_[(k*nj*ni)+((j-1)*ni)+i+1]+(1.0 - fracy_[j-1])*vis_[(k*nj*ni)+(j*ni)+i+1]);
                        visvw = vis_[(k*nj*ni)+(j*ni)+i-1]*vis_[(k*nj*ni)+((j-1)*ni)+i-1]/(fracy_[j-1]*vis_[(k*nj*ni)+((j-1)*ni)+i-1]+(1.0 - fracy_[j-1])*vis_[(k*nj*ni)+(j*ni)+i-1]);
                        if(i == nim2)
                        {
                            vise = visv;
                        }
                        else
                        {
                            vise = visv*visve/((1.0 - fracx_[i])*visve + fracx_[i]*visv);
                        }//end if

                        if(i == 1)
                        {
                            visw = visv;
                        }
                        else
                        {
                            visw = visv*visvw/(fracx_[i-1]*visvw + (1.0 - fracx_[i-1])*visv);
                        }//end if
                        
                        visvt = vis_[((k+1)*nj*ni)+(j*ni)+i]*vis_[((k+1)*nj*ni)+((j-1)*ni)+i]/(fracy_[j-1]*vis_[((k+1)*nj*ni)+((j-1)*ni)+i]+(1.0 - fracy_[j-1])*vis_[((k+1)*nj*ni)+(j*ni)+i]);
                        visvb = vis_[((k-1)*nj*ni)+(j*ni)+i]*vis_[((k-1)*nj*ni)+((j-1)*ni)+i]/(fracy_[j-1]*vis_[((k-1)*nj*ni)+((j-1)*ni)+i]+(1.0 - fracy_[j-1])*vis_[((k-1)*nj*ni)+(j*ni)+i]);

                        if(k == nkm2)
                        {
                            vist = visv;
                        }
                        else
                        {
                            vist = visv*visvt/((1.0 - fracz_[k])*visvt + fracz_[k]*visv);
                        }//end if

                        if(k == 1)
                        {
                            visb = visv;
                        }
                        else
                        {
                            visb = visv*visvb/(fracz_[k-1]*visvb + (1.0 - fracz_[k-1])*visv);
                        }//end if

                        //----diffusion coefficients----------------------------
                        dn = visn*areaik_[(k*ni)+i]/(yv_[j+1] - yv_[j]);
                        ds = viss*areaik_[(k*ni)+i]/(yv_[j] - yv_[j-1]);
                        de = vise*areavjk_[(k*nj)+j]*dxpwinv_[i+1];
                        dw = visw*areavjk_[(k*nj)+j]*dxpwinv_[i];
                        dt = vist*areavij_[(j*ni)+i]*dzpbinv_[k+1];
                        db = visb*areavij_[(j*ni)+i]*dzpbinv_[k];

                        //----matrix coefficients-------------------------------
                        an_[(k*nj*ni)+(j*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                        as_[(k*nj*ni)+(j*ni)+i] = ds*discFunct(abs(fs)/ds) + max(0.0,fs);
                        ae_[(k*nj*ni)+(j*ni)+i] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                        aw_[(k*nj*ni)+(j*ni)+i] = dw*discFunct(abs(fw)/dw) + max(0.0,fw);
                        at_[(k*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                        ab_[(k*nj*ni)+(j*ni)+i] = db*discFunct(abs(fb)/db) + max(0.0,fb);
                        apnot_[(k*nj*ni)+(j*ni)+i] = rho_[(k*nj*ni)+(j*ni)+i]*volume_v_[(k*nj*ni)+(j*ni)+i]/delt;

                        //----su and sp-----------------------------------------
                        delf = fn - fs + fe - fw + ft - fb;
                        cp0 = max(0.0,delf);
                        cp1 = min(0.0,delf);
                        su_[(k*nj*ni)+(j*ni)+i] = -cp1*vvel_[(k*nj*ni)+(j*ni)+i];
                        su_[(k*nj*ni)+(j*ni)+i] += areaik_[(k*ni)+i]*(pressure_[(k*nj*ni)+((j-1)*ni)+i] - pressure_[(k*nj*ni)+(j*ni)+i]);
                        sp_[(k*nj*ni)+(j*ni)+i] = -cp0;
                        su_[(k*nj*ni)+(j*ni)+i] += apnot_[(k*nj*ni)+(j*ni)+i]*vnot_[(k*nj*ni)+(j*ni)+i];
                        
                        dudyp = (uvel_[(k*nj*ni)+(j*ni)+i+1] - uvel_[(k*nj*ni)+((j-1)*ni)+i+1])*dypsinv_[j];
                        dudym = (uvel_[(k*nj*ni)+(j*ni)+i] - uvel_[(k*nj*ni)+((j-1)*ni)+i])*dypsinv_[j];
                        su_[(k*nj*ni)+(j*ni)+i] += (vise*dudyp - visw*dudym)*areavjk_[(k*nj)+j];
                        dvdyp = (vvel_[(k*nj*ni)+((j+1)*ni)+i] - vvel_[(k*nj*ni)+(j*ni)+i])/(yv_[j+1] - yv_[j]);
                        dvdym = (vvel_[(k*nj*ni)+(j*ni)+i] - vvel_[(k*nj*ni)+((j-1)*ni)+i])/(yv_[j] - yv_[j-1]);
                        su_[(k*nj*ni)+(j*ni)+i] += (visn*dvdyp - viss*dvdym)*areaik_[(k*ni)+i];
                        dwdyp = (wvel_[((k+1)*nj*ni)+(j*ni)+i] - wvel_[((k+1)*nj*ni)+((j-1)*ni)+i])*dypsinv_[j];
                        dwdym = (wvel_[(k*nj*ni)+(j*ni)+i] - wvel_[(k*nj*ni)+((j-1)*ni)+i])*dypsinv_[j];
                        su_[(k*nj*ni)+(j*ni)+i] += (vist*dwdyp - visb*dwdym)*areavij_[(j*ni)+i];
                    }//end for(i)
                }//end for(j)
            }//end for(k)
            break;
        }//end case(2)
        case 3:
        {
            #pragma omp parallel for \
                    private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, vis_w, viswn, visws, visn, viss, viswe, visww, vise, visw,vist, visb, dn, ds, de, dw, dt, db, delf, cp0, cp1, dudzp, dudzm, dvdzp, dvdzm, dwdzp, dwdzm) \
                    collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        vn = vvel_[(k*nj*ni)+((j+1)*ni)+i]*(1.0 - fracz_[k-1]) + vvel_[((k-1)*nj*ni)+((j+1)*ni)+i]*fracz_[k-1];
                        vs = vvel_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracz_[k-1]) + vvel_[((k-1)*nj*ni)+(j*ni)+i]*fracz_[k-1];
                        ue = uvel_[(k*nj*ni)+(j*ni)+i+1]*(1.0 - fracz_[k-1]) + uvel_[((k-1)*nj*ni)+(j*ni)+i+1]*fracz_[k-1];
                        uw = uvel_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracz_[k-1]) + uvel_[((k-1)*nj*ni)+(j*ni)+i]*fracz_[k-1];
                        wt = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k+1)*nj*ni)+(j*ni)+i])*0.5;
                        wb = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k-1)*nj*ni)+(j*ni)+i])*0.5;
                        
                        //----calculate convection coefficients-----------------
                        fn = vn*rho_[(k*nj*ni)+(j*ni)+i]*areawik_[(k*ni)+i];
                        fs = vs*rho_[(k*nj*ni)+(j*ni)+i]*areawik_[(k*ni)+i];
                        fe = ue*rho_[(k*nj*ni)+(j*ni)+i]*areawjk_[(k*nj)+j];
                        fw = uw*rho_[(k*nj*ni)+(j*ni)+i]*areawjk_[(k*nj)+j];
                        ft = wt*rho_[(k*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i];
                        fb = wb*rho_[(k*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i];

                        //----viscosity at cv faces-----------------------------
                        vis_w = vis_[(k*nj*ni)+(j*ni)+i]*vis_[((k-1)*nj*ni)+(j*ni)+i]/(fracz_[k-1]*vis_[((k-1)*nj*ni)+(j*ni)+i] + 
                                    (1.0 - fracz_[k-1])*vis_[(k*nj*ni)+(j*ni)+i]);
                        viswn = vis_[(k*nj*ni)+((j+1)*ni)+i]*vis_[((k-1)*nj*ni)+((j+1)*ni)+i]/(fracz_[k-1]*vis_[((k-1)*nj*ni)+((j+1)*ni)+i] + 
                                    (1.0 - fracz_[k-1])*vis_[(k*nj*ni)+((j+1)*ni)+i]);
                        visws = vis_[(k*nj*ni)+((j-1)*ni)+i]*vis_[((k-1)*nj*ni)+((j-1)*ni)+i]/(fracz_[k-1]*vis_[((k-1)*nj*ni)+((j-1)*ni)+i] + 
                                    (1.0 - fracz_[k-1])*vis_[(k*nj*ni)+((j-1)*ni)+i]);
                        if(j == njm2)
                        {
                            visn = vis_w;
                        }
                        else
                        {
                            visn = vis_w*viswn/((1.0-fracy_[j])*viswn+fracy_[j]*vis_w);
                        }//end if

                        if(j == 1)
                        {
                            viss = vis_w;
                        }
                        else
                        {
                            viss = vis_w*visws/(fracy_[j-1]*visws+(1.0-fracy_[j-1])*vis_w);
                        }//end if

                        viswe = vis_[(k*nj*ni)+(j*ni)+i+1]*vis_[((k-1)*nj*ni)+(j*ni)+i+1]/(fracz_[k-1]*vis_[((k-1)*nj*ni)+(j*ni)+i+1] + 
                                    (1.0 - fracz_[k-1])*vis_[(k*nj*ni)+(j*ni)+i+1]);
                        visww = vis_[(k*nj*ni)+(j*ni)+i-1]*vis_[((k-1)*nj*ni)+(j*ni)+i-1]/(fracz_[k-1]*vis_[((k-1)*nj*ni)+(j*ni)+i-1] + 
                                    (1.0 - fracz_[k-1])*vis_[(k*nj*ni)+(j*ni)+i-1]);
                        if(i == nim2)
                        {
                            vise = vis_w;
                        }
                        else
                        {
                            vise = vis_w*viswe/((1.0 - fracx_[i])*viswe + fracx_[i]*vis_w);
                        }//end if

                        if(i == 1)
                        {
                            visw = vis_w;
                        }
                        else
                        {
                            visw = vis_w*visww/(fracx_[i-1]*visww + (1.0 - fracx_[i-1])*vis_w);
                        }//end if

                        vist = vis_[(k*nj*ni)+(j*ni)+i];
                        visb = vis_[((k-1)*nj*ni)+(j*ni)+i];

                        //----diffusion coefficients----------------------------
                        dn = visn*areawik_[(k*ni)+i]*dypsinv_[j+1];
                        ds = viss*areawik_[(k*ni)+i]*dypsinv_[j];
                        de = vise*areawjk_[(k*nj)+j]*dxpwinv_[i+1];
                        dw = visw*areawjk_[(k*nj)+j]*dxpwinv_[i];
                        dt = vist*areaij_[(j*ni)+i]/(zw_[k+1] - zw_[k]);
                        db = visb*areaij_[(j*ni)+i]/(zw_[k] - zw_[k-1]);

                        //----matrix coefficients-------------------------------
                        an_[(k*nj*ni)+(j*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                        as_[(k*nj*ni)+(j*ni)+i] = ds*discFunct(abs(fs)/ds) + max(0.0,fs);
                        ae_[(k*nj*ni)+(j*ni)+i] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                        aw_[(k*nj*ni)+(j*ni)+i] = dw*discFunct(abs(fw)/dw) + max(0.0,fw);
                        at_[(k*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                        ab_[(k*nj*ni)+(j*ni)+i] = db*discFunct(abs(fb)/db) + max(0.0,fb);
                        apnot_[(k*nj*ni)+(j*ni)+i] = rho_[(k*nj*ni)+(j*ni)+i]*volume_w_[(k*nj*ni)+(j*ni)+i]/delt;

                        //----su and sp-----------------------------------------
                        delf = fn - fs + fe - fw + ft - fb;
                        cp0 = max(0.0,delf);
                        cp1 = min(0.0,delf);
                        su_[(k*nj*ni)+(j*ni)+i] = -cp1*wvel_[(k*nj*ni)+(j*ni)+i];
                        su_[(k*nj*ni)+(j*ni)+i] += areaij_[(j*ni)+i]*(pressure_[((k-1)*nj*ni)+(j*ni)+i] - pressure_[(k*nj*ni)+(j*ni)+i]);
                        sp_[(k*nj*ni)+(j*ni)+i] = -cp0;
                        su_[(k*nj*ni)+(j*ni)+i] += apnot_[(k*nj*ni)+(j*ni)+i]*wnot_[(k*nj*ni)+(j*ni)+i];
                        
                        dudzp = (uvel_[(k*nj*ni)+(j*ni)+i+1] - uvel_[((k-1)*nj*ni)+(j*ni)+i+1])*dzpbinv_[k];
                        dudzm = (uvel_[(k*nj*ni)+(j*ni)+i] - uvel_[((k-1)*nj*ni)+(j*ni)+i])*dzpbinv_[k];
                        su_[(k*nj*ni)+(j*ni)+i] += (vise*dudzp-visw*dudzm)*areawjk_[(k*nj)+j];
                        dvdzp = (vvel_[(k*nj*ni)+((j+1)*ni)+i] - vvel_[((k-1)*nj*ni)+((j+1)*ni)+i])*dzpbinv_[k];
                        dvdzm = (vvel_[(k*nj*ni)+(j*ni)+i] - vvel_[((k-1)*nj*ni)+(j*ni)+i])*dzpbinv_[k];
                        su_[(k*nj*ni)+(j*ni)+i] +=  (visn*dvdzp - viss*dvdzm)*areawik_[(k*ni)+i];
                        dwdzp = (wvel_[((k+1)*nj*ni)+(j*ni)+i] - wvel_[(k*nj*ni)+(j*ni)+i])/(zw_[k+1]-zw_[k]);
                        dwdzm = (wvel_[(k*nj*ni)+(j*ni)+i] - wvel_[((k-1)*nj*ni)+(j*ni)+i])/(zw_[k]-zw_[k-1]);
                        su_[(k*nj*ni)+(j*ni)+i] += (vist*dwdzp - visb*dwdzm)*areaij_[(j*ni)+i];
                    }//end for(i)
                }//end for(j)
            }//end for(k)
            break;
        }//end case(3)
        case 4:
        {
            resorm_ = 0.0;
            #pragma omp parallel for private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, delf) \
                                     reduction(+: resorm_) \
                                     collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        //----main coefficients---------------------------------
                        an_[(k*nj*ni)+(j*ni)+i] = areaik_[(k*ni)+i]*dvy_[(k*nj*ni)+((j+1)*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        as_[(k*nj*ni)+(j*ni)+i] = areaik_[(k*ni)+i]*dvy_[(k*nj*ni)+(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        ae_[(k*nj*ni)+(j*ni)+i] = areajk_[(k*nj)+j]*dux_[(k*nj*ni)+(j*ni)+i+1]*rho_[(k*nj*ni)+(j*ni)+i];
                        aw_[(k*nj*ni)+(j*ni)+i] = areajk_[(k*nj)+j]*dux_[(k*nj*ni)+(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        at_[(k*nj*ni)+(j*ni)+i] = areaij_[(j*ni)+i]*dwz_[((k+1)*nj*ni)+(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        ab_[(k*nj*ni)+(j*ni)+i] = areaij_[(j*ni)+i]*dwz_[(k*nj*ni)+(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];

                        //----velocities at cv faces----------------------------
                        vn = vvel_[(k*nj*ni)+((j+1)*ni)+i];
                        vs = vvel_[(k*nj*ni)+(j*ni)+i];
                        ue = uvel_[(k*nj*ni)+(j*ni)+i+1];
                        uw = uvel_[(k*nj*ni)+(j*ni)+i];
                        wt = wvel_[((k+1)*nj*ni)+(j*ni)+i];
                        wb = wvel_[(k*nj*ni)+(j*ni)+i];

                        fn = vn*areaik_[(k*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        fs = vs*areaik_[(k*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        fe = ue*areajk_[(k*nj)+j]*rho_[(k*nj*ni)+(j*ni)+i];
                        fw = uw*areajk_[(k*nj)+j]*rho_[(k*nj*ni)+(j*ni)+i];
                        ft = wt*areaij_[(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];
                        fb = wb*areaij_[(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i];

                        delf = fn - fs + fe - fw + ft - fb;
                        sp_[(k*nj*ni)+(j*ni)+i] = 0.0;
                        su_[(k*nj*ni)+(j*ni)+i] = -delf;

                        resorm_ += abs(delf);
                    }//end for(i)
                }//end for(j)
            }//end for(k)
            break;
        }//end case(4)
        case 5:
        {
            #pragma omp parallel private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, difn, dife, dift, dn, ds, de, dw, dt, db) 
            {
                #pragma omp for collapse(3) nowait
                for (k=1; k<nkm1; k++) 
                {
                    for (j=1; j<njm1; j++)
                    {
                        for (i=1; i<nim1; i++)
                        {
                            vn = vvel_[(k*nj*ni)+((j+1)*ni)+i];
                            ue = uvel_[(k*nj*ni)+(j*ni)+i+1];
                            wt = wvel_[((k+1)*nj*ni)+(j*ni)+i];
                            
                            //----convection coefficients-----------------------
                            double densbydt = rho_[(k*nj*ni)+(j*ni)+i]/delt;
                            fn = rho_[(k*nj*ni)+((j+1)*ni)+i]*vn*areaik_[(k*ni)+i];
                            fe = rho_[(k*nj*ni)+(j*ni)+i+1]*ue*areajk_[(k*nj)+j];
                            ft = rho_[((k+1)*nj*ni)+(j*ni)+i]*wt*areaij_[(j*ni)+i];

                            //----diffusion coefficients (calculate with harmonic mean)
                            difn = diff_[(k*nj*ni)+(j*ni)+i]*diff_[(k*nj*ni)+((j+1)*ni)+i]/
                                    ((1.0 - fracy_[j])*diff_[(k*nj*ni)+((j+1)*ni)+i] + fracy_[j]*diff_[(k*nj*ni)+(j*ni)+i]);
                            dife = diff_[(k*nj*ni)+(j*ni)+i]*diff_[(k*nj*ni)+(j*ni)+i+1]/
                                    ((1.0 - fracx_[i])*diff_[(k*nj*ni)+(j*ni)+i+1] + fracx_[i]*diff_[(k*nj*ni)+(j*ni)+i]);
                            dift = diff_[(k*nj*ni)+(j*ni)+i]*diff_[((k+1)*nj*ni)+(j*ni)+i]/
                                    ((1.0 - fracz_[k])*diff_[((k+1)*nj*ni)+(j*ni)+i] + fracz_[k]*diff_[(k*nj*ni)+(j*ni)+i]);

                            dn = difn*areaik_[(k*ni)+i]*dypsinv_[j+1];
                            de = dife*areajk_[(k*nj)+j]*dxpwinv_[i+1];
                            dt = dift*areaij_[(j*ni)+i]*dzpbinv_[k+1];

                            //----matrix coefficients---------------------------
                            an_[(k*nj*ni)+(j*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                            as_[(k*nj*ni)+((j+1)*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,fn);
                            ae_[(k*nj*ni)+(j*ni)+i] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                            aw_[(k*nj*ni)+(j*ni)+i+1] = de*discFunct(abs(fe)/de) + max(0.0,fe);
                            at_[(k*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                            ab_[((k+1)*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,ft);
                            apnot_[(k*nj*ni)+(j*ni)+i] = densbydt*volume_[(k*nj*ni)+(j*ni)+i];

                            // sp here is different from BOOK <numerical heat transfer>, equal to sp*deltaV in the book
                            sp_[(k*nj*ni)+(j*ni)+i] = 0.0;

                            // su here is different from BOOK <numerical heat transfer>, equal to b in the book
                            su_[(k*nj*ni)+(j*ni)+i] = apnot_[(k*nj*ni)+(j*ni)+i]*hnot_[(k*nj*ni)+(j*ni)+i];
                        }//end for(i)
                    }//end for(j)
                }//end for(k)

                //----j=1 & j=njm2
                #pragma omp for collapse(2) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        vs = vvel_[(k*nj*ni)+(1*ni)+i];
                        fs = rho_[(k*nj*ni)+(1*ni)+i]*vs*areaik_[(k*ni)+i];
                        ds = diff_[(k*nj*ni)+(1*ni)+i]*areaik_[(k*ni)+i]*dypsinv_[1];
                        as_[(k*nj*ni)+(1*ni)+i] = ds*discFunct(abs(fs)/ds) + max(0.0,fs);

                        vn = vvel_[(k*nj*ni)+(njm1*ni)+i];
                        fn = rho_[(k*nj*ni)+(njm1*ni)+i]*vn*areaik_[(k*ni)+i];
                        dn = diff_[(k*nj*ni)+(njm2*ni)+i]*areaik_[(k*ni)+i]*dypsinv_[njm1];
                        an_[(k*nj*ni)+(njm2*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                    }//end for(i)
                }//end for(k)

                //----i=1 & i=nim2
                #pragma omp for collapse(2) nowait 
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {
                        uw = uvel_[(k*nj*ni)+(j*ni)+1];
                        fw = rho_[(k*nj*ni)+(j*ni)+1]*uw*areajk_[(k*nj)+j];
                        dw = diff_[(k*nj*ni)+(j*ni)+1]*areajk_[(k*nj)+j]*dxpwinv_[1];
                        aw_[(k*nj*ni)+(j*ni)+1] = dw*discFunct(abs(fw)/dw) + max(0.0,fw);

                        ue = uvel_[(k*nj*ni)+(j*ni)+nim1];
                        fe = rho_[(k*nj*ni)+(j*ni)+nim1]*ue*areajk_[(k*nj)+j];
                        de = diff_[(k*nj*ni)+(j*ni)+nim2]*areajk_[(k*nj)+j]*dxpwinv_[nim1];
                        ae_[(k*nj*ni)+(j*ni)+nim2] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                    }//end for(j)
                }//end for(k)

                //----k=1 & k=nkm2
                #pragma omp for collapse(2) nowait 
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        wb = wvel_[(1*nj*ni)+(j*ni)+i];
                        fb = rho_[(1*nj*ni)+(j*ni)+i]*wb*areaij_[(j*ni)+i];
                        db = diff_[(1*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i]*dzpbinv_[1];
                        ab_[(1*nj*ni)+(j*ni)+i] = db*discFunct(abs(fb)/db) + max(0.0,fb);

                        wt = wvel_[(nkm1*nj*ni)+(j*ni)+i];
                        ft = rho_[(nkm1*nj*ni)+(j*ni)+i]*wt*areaij_[(j*ni)+i];
                        dt = diff_[(nkm2*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i]*dzpbinv_[nkm1];
                        at_[(nkm2*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                    }//end for(i)
                }//end for(j)
            }//end omp parallel

            break;
        }//end case(5)
        case 6:
        {
            // grab some stuff from the domain
            double tsolid = domainMgr_->tsolid_;
            double diffspeciesl = domainMgr_->massdiff_;

            // define some local variables
            vector<double> diffspecies(nk*nj*ni);
            
            for (int loc = 0; loc<ni*nj*nk; loc++)
            {
                if(temp_[loc] <= tsolid)
                {
            	    diffspecies[loc] = 1e-11;
                }
                else
                {
                    diffspecies[loc] = diffspeciesl*rho_[loc];
                }//end if
            }//end for(loc)

            #pragma omp parallel private(vn, vs, ue, uw, wt, wb, fn, fs, fe, fw, ft, fb, dn, ds, de, dw, dt, db) 
            {
                #pragma omp for collapse(3) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {  
                        for (i=1; i<nim1; i++)
                        {
                            vn = vvel_[(k*nj*ni)+((j+1)*ni)+i]*2.0;
                            ue = uvel_[(k*nj*ni)+(j*ni)+i+1]*2.0;
                            wt = wvel_[((k+1)*nj*ni)+(j*ni)+i]*2.0;
                            
                	    //----convection coefficients-----------------------
                            fn = rho_[(k*nj*ni)+((j+1)*ni)+i]*vn*areaik_[(k*ni)+i];
                            fe = rho_[(k*nj*ni)+(j*ni)+i+1]*ue*areajk_[(k*nj)+j];
                            ft = rho_[((k+1)*nj*ni)+(j*ni)+i]*wt*areaij_[(j*ni)+i];

                	    //----diffusion coefficients------------------------
                            dn = diffspecies[(k*nj*ni)+(j*ni)+i]*areaik_[(k*ni)+i]*dypsinv_[j+1];
                            de = diffspecies[(k*nj*ni)+(j*ni)+i]*areajk_[(k*nj)+j]*dxpwinv_[i+1];
                            dt = diffspecies[(k*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i]*dzpbinv_[k+1];

                            //----matrix coefficients---------------------------
                            double densbydt = rho_[(k*nj*ni)+(j*ni)+i]/delt;
                            an_[(k*nj*ni)+(j*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,-fn);
                            as_[(k*nj*ni)+((j+1)*ni)+i] = dn*discFunct(abs(fn)/dn) + max(0.0,fn);
                            ae_[(k*nj*ni)+(j*ni)+i] = de*discFunct(abs(fe)/de) + max(0.0,-fe);
                            aw_[(k*nj*ni)+(j*ni)+i+1] = de*discFunct(abs(fe)/de) + max(0.0,fe);
                            at_[(k*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,-ft);
                            ab_[((k+1)*nj*ni)+(j*ni)+i] = dt*discFunct(abs(ft)/dt) + max(0.0,ft);
                            apnot_[(k*nj*ni)+(j*ni)+i] = densbydt*volume_[(k*nj*ni)+(j*ni)+i];

                            // sp here is different from BOOK <numerical heat transfer>, equal to sp*deltaV in the book
                            sp_[(k*nj*ni)+(j*ni)+i] = 0.0;

                            // su here is different from BOOK <numerical heat transfer>, equal to b in the book
                            su_[(k*nj*ni)+(j*ni)+i] = apnot_[(k*nj*ni)+(j*ni)+i]*avgconcentrationnot_[(k*nj*ni)+(j*ni)+i];
                        }//end for(i)
                    }//end for(j)
                }//end for(k)	

                //----j=1
                #pragma omp for collapse(2) nowait 
                for (k=1; k<nkm1; k++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        vs = vvel_[(k*nj*ni)+(1*ni)+i]*2;
                        fs = rho_[(k*nj*ni)+(1*ni)+i]*vs*areaik_[(k*ni)+i];
                        ds = diffspecies[(k*nj*ni)+(1*ni)+i]*areaik_[(k*ni)+i]*dypsinv_[1];
                        as_[(k*nj*ni)+(1*ni)+i] = ds*discFunct(abs(fs)/ds) + max(0.0,fs);
                    }//end for(i)
                }//end for(k)

                //----i=1
                #pragma omp for collapse(2) nowait 
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {
                        uw = uvel_[(k*nj*ni)+(j*ni)+1]*2;
                        fw = rho_[(k*nj*ni)+(j*ni)+1]*uw*areajk_[(k*nj)+j];
                        dw = diffspecies[(k*nj*ni)+(j*ni)+1]*areajk_[(k*nj)+j]*dxpwinv_[1];
                        aw_[(k*nj*ni)+(j*ni)+1] = dw*discFunct(abs(fw)/dw) + max(0.0,fw);
                    }//end for(j)
                }//end for(k)

                //----k=1
                #pragma omp for collapse(2) nowait 
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        wb = wvel_[(1*nj*ni)+(j*ni)+i]*2;
                        fb = rho_[(1*nj*ni)+(j*ni)+i]*wb*areaij_[(j*ni)+i];
                        db = diffspecies[(1*nj*ni)+(j*ni)+i]*areaij_[(j*ni)+i]*dzpbinv_[1];
                        ab_[(1*nj*ni)+(j*ni)+i] = db*discFunct(abs(fb)/db) + max(0.0,fb);
                    }//end for(i)
                }//end for(j)
            }//end omp parallel

            break;
        }//end case(6)
    }//end switch
}//end discretizeEquations

//////////////////////////////////////////////////////
//		 calculateSource                    //
//////////////////////////////////////////////////////
void 
CFDSolverManager::calculateSource(int &ivar)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;
    int nkm2 = domainMgr_->nkm2_;
    bool steady = domainMgr_->steady_;
    double tsolid = domainMgr_->tsolid_;
    double scanvel = domainMgr_->scanvel_;
    double a_oxygen = domainMgr_->a_oxygen_;
    double darcyKo = domainMgr_->darcyKo_;
    
    // grab some stuff from the bc
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;
     
    // define some local variables
    int i, j, k;
    double fraclu, fraclv, fraclw, tw;
    double term, term1;
    double tulc, tvlc, twlc, rhoscan;


    // calculate the source
    switch(ivar)
    {
        case 1:
        {
            double urfu = domainMgr_->urfu_;
            #pragma omp parallel
            {
                #pragma omp for private(fraclu, term, term1, rhoscan) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            fraclu = fracl_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracx_[i-1]) + fracl_[(k*nj*ni)+(j*ni)+i-1]*fracx_[i-1];
                            rhoscan = rho_[(k*nj*ni)+(j*ni)+i]*scanvel;
                            if(fraclu > 0)  
                            {
                                //----mushy zone
                                term = darcyKo*((1.0 - fraclu)*(1.0 - fraclu))/((fraclu*fraclu*fraclu) + small_);
                                sp_[(k*nj*ni)+(j*ni)+i] -= term*volume_u_[(k*nj*ni)+(j*ni)+i];
                               
                                //----scanning velocity
                                if(steady)  
                                {
                                    term1 = rhoscan*areajk_[(k*nj)+j]*(uvel_[(k*nj*ni)+(j*ni)+i-1] - uvel_[(k*nj*ni)+(j*ni)+i]);
                                    su_[(k*nj*ni)+(j*ni)+i] += term1;
                                }//end if
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
        
                //----k=nkm1
                #pragma omp for collapse(2) nowait
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        su_[(nkm2*nj*ni)+(j*ni)+i] += at_[(nkm2*nj*ni)+(j*ni)+i]*uvel_[(nkm1*nj*ni)+(j*ni)+i];
                        sp_[(nkm2*nj*ni)+(j*ni)+i] -= at_[(nkm2*nj*ni)+(j*ni)+i];
                        at_[(nkm2*nj*ni)+(j*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(j)
        
                #pragma omp for private(tulc) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                                                      ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                                                      at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] + 
                                                      apnot_[(k*nj*ni)+(j*ni)+i] - sp_[(k*nj*ni)+(j*ni)+i];

                            dux_[(k*nj*ni)+(j*ni)+i] = areajk_[(k*nj)+j]/ap_[(k*nj*ni)+(j*ni)+i]; 
        
                            //----under-relaxation
                            ap_[(k*nj*ni)+(j*ni)+i] /= urfu;
                            su_[(k*nj*ni)+(j*ni)+i] += (1.0 - urfu)*ap_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+(j*ni)+i];
                            dux_[(k*nj*ni)+(j*ni)+i] *= urfu;
        
                            //----zero velocity
                            tulc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+(j*ni)+i-1]);
                            if(tulc <= tsolid)  
                            {
                                su_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                an_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                as_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ae_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                aw_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                at_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ab_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ap_[(k*nj*ni)+(j*ni)+i] = large_;
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
            }//end omp parallel

            break;
        }//end case(1)
        case 2:
        {
            double urfv = domainMgr_->urfv_;
            #pragma omp parallel
            {
                #pragma omp for private(fraclv, term, term1, rhoscan) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            fraclv = fracl_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracy_[j-1]) + fracl_[(k*nj*ni)+((j-1)*ni)+i] * fracy_[j-1];
                            rhoscan = rho_[(k*nj*ni)+(j*ni)+i]*scanvel;
                            if(fraclv > 0)  
                            {
                                //----mushy zone
                                term = darcyKo*((1.0 - fraclv)*(1.0 - fraclv))/((fraclv*fraclv*fraclv) + small_);
                                sp_[(k*nj*ni)+(j*ni)+i] -= term*volume_v_[(k*nj*ni)+(j*ni)+i];

                                //----scanning velocity
                                if(steady)  
                                {
                                    term1 = rhoscan*areavjk_[(k*nj)+j]*(vvel_[(k*nj*ni)+(j*ni)+i-1] - vvel_[(k*nj*ni)+(j*ni)+i]);
                                    su_[(k*nj*ni)+(j*ni)+i] += term1;
                                }//end if
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)

                //----k=nk
                #pragma omp for collapse(2) nowait
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        su_[(nkm2*nj*ni)+(j*ni)+i] += at_[(nkm2*nj*ni)+(j*ni)+i]*vvel_[(nkm1*nj*ni)+(j*ni)+i];
                        sp_[(nkm2*nj*ni)+(j*ni)+i] -= at_[(nkm2*nj*ni)+(j*ni)+i];
                        at_[(nkm2*nj*ni)+(j*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(j)

                #pragma omp for private(tvlc) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                                                      ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                                                      at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] + 
                                                      apnot_[(k*nj*ni)+(j*ni)+i] - sp_[(k*nj*ni)+(j*ni)+i];
                            dvy_[(k*nj*ni)+(j*ni)+i] = areaik_[(k*ni)+i]/ap_[(k*nj*ni)+(j*ni)+i];

                            //----under-relaxation
                            ap_[(k*nj*ni)+(j*ni)+i] /= urfv;
                            su_[(k*nj*ni)+(j*ni)+i] += (1.0 - urfv)*ap_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+(j*ni)+i];
                            dvy_[(k*nj*ni)+(j*ni)+i] *= urfv;

                            //----zero velocity
                            tvlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+((j-1)*ni)+i]);
                            if(tvlc <= tsolid)  
                            {
                                su_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                an_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                as_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ae_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                aw_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                at_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ab_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ap_[(k*nj*ni)+(j*ni)+i] = large_;
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
            }//end omp parallel

            break;
        }//end case(2)
        case 3:
        {
            double urfw = domainMgr_->urfw_;
            #pragma omp parallel
            {
                #pragma omp for private(fraclw, tw, term, term1, rhoscan) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            fraclw = fracl_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracz_[k-1])+fracl_[((k-1)*nj*ni)+(j*ni)+i]*fracz_[k-1];
                            rhoscan = rho_[(k*nj*ni)+(j*ni)+i]*scanvel;
                            double boufac = rho_[(k*nj*ni)+(j*ni)+i]*meshObj_->grav_*domainMgr_->beta_; 
                            if(fraclw > 0)  
                            {
                                //----mushy zone
                                term = darcyKo*((1.0 - fraclw)*(1.0 - fraclw))/((fraclw*fraclw*fraclw) + small_);
                                sp_[(k*nj*ni)+(j*ni)+i] -= term*volume_w_[(k*nj*ni)+(j*ni)+i];

                                //----scanning velocity
                                if(steady)  
                                {
                                    term1 = rhoscan*areawjk_[(k*nj)+j]*(wvel_[(k*nj*ni)+(j*ni)+i-1] - wvel_[(k*nj*ni)+(j*ni)+i]);
                                    su_[(k*nj*ni)+(j*ni)+i] += term1;
                                }//end if

                                //----buoyancy
                                tw = temp_[(k*nj*ni)+(j*ni)+i]*(1.0 - fracz_[k-1]) + temp_[((k-1)*nj*ni)+(j*ni)+i]*fracz_[k-1];
                                su_[(k*nj*ni)+(j*ni)+i] += boufac*volume_w_[(k*nj*ni)+(j*ni)+i]*(tw-tsolid);
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)

                #pragma omp for private(twlc) collapse(3) nowait
                for (k=kstat; k<nkm1; k++)
                {
                    for (j=jstat; j<=jend; j++)
                    {
                        for (i=istatp1; i<=iendm1; i++)
                        {
                            ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                                                      ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                                                      at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] + 
                                                      apnot_[(k*nj*ni)+(j*ni)+i] - sp_[(k*nj*ni)+(j*ni)+i];
                            dwz_[(k*nj*ni)+(j*ni)+i] = areaij_[(j*ni)+i]/ap_[(k*nj*ni)+(j*ni)+i];

                            //----under-relaxation
                            ap_[(k*nj*ni)+(j*ni)+i] /= urfw;
                            su_[(k*nj*ni)+(j*ni)+i] += (1.0 - urfw)*ap_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+(j*ni)+i];
                            dwz_[(k*nj*ni)+(j*ni)+i] *= urfw;

                            //----zero velocity
                            twlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[((k-1)*nj*ni)+(j*ni)+i]);
                            if(twlc <= tsolid)  
                            {
                                su_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                an_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                as_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ae_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                aw_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                at_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ab_[(k*nj*ni)+(j*ni)+i] = 0.0;
                                ap_[(k*nj*ni)+(j*ni)+i] = large_;
                            }//end if
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
            }//end omp parallel

            break;
        }//end case(3)
        case 4:
        {
            #pragma omp parallel for collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                                                  ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                                                  at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] - 
                                                  sp_[(k*nj*ni)+(j*ni)+i];

                        //----zero pressure
                        if(temp_[(k*nj*ni)+(j*ni)+i] <= tsolid) 
                        {
                            su_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            an_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            as_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            ae_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            aw_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            at_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            ab_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            ap_[(k*nj*ni)+(j*ni)+i] = large_;
                        }//end if
                    }//end for(i)
                }//end for(j)
            }//end for(k)

            break;
        }//end case(4)
        case 5:
        {
            // grab some stuff from the domain
            int nim1 = domainMgr_->nim1_;
            int njm1 = domainMgr_->njm1_;
            int nim2 = domainMgr_->nim2_;
            int njm2 = domainMgr_->njm2_;
            double urfh = domainMgr_->urfh_;
            double hlatnt = domainMgr_->hlatnt_;
            double avolfact = domainMgr_->avolfact_;
            double avolpow = domainMgr_->avolpow_;
            double apowseta = domainMgr_->apowseta_;
            double heatthick = domainMgr_->heatthick_;
            double heatrb = domainMgr_->heatrb_;
             
            // grab some stuff from the bc
            double beamposx = bcMgr_->beamposx_;
            double beamposy = bcMgr_->beamposy_;

            //----source term
            #pragma omp parallel
            {
                double volht, flew, flns, fltb, variable1;

                #pragma omp for collapse(3) private(volht, flew, flns, fltb, variable1) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {
                        for (i=1; i<nim1; i++)
                        {
                            if(z_[nkm1] - z_[k] <= heatthick)  
                            {
                                sourceinput_[(k*nj*ni)+(j*ni)+i] = avolpow*avolfact/M_PI/(heatrb*heatrb)/heatthick* 
                                                     apowseta*exp(-avolfact/(heatrb*heatrb)*
                                                    (((beamposx-x_[i])*(beamposx-x_[i])) + 
                                                    ((beamposy-y_[j])*(beamposy-y_[j]))))*(double)bcMgr_->powerindicator_;

                                su_[(k*nj*ni)+(j*ni)+i] += volume_[(k*nj*ni)+(j*ni)+i]*sourceinput_[(k*nj*ni)+(j*ni)+i];
                            }
                            else
                            {
                                sourceinput_[(k*nj*ni)+(j*ni)+i] = 0.0;
                            }//end if

                            // source terms due to latent heat (two terms)
                            variable1 = hlatnt*rho_[(k*nj*ni)+(j*ni)+i]/meshObj_->delt_;
                            volht = volume_[(k*nj*ni)+(j*ni)+i]*variable1;
                            su_[(k*nj*ni)+(j*ni)+i] -= volht*(fracl_[(k*nj*ni)+(j*ni)+i] - fraclnot_[(k*nj*ni)+(j*ni)+i]);

                            flew = areajk_[(k*nj)+j]*(max(uvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i-1] - 
                                                      max(-uvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i] +
                                                      max(-uvel_[(k*nj*ni)+(j*ni)+i+1],0.0)*fracl_[(k*nj*ni)+(j*ni)+i+1] - 
                                                      max(uvel_[(k*nj*ni)+(j*ni)+i+1],0.0)*fracl_[(k*nj*ni)+(j*ni)+i]);

                            flns = areaik_[(k*ni)+i]*(max(vvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+((j-1)*ni)+i] - 
                                                      max(-vvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i] +
                                                      max(-vvel_[(k*nj*ni)+((j+1)*ni)+i],0.0)*fracl_[(k*nj*ni)+((j+1)*ni)+i] - 
                                                      max(vvel_[(k*nj*ni)+((j+1)*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i]);

                            fltb = areaij_[(j*ni)+i]*(max(wvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[((k-1)*nj*ni)+(j*ni)+i] - 
                                                      max(-wvel_[(k*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i] +
                                                      max(-wvel_[((k+1)*nj*ni)+(j*ni)+i],0.0)*fracl_[((k+1)*nj*ni)+(j*ni)+i] - 
                                                      max(wvel_[((k+1)*nj*ni)+(j*ni)+i],0.0)*fracl_[(k*nj*ni)+(j*ni)+i]);

                            su_[(k*nj*ni)+(j*ni)+i] += rho_[(k*nj*ni)+(j*ni)+i]*hlatnt*(flew + flns + fltb);

                        }//end for(i)
                    }//end for(j)
                }//end for(k)

                
                //----k=1 & k=nkm2   transfer discrete coefficients at boundary nodes to source term, easy to solve
                #pragma omp for collapse(2) nowait
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        su_[(1*nj*ni)+(j*ni)+i] += ab_[(1*nj*ni)+(j*ni)+i]*enthalpy_[(0*nj*ni)+(j*ni)+i];
                        sp_[(1*nj*ni)+(j*ni)+i] -= ab_[(1*nj*ni)+(j*ni)+i];
                        ab_[(1*nj*ni)+(j*ni)+i] = 0.0;

                        su_[(nkm2*nj*ni)+(j*ni)+i] += at_[(nkm2*nj*ni)+(j*ni)+i]*enthalpy_[(nkm1*nj*ni)+(j*ni)+i];
                        sp_[(nkm2*nj*ni)+(j*ni)+i] -= at_[(nkm2*nj*ni)+(j*ni)+i];
                        at_[(nkm2*nj*ni)+(j*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(j)

                //----j=1 & j=njm2
                #pragma omp for collapse(2) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        su_[(k*nj*ni)+(1*ni)+i] += as_[(k*nj*ni)+(1*ni)+i]*enthalpy_[(k*nj*ni)+(0*ni)+i];
                        sp_[(k*nj*ni)+(1*ni)+i] -= as_[(k*nj*ni)+(1*ni)+i];
                        as_[(k*nj*ni)+(1*ni)+i] = 0.0;
    
                        su_[(k*nj*ni)+(njm2*ni)+i] += an_[(k*nj*ni)+(njm2*ni)+i]*enthalpy_[(k*nj*ni)+(njm1*ni)+i];
                        sp_[(k*nj*ni)+(njm2*ni)+i] -= an_[(k*nj*ni)+(njm2*ni)+i];
                        an_[(k*nj*ni)+(njm2*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(k)
                
                //---i=1 & i=nim2
                #pragma omp for collapse(2) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {
                        su_[(k*nj*ni)+(j*ni)+1] += aw_[(k*nj*ni)+(j*ni)+1]*enthalpy_[(k*nj*ni)+(j*ni)+0];
                        sp_[(k*nj*ni)+(j*ni)+1] -= aw_[(k*nj*ni)+(j*ni)+1];
                        aw_[(k*nj*ni)+(j*ni)+1] = 0.0;
    
                        su_[(k*nj*ni)+(j*ni)+nim2] += ae_[(k*nj*ni)+(j*ni)+nim2]*enthalpy_[(k*nj*ni)+(j*ni)+nim1];
                        sp_[(k*nj*ni)+(j*ni)+nim2] -= ae_[(k*nj*ni)+(j*ni)+nim2];
                        ae_[(k*nj*ni)+(j*ni)+nim2] = 0.0;
                    }//end for(j)
                }//end for(k)
    
                // calculate ap, the defination of sp is different from <numerical heat transfer>, equals to sp*deltaV in that BOOK
                #pragma omp for collapse(3) nowait
                for (k=1; k<nkm1; k++) 
                {
                    for (j=1; j<njm1; j++)
                    {
                        for (i=1; i<nim1; i++)
                        {
                            ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                                                      ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                                                      at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] + 
                                                      apnot_[(k*nj*ni)+(j*ni)+i] - sp_[(k*nj*ni)+(j*ni)+i];

                            //----under-relaxation (restrict the change rate of enthalpy at each iteration)
                            ap_[(k*nj*ni)+(j*ni)+i] /= urfh;
                            su_[(k*nj*ni)+(j*ni)+i] += (1.0 - urfh)*ap_[(k*nj*ni)+(j*ni)+i]*enthalpy_[(k*nj*ni)+(j*ni)+i];
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
            }//end omp parallel

            break;
        }//end case(5)
        case 6:
        {
            // grab some stuff from the domain
            int nim1 = domainMgr_->nim1_;
            int njm1 = domainMgr_->njm1_;
            int nim2 = domainMgr_->nim2_;
            int njm2 = domainMgr_->njm2_;
            double urfc = domainMgr_->urfc_;

            //----k=1 & k=nkm2  transfer discrete coefficients at boundary nodes to source term, easy to solve
            #pragma omp parallel 
            {
                //----k=1 & i=nkm2
                #pragma omp for collapse(2) nowait
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        su_[(1*nj*ni)+(j*ni)+i] += ab_[(1*nj*ni)+(j*ni)+i]*avgconcentration_[(0*nj*ni)+(j*ni)+i];
                        sp_[(1*nj*ni)+(j*ni)+i] -= ab_[(1*nj*ni)+(j*ni)+i];
                        ab_[(1*nj*ni)+(j*ni)+i] = 0.0;

                        su_[(nkm2*nj*ni)+(j*ni)+i] += at_[(nkm2*nj*ni)+(j*ni)+i]*avgconcentration_[(nkm1*nj*ni)+(j*ni)+i];
                        sp_[(nkm2*nj*ni)+(j*ni)+i] -= at_[(nkm2*nj*ni)+(j*ni)+i];
                        at_[(nkm2*nj*ni)+(j*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(j)

                //----j=1 & j=njm2
                #pragma omp for collapse(2) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        su_[(k*nj*ni)+(1*ni)+i] += as_[(k*nj*ni)+(1*ni)+i]*avgconcentration_[(k*nj*ni)+(0*ni)+i];
                        sp_[(k*nj*ni)+(1*ni)+i] -= as_[(k*nj*ni)+(1*ni)+i];
                        as_[(k*nj*ni)+(1*ni)+i] = 0.0;
                
                        su_[(k*nj*ni)+(njm2*ni)+i] += an_[(k*nj*ni)+(njm2*ni)+i]*avgconcentration_[(k*nj*ni)+(njm1*ni)+i];
                        sp_[(k*nj*ni)+(njm2*ni)+i] -= an_[(k*nj*ni)+(njm2*ni)+i];
                        an_[(k*nj*ni)+(njm2*ni)+i] = 0.0;
                    }//end for(i)
                }//end for(k)
                
                //---i=1 & i=nim2
                #pragma omp for collapse(2) nowait
                for (k=1; k<nkm1; k++)
                {
                    for (j=1; j<njm1; j++)
                    {
                        su_[(k*nj*ni)+(j*ni)+1] += aw_[(k*nj*ni)+(j*ni)+1]*avgconcentration_[(k*nj*ni)+(j*ni)+0];
                        sp_[(k*nj*ni)+(j*ni)+1] -= aw_[(k*nj*ni)+(j*ni)+1];
                        aw_[(k*nj*ni)+(j*ni)+1] = 0.0;
                
                        su_[(k*nj*ni)+(j*ni)+nim2] += ae_[(k*nj*ni)+(j*ni)+nim2]*avgconcentration_[(k*nj*ni)+(j*ni)+nim1];
                        sp_[(k*nj*ni)+(j*ni)+nim2] -= ae_[(k*nj*ni)+(j*ni)+nim2];
                        ae_[(k*nj*ni)+(j*ni)+nim2] = 0.0;
                    }//end for(j)
                }//end for(k)
                
                // calculate ap, the defination of sp is different from <numerical heat transfer>, equals to sp*deltaV in that BOOK
                #pragma omp for collapse(3) nowait
                for (k=1; k<nkm1; k++) 
                {
                   for (j=1; j<njm1; j++)
                   {
                       for (i=1; i<nim1; i++)
                       {
                            ap_[(k*nj*ni)+(j*ni)+i] = an_[(k*nj*ni)+(j*ni)+i] + as_[(k*nj*ni)+(j*ni)+i] + 
                           			          ae_[(k*nj*ni)+(j*ni)+i] + aw_[(k*nj*ni)+(j*ni)+i] + 
                           			          at_[(k*nj*ni)+(j*ni)+i] + ab_[(k*nj*ni)+(j*ni)+i] + 
                           			          apnot_[(k*nj*ni)+(j*ni)+i] - sp_[(k*nj*ni)+(j*ni)+i];

                            //----under-relaxation (restrict the change rate of avgconcentration at_ each iteration)
                            ap_[(k*nj*ni)+(j*ni)+i] /= urfc;
                            su_[(k*nj*ni)+(j*ni)+i] += (1.0 - urfc)*ap_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[(k*nj*ni)+(j*ni)+i];
                        }//end for(i)
                    }//end for(j)
                }//end for(k)
            }//end omp parallel

            break;
        }//end case(6)
    }//end switch
}//end calculateSource

//////////////////////////////////////////////////////
//		 calculateResiduals                 //
//////////////////////////////////////////////////////
void 
CFDSolverManager::calculateResiduals(int &ivar)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    double dens = domainMgr_->dens_;
    
    // grab some stuff from the bc
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;
    double width = bcMgr_->width_;

    // define local variables
    int i, j, k;
    double resor; 
    
    // solve for the residual of the respective equation
    switch(ivar)
    {
        case 1:
        {
            double umaxt = 0.0;
            double sumd = 0.0;
            #pragma omp parallel for private(resor) reduction(+: sumd) collapse(3)
            for (k=kstat; k<nkm1; k++)
            {
                for (j=jstat; j<=jend; j++)
                {
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        resor = an_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+((j+1)*ni)+i] + 
                                as_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+((j-1)*ni)+i] + 
                                ae_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+(j*ni)+i+1] + 
                                aw_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+(j*ni)+i-1] +
                                at_[(k*nj*ni)+(j*ni)+i]*uvel_[((k+1)*nj*ni)+(j*ni)+i] + 
                                ab_[(k*nj*ni)+(j*ni)+i]*uvel_[((k-1)*nj*ni)+(j*ni)+i] +
                                su_[(k*nj*ni)+(j*ni)+i] - ap_[(k*nj*ni)+(j*ni)+i]*uvel_[(k*nj*ni)+(j*ni)+i];        
                        sumd += abs(resor);
                    }//end for(i)
                }//end for(j)
            }//end for(k)

            for (j=jstat; j<=jend; j++)
            {
                for (i=istatp1; i<=iendm1; i++)
                    umaxt = max(umaxt,abs(uvel_[(nkm1*nj*ni)+(j*ni)+i]));
            } 
              // reference momentum
            refmom_ = 0.25*M_PI*(width*width)*dens*(umaxt*umaxt);
              // normalized residual
            resoru_ = sumd/refmom_;

            break;
        }//end case(1)
        case 2:
        {
            double sumd = 0.0;
            #pragma omp parallel for private(resor) reduction(+: sumd) collapse(3)
            for (k=kstat; k<nkm1; k++) 
            {
                for (j=jstat; j<=jend; j++) 
                {
                    for (i=istatp1; i<=iendm1; i++) 
                    {
                        resor = an_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+((j+1)*ni)+i] + 
                                as_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+((j-1)*ni)+i] + 
                                ae_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+(j*ni)+i+1] + 
                                aw_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+(j*ni)+i-1] +
                                at_[(k*nj*ni)+(j*ni)+i]*vvel_[((k+1)*nj*ni)+(j*ni)+i] + 
                                ab_[(k*nj*ni)+(j*ni)+i]*vvel_[((k-1)*nj*ni)+(j*ni)+i] + 
                                su_[(k*nj*ni)+(j*ni)+i] - ap_[(k*nj*ni)+(j*ni)+i]*vvel_[(k*nj*ni)+(j*ni)+i];

                        sumd += abs(resor);
                    }//end for(i)
                }//end for(j)
            }//end for(k)

            resorv_ = sumd/refmom_;
            break;
        }//end case(2)
        case 3:
        {
            double sumd = 0.0;
            #pragma omp parallel for private(resor) reduction(+: sumd) collapse(3)
            for (k=kstat; k<nkm1; k++) 
            {
                for (j=jstat; j<=jend; j++) 
                {
                    for (i=istatp1; i<=iendm1; i++) 
                    {
                        resor = an_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+((j+1)*ni)+i] + 
                                as_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+((j-1)*ni)+i] + 
                                ae_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+(j*ni)+i+1] + 
                                aw_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+(j*ni)+i-1] +
                                at_[(k*nj*ni)+(j*ni)+i]*wvel_[((k+1)*nj*ni)+(j*ni)+i] + 
                                ab_[(k*nj*ni)+(j*ni)+i]*wvel_[((k-1)*nj*ni)+(j*ni)+i] + 
                                su_[(k*nj*ni)+(j*ni)+i] - ap_[(k*nj*ni)+(j*ni)+i]*wvel_[(k*nj*ni)+(j*ni)+i];

                        sumd += abs(resor);
                    }//end for(i)
                }//end for(j)
            }//end for(k)

            resorw_ = sumd/refmom_;
            break;
        }//end case(3)
        case 4:
        {
            // normalized mass source
            double dtpvar;
            double denom = 0.0;

            #pragma omp parallel for private(dtpvar) reduction(+: denom) collapse(3)
            for (k=kstat; k<nkm1; k++) 
            {
                for (j=jstat; j<=jend; j++) 
                {
                    for (i=istatp1; i<=iendm1; i++) 
                    {
                        dtpvar = (abs(uvel_[(k*nj*ni)+(j*ni)+i]) + abs(uvel_[(k*nj*ni)+(j*ni)+i+1]))*areajk_[(k*nj)+j] + 
                                 (abs(vvel_[(k*nj*ni)+(j*ni)+i]) + abs(vvel_[(k*nj*ni)+((j+1)*ni)+i]))*areaik_[(k*ni)+i] + 
                                 (abs(wvel_[(k*nj*ni)+(j*ni)+i]) + abs(wvel_[((k+1)*nj*ni)+(j*ni)+i]))*areaij_[(j*ni)+i];

                        denom += 0.5*abs(dtpvar);
                    }//end for(i)
                }//end for(j)
            }//end for(k)

            denom *= dens;
            resorm_ /= (denom + small_);
            break;
        }//end case(4)
        case 5:
        {
            double sumh = 0.0;
            double sumd = 0.0;

            #pragma omp parallel for private(resor) reduction(+: sumd, sumh) collapse(3)
            for (k=1; k<nkm1; k++) 
            {
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        resor = (an_[(k*nj*ni)+(j*ni)+i]*enthalpy_[(k*nj*ni)+((j+1)*ni)+i] + 
                                 as_[(k*nj*ni)+(j*ni)+i]*enthalpy_[(k*nj*ni)+((j-1)*ni)+i] +
                                 ae_[(k*nj*ni)+(j*ni)+i]*enthalpy_[(k*nj*ni)+(j*ni)+i+1] + 
                                 aw_[(k*nj*ni)+(j*ni)+i]*enthalpy_[(k*nj*ni)+(j*ni)+i-1] +
                                 at_[(k*nj*ni)+(j*ni)+i]*enthalpy_[((k+1)*nj*ni)+(j*ni)+i] + 
                                 ab_[(k*nj*ni)+(j*ni)+i]*enthalpy_[((k-1)*nj*ni)+(j*ni)+i] + 
                                 su_[(k*nj*ni)+(j*ni)+i])/ap_[(k*nj*ni)+(j*ni)+i] - enthalpy_[(k*nj*ni)+(j*ni)+i];

                        sumd += abs(resor);
                        sumh += abs(enthalpy_[(k*nj*ni)+(j*ni)+i]);
                    }//end for(i)
                }//end for(j)
            }//end for(k)
 
            resorh_ = sumd/(sumh + small_);
            break;
        }//end case(5)
        case 6:
        {
            double sumc = 0.0;
            double sumd = 0.0;

            #pragma omp parallel for private(resor) reduction(+: sumd, sumc) collapse(3)
            for (k=1; k<nkm1; k++) 
            {
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        resor = (an_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[(k*nj*ni)+((j+1)*ni)+i] + 
                                 as_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[(k*nj*ni)+((j-1)*ni)+i] +
                                 ae_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[(k*nj*ni)+(j*ni)+i+1] + 
                                 aw_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[(k*nj*ni)+(j*ni)+i-1] +
                                 at_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[((k+1)*nj*ni)+(j*ni)+i] + 
                                 ab_[(k*nj*ni)+(j*ni)+i]*avgconcentration_[((k-1)*nj*ni)+(j*ni)+i] + 
                                 su_[(k*nj*ni)+(j*ni)+i])/ap_[(k*nj*ni)+(j*ni)+i] - avgconcentration_[(k*nj*ni)+(j*ni)+i];

                                sumd += abs(resor);
                                sumc += abs(avgconcentration_[(k*nj*ni)+(j*ni)+i]);
                    }//end for(i)
                 }//end for(j)
            }//end for(k)
 
            resorc_ = sumd/(sumc + small_);
            break;
        }//end case(6)
    }//end switch
}//end calculateResiduals

//////////////////////////////////////////////////////
//	        bareplateMaterialProps	            //
//////////////////////////////////////////////////////
void
CFDSolverManager::bareplateMaterialProps()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;
    double tliquid = domainMgr_->tliquid_;
    double tsolid = domainMgr_->tsolid_;
    double ystart = domainMgr_->ystart_;
    double alasrb = domainMgr_->alasrb_;
    double viscos = domainMgr_->viscos_;
    double dens = domainMgr_->dens_;
    double denl = domainMgr_->denl_;
    double cp[3] = {domainMgr_->acpa_,
                    domainMgr_->acpb_,
                    domainMgr_->acpl_};
    double cond[5] = {domainMgr_->thconsa_,
                      domainMgr_->thconsb_,
                      domainMgr_->thconsc_,
                      domainMgr_->thconla_,
                      domainMgr_->thconlb_};

    // define some local variables
    double diffs, diffl;

    #pragma omp parallel for private(diffs, diffl)
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        //----temperature is larger than liquidus 
        diffl = (cond[3]*temp_[loc] + cond[4])/cp[2]; 
        diff_[loc] = diffl;
        rho_[loc] = denl;
        vis_[loc] = viscos;            
        if(temp_[loc] >= tliquid) 
            continue;         
        
        //----temperature is less than solidus
        diffs = (cond[0]*(temp_[loc]*temp_[loc]) + 
                cond[1]*temp_[loc] + cond[2])/(cp[0]*temp_[loc] + cp[1]);
        diff_[loc] = diffs; 
        rho_[loc] = dens; 
        vis_[loc] = 1e10;
        if(temp_[loc] <= tsolid) 
            continue;  

        //----temperature is between liquidus and solidus           
        diff_[loc] = fracl_[loc]*diffl + (1.0 - fracl_[loc])*diffs;
        rho_[loc] = fracl_[loc]*denl + (1.0 - fracl_[loc])*dens;
        vis_[loc] = viscos;
    }//end for(loc)
}//end bareplateMaterialProps

//////////////////////////////////////////////////////
//	        powderMaterialProps	            //
//////////////////////////////////////////////////////
void
CFDSolverManager::powderMaterialProps()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    double tliquid = domainMgr_->tliquid_;
    double tsolid = domainMgr_->tsolid_;
    double ystart = domainMgr_->ystart_;
    double alasrb = domainMgr_->alasrb_;
    double viscos = domainMgr_->viscos_;
    double dens = domainMgr_->dens_;
    double denl = domainMgr_->denl_;
    double cp[3] = {domainMgr_->acpa_,
                    domainMgr_->acpb_,
                    domainMgr_->acpl_};
    double cond[5] = {domainMgr_->thconsa_,
                      domainMgr_->thconsb_,
                      domainMgr_->thconsc_,
                      domainMgr_->thconla_,
                      domainMgr_->thconlb_};

    // material properties for powder
    double powderden = domainMgr_->powderden_;
    double powdercp[2] = {domainMgr_->powdercpa_,
                          domainMgr_->powdercpb_};
    double powdercond[2] = {domainMgr_->powderthcona_,
                            domainMgr_->powderthconb_};

    // define some local variables
    double diffs, diffl, cps, conds, diffp, rhos;

    #pragma omp parallel for private(diffs, diffl, cps, conds, diffp, rhos)
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        //----temperature is larger than liquidus 
        diffl = (cond[3]*temp_[loc] + cond[4])/cp[2]; 
        diff_[loc] = diffl;
        rho_[loc] = denl;
        vis_[loc] = viscos;            
        if(temp_[loc] >= tliquid) 
            continue;         
        
        //----temperature is less than solidus
        // *changes made to account for effective powder properties*
        // calculate specific heat and conducitivty of the base-plate solid
        cps = cp[0]*temp_[loc] + cp[1];
        conds = cond[0]*(temp_[loc]*temp_[loc]) + 
                cond[1]*temp_[loc] + cond[2];

        // calculate diffusivity of the powder
        diffp = (powdercond[0]*temp_[loc] + powdercond[1])/
                (powdercp[0]*temp_[loc] + powdercp[1]);

        // calculate effective diffusivty and density of the material using linear approximation
        diffs = (1.0 - csfrac_[loc])*diffp + csfrac_[loc]*conds/cps;
        rhos = (1.0 - csfrac_[loc])*powderden + csfrac_[loc]*dens;

        // apply effective material properties
        diff_[loc] = diffs; 
        rho_[loc] = rhos; 
        vis_[loc] = 1e10;
        if(temp_[loc] <= tsolid) 
            continue;  

        //----temperature is between liquidus and solidus           
        diff_[loc] = fracl_[loc]*diffl + (1.0 - fracl_[loc])*diffs;
        rho_[loc] = fracl_[loc]*denl + (1.0 - fracl_[loc])*rhos;
        vis_[loc] = viscos;
    }//end for(loc)
}//end powderMaterialProps

//////////////////////////////////////////////////////
//		enhanceConvergenceSpeed             //
//////////////////////////////////////////////////////
void 
CFDSolverManager::enhanceConvergenceSpeed(int &ivar)
{
    // first check if this is needed
    double *phi;
    if(ivar == 5)
    {
        phi = &enthalpy_[0];
    }
    else if(ivar == 6)
    {
        phi = &avgconcentration_[0];
    }
    else
    {
        return;
    }//end if(ivar)
     
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    int nim2 = domainMgr_->nim2_;
    
    // define local variables
    int i, j, k;
    double denom;
    vector <double> bl(ni,0.0), blp(ni,0.0), blm(ni,0.0), blc(ni,0.0), 
                    delphi(ni,0.0), pib(ni,0.0), qib(ni,0.0);


    // in order to do reduction on vectors in C++, it must first be declared by user
    #pragma omp declare reduction(vec_double_plus : vector<double> : \
                                    std::transform(omp_out.begin(), omp_out.end(), \
                                    omp_in.begin(), omp_out.begin(), plus<double>())) \
                                    initializer(omp_priv = omp_orig)

    #pragma omp parallel for reduction(vec_double_plus: bl, blp, blm, blc) collapse(3) 
    for (k=1; k<nkm1; k++)
    {
        for (j=1; j<njm1; j++)
        {
            for (i=1; i<nim1; i++)
            {
                bl[i] += ap_[(k*nj*ni)+(j*ni)+i] - an_[(k*nj*ni)+(j*ni)+i] - as_[(k*nj*ni)+(j*ni)+i] - 
                         at_[(k*nj*ni)+(j*ni)+i] - ab_[(k*nj*ni)+(j*ni)+i];

                blp[i] += ae_[(k*nj*ni)+(j*ni)+i];
                blm[i] += aw_[(k*nj*ni)+(j*ni)+i];

                blc[i] += ae_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+(j*ni)+i+1] + 
                          aw_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+(j*ni)+i-1] + 
                          an_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j+1)*ni)+i] + 
                          as_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j-1)*ni)+i] + 
                          at_[(k*nj*ni)+(j*ni)+i]*phi[((k+1)*nj*ni)+(j*ni)+i] + 
                          ab_[(k*nj*ni)+(j*ni)+i]*phi[((k-1)*nj*ni)+(j*ni)+i] + 
                          su_[(k*nj*ni)+(j*ni)+i] - ap_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+(j*ni)+i];
            }//end for(i)
        }//end for(j)
    }//end for(k)
    
    pib[1] = blp[1]/bl[1];
    qib[1] = blc[1]/bl[1];

    for (i=2; i<nim1; i++)
    {
        denom = bl[i]-blm[i]*pib[i-1];
        pib[i] = blp[i]/denom;
        qib[i] = (blc[i] + blm[i]*qib[i-1])/denom;

    }//end for(i)
    delphi[nim2] = qib[nim2];

    for (i = nim2-1; i>0; i--) 
    {
        delphi[i] = pib[i]*delphi[i+1] + qib[i];
    }//end for(i)

    #pragma omp parallel for collapse(3)
    for (k=1; k<nkm1; k++)
    {
        for (j=1; j<njm1; j++)
        {
            for (i=1; i<nim1; i++)
            {
                phi[(k*nj*ni)+(j*ni)+i] += delphi[i];
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end enhanceConvergenceSpeed

//////////////////////////////////////////////////////
//		 solveLiquidDomain                  //
//////////////////////////////////////////////////////
void 
CFDSolverManager::solveLiquidDomain(int &ivar)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm2 = domainMgr_->nkm2_;

    // grab some stuff from the bc
    int istat = bcMgr_->istat_;
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;

    // define local variables
    int i, j, k, ksweep, jsweep, isweep;
    double d, denom;
    double *phi;
    double pr[ni], qr[ni];
    
    switch(ivar)
    {
        case 1:
        {
            phi = &uvel_[0];
            break;
        }//end case(1)
        case 2:
        {
            phi = &vvel_[0];
            break;
        }//end case(2)
        case 3:
        {
            phi = &wvel_[0];
            break;
        }//end case(3)
        case 4:
        {
            phi = &pp_[0];
            break;
        }//end case(4)
    }//end switch

    // TDMA
    for (ksweep=0; ksweep<2; ksweep++)
    {
        for (k=nkm2; k>=kstat; k--)
        {
            for (jsweep=0; jsweep<2; jsweep++)
            {
                #pragma omp parallel for private(pr, qr, i, d, denom)
                for (j=jstat; j<=jend; j++)
                {
                    i = istat;
                    pr[i] = 0.0;
                    qr[i] = phi[(k*nj*ni)+(j*ni)+i];
    
                    for (i=istatp1; i<=iendm1; i++)
                    {
                        d = at_[(k*nj*ni)+(j*ni)+i]*phi[((k+1)*nj*ni)+(j*ni)+i] + 
                            ab_[(k*nj*ni)+(j*ni)+i]*phi[((k-1)*nj*ni)+(j*ni)+i] +
                            an_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j+1)*ni)+i] + 
                            as_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j-1)*ni)+i] +
                            su_[(k*nj*ni)+(j*ni)+i];

                        denom = ap_[(k*nj*ni)+(j*ni)+i] - aw_[(k*nj*ni)+(j*ni)+i]*pr[i-1];
    
                        if(abs(denom) <= 1e-12) 
                            denom += small_;     //avoid divide zero
    
                        pr[i] = ae_[(k*nj*ni)+(j*ni)+i]/denom;
                        qr[i] = (d + aw_[(k*nj*ni)+(j*ni)+i]*qr[i-1])/denom;
                    }//end for(i)
                    
                    // back
                    for (i=iendm1; i>=istatp1; i--)
                    {
                        phi[(k*nj*ni)+(j*ni)+i] = pr[i]*phi[(k*nj*ni)+(j*ni)+i+1] + qr[i];
                    }//end for(i)
                }//end for(j)
            }//end for(jsweep)
        }//end for(k)
    }//end for(ksweep)
}//solveLiquidDomain

//////////////////////////////////////////////////////
//		 solveEntireDomain                  //
//////////////////////////////////////////////////////
void 
CFDSolverManager::solveEntireDomain(int &ivar)
{
    // first check if this is needed
    double *phi;
    if(ivar == 5)
    {
        phi = &enthalpy_[0];
    }
    else if(ivar == 6)
    {
        phi = &avgconcentration_[0];
    }
    else
    {
        return;
    }//end if(ivar)
    
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nim2 = domainMgr_->nim2_;
    int nkm2 = domainMgr_->nkm2_;

    // define local variables
    int i, j, k, ksweep, jsweep, isweep;
    double d, denom;
    double pr[ni], qr[ni];
    
    // TDMA
    for (ksweep=0; ksweep<2; ksweep++) //scan k direction twice
    {
        for (k=nkm2; k>0; k--)
        {
            for (jsweep=0; jsweep<2; jsweep++) //scan j direction twice
            {
                #pragma omp parallel for private(pr, qr, i, d, denom)
                for (j=1; j<njm1; j++)
                {
                    pr[0] = 0.0;
                    qr[0] = phi[(k*nj*ni)+(j*ni)+0];
    
                    for (i=1; i<nim1; i++)
                    {
                        d = at_[(k*nj*ni)+(j*ni)+i]*phi[((k+1)*nj*ni)+(j*ni)+i] + 
                            ab_[(k*nj*ni)+(j*ni)+i]*phi[((k-1)*nj*ni)+(j*ni)+i] + 
                            an_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j+1)*ni)+i] + 
                            as_[(k*nj*ni)+(j*ni)+i]*phi[(k*nj*ni)+((j-1)*ni)+i] + 
                            su_[(k*nj*ni)+(j*ni)+i];

                        denom = ap_[(k*nj*ni)+(j*ni)+i] - aw_[(k*nj*ni)+(j*ni)+i]*pr[i-1];
                        pr[i] = ae_[(k*nj*ni)+(j*ni)+i]/denom;
                        qr[i] = (d + aw_[(k*nj*ni)+(j*ni)+i]*qr[i-1])/denom;
                    }//end for(i)
                     
                    // back 
                    for (i=nim2; i>0; i--)
                    {
                        phi[(k*nj*ni)+(j*ni)+i]=pr[i]*phi[(k*nj*ni)+(j*ni)+i+1] + qr[i];
                    }//end for(i)
                }//end for(j)
            }//end for(jsweep)
        }//end for(k)
    }//end for(ksweep)
}//end solveEntireDomain

//////////////////////////////////////////////////////
//		 convertEnthalpyToTemp              //
//////////////////////////////////////////////////////
void 
CFDSolverManager::convertEnthalpyToTemp()
{
    // grab some stuff from the domain manager
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    double acpa = domainMgr_->acpa_;
    double acpb = domainMgr_->acpb_;
    double acpl = domainMgr_->acpl_;
    double hlcal = domainMgr_->hlcal_;
    double hsmelt = domainMgr_->hsmelt_;
    double tsolid = domainMgr_->tsolid_;
    double tliquid = domainMgr_->tliquid_;

    // define some local variables
    double deltemp = tliquid - tsolid;

    // convert the enthalpy into temperature
    #pragma omp parallel for
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        if(enthalpy_[loc] >= hlcal)
        {
            fracl_[loc] = 1.0;
            temp_[loc] = (enthalpy_[loc]-hlcal)/acpl + tliquid;
        }
        else if(enthalpy_[loc] <= hsmelt)
        {
            fracl_[loc] = 0.0;
            //temp_[loc] = tsolid - (hsmelt-enthalpy_[loc])/acp;
            temp_[loc] = (sqrt(acpb*acpb + 2*acpa*enthalpy_[loc]) - acpb)/acpa;
        }
        else
        {
            // mushy zone has linear assumption
            fracl_[loc] = (enthalpy_[loc] - hsmelt)/(hlcal - hsmelt);
            temp_[loc] = deltemp*fracl_[loc] + tsolid;
        }//end if
    }//end for(loc)
}//end convertEnthalpyToTemp


//////////////////////////////////////////////////////
//		 cleanVelocities                    //
//////////////////////////////////////////////////////
void 
CFDSolverManager::cleanVelocities()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;
    double tsolid = domainMgr_->tsolid_;
    
    // grab some stuff from the bc
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;
    
    // define some local variables
    int i, j, k;
    double tulc, tvlc, twlc;
    
    #pragma omp parallel for private(tulc, tvlc, twlc) collapse(3)
    for (k=kstat; k<nkm1; k++)
    {
        for (j=jstat; j<=jend; j++)
        {
            for (i=istatp1; i<=iendm1; i++)
            {
                tulc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+(j*ni)+i+1]);
                tvlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+((j+1)*ni)+i]);
                twlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[((k+1)*nj*ni)+(j*ni)+i]);

                if(tulc <= tsolid) 
                    uvel_[(k*nj*ni)+(j*ni)+i+1] = 0.0;

                if(tvlc <= tsolid) 
                    vvel_[(k*nj*ni)+((j+1)*ni)+i] = 0.0;

                if(twlc <= tsolid) 
                    wvel_[((k+1)*nj*ni)+(j*ni)+i] = 0.0;
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end cleanVelocities 

//////////////////////////////////////////////////////
//		updateSolutions                     //
//////////////////////////////////////////////////////
void 
CFDSolverManager::updateSolutions()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    double tsolid = domainMgr_->tsolid_;

    // grab some stuff from the BC
    double trackindx = (double)bcMgr_-> trackindx_;

    // push back the velocities
    #pragma omp parallel for
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        // zero out velocities in solid
        if (temp_[loc] <= tsolid)
        {
            uvel_[loc] = 0.0;
            vvel_[loc] = 0.0;
            wvel_[loc] = 0.0;
        }//end if

        // reserve previous time values
        unot_[loc] = uvel_[loc];
        vnot_[loc] = vvel_[loc];
        wnot_[loc] = wvel_[loc];
        hnot_[loc] = enthalpy_[loc];
        fraclnot_[loc] = fracl_[loc];

        if(meshObj_->powderbed_)
        {
            csfrac_[loc] = max(csfrac_[loc],fracl_[loc]);
            solfrac_[loc] = max(solfrac_[loc],fracl_[loc]*trackindx);
        }//end if
        if(meshObj_->species_)
        {
            avgconcentrationnot_[loc] = avgconcentration_[loc];
            concentration_[loc] = avgconcentration_[loc]/
                               (fracl_[loc] + domainMgr_->kp_*(1.0 - fracl_[loc]));
        }//end if
    }//end for(loc)
}//end updateSolutions

//////////////////////////////////////////////////////
//	        correctPressure	            	    //
//////////////////////////////////////////////////////
void
CFDSolverManager::correctPressure(int &ivar)
{
    // exit if pp is not the integration variable
    if(ivar != 4) 
        return;
     
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;
    double tsolid = domainMgr_->tsolid_;
    double urfp = domainMgr_->urfp_;
     
    // grab some stuff from the bc
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    int kstat = bcMgr_->kstat_;

    // define some local variables 
    int i, j, k;
    double tulc, tvlc, twlc;

    // calculate the pressure correction (for conservation of mass)
    #pragma omp parallel for private(tulc, tvlc, twlc) collapse(3)
    for (k=kstat; k<nkm1; k++)
    {
        for (j=jstat; j<=jend; j++)
        {
            for (i=istatp1; i<=iendm1; i++)
            {
                tulc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+(j*ni)+i-1]);
                if(tulc > tsolid)
                    uvel_[(k*nj*ni)+(j*ni)+i] += dux_[(k*nj*ni)+(j*ni)+i]*
                                                (pp_[(k*nj*ni)+(j*ni)+i-1] - pp_[(k*nj*ni)+(j*ni)+i]);

                tvlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[(k*nj*ni)+((j-1)*ni)+i]);
                if(tvlc > tsolid)
                    vvel_[(k*nj*ni)+(j*ni)+i] += dvy_[(k*nj*ni)+(j*ni)+i]*
                                                (pp_[(k*nj*ni)+((j-1)*ni)+i] - pp_[(k*nj*ni)+(j*ni)+i]);

                twlc = min(temp_[(k*nj*ni)+(j*ni)+i], temp_[((k-1)*nj*ni)+(j*ni)+i]);
                if(twlc > tsolid)
                    wvel_[(k*nj*ni)+(j*ni)+i] += dwz_[(k*nj*ni)+(j*ni)+i]*
                                                (pp_[((k-1)*nj*ni)+(j*ni)+i] - pp_[(k*nj*ni)+(j*ni)+i]);
            }//ennd for(i)
        }//end for(j)
    }//end for(k)

    for (k=kstat; k<nkm1; k++)
    {
        for (j=jstat; j<=jend; j++)
        {
            for (i=istatp1; i<=iendm1; i++)
            {
                if(temp_[(k*nj*ni)+(j*ni)+i] > tsolid)
                {
                    pressure_[(k*nj*ni)+(j*ni)+i] += urfp*pp_[(k*nj*ni)+(j*ni)+i];
                    pp_[(k*nj*ni)+(j*ni)+i] = 0.0;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end correctPressure

//////////////////////////////////////////////////////
//	        getSolidificationParameters    	    //
//////////////////////////////////////////////////////
void
CFDSolverManager::getSolidificationParameters()
{
    // exit if there is no liquid region
    if(bcMgr_->tpeak_ <= domainMgr_->tsolid_)
        return;
    
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;
    double scanvel = domainMgr_->scanvel_;
     
    // define some local variables
    int i, j, k;
    double gradx, grady, gradz;
    
    // calculate temperature gradient and solidification growth rate
    for (k = 0; k<nkm1; k++)
    {
        for (j = 0; j<nj; j++)
        {
            for (i = 0; i<ni; i++)
            {
                if(temp_[(k*nj*ni)+(j*ni)+i] >= domainMgr_->tsolid_)
                {
                    gradx = (temp_[(k*nj*ni)+(j*ni)+i] - temp_[(k*nj*ni)+(j*ni)+i-2])/(x_[i] - x_[i-1])/2;
                    grady = (temp_[(k*nj*ni)+((j+1)*ni)+i-1] - temp_[(k*nj*ni)+((j-1)*ni)+i-1])/(y_[j] - y_[j-1])/2;
                    gradz = (temp_[((k+1)*nj*ni)+(j*ni)+i-1] - temp_[((k-1)*nj*ni)+(j*ni)+i-1])/(z_[k] - z_[k-1])/2;
                    grad_[(k*nj*ni)+(j*ni)+i-1] = sqrt(gradx*gradx + grady*grady + gradz*gradz);
                    
                    //if(grad_[(k*nj*ni)+(j*ni)+i-1] >= 1e8) 
                    //    grad_[(k*nj*ni)+(j*ni)+i-1] = 1e8;
                    
                    rate_[(k*nj*ni)+(j*ni)+i-1] = scanvel*gradx/grad_[(k*nj*ni)+(j*ni)+i-1];
                    break;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end getSolidificationParameters

//////////////////////////////////////////////////////
//	        remapSolutions	                    //
//////////////////////////////////////////////////////
void
CFDSolverManager::remapSolutions(const vector <double> &oldx, const vector <double> &oldy, const vector <double> &oldz)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    double newXCoordinate[ni*nj*nk], newYCoordinate[ni*nj*nk], newZCoordinate[ni*nj*nk], oldXCoordinate[ni*nj*nk], oldYCoordinate[ni*nj*nk], oldZCoordinate[ni*nj*nk];
    int dimensionOfDomain[] = {ni, nj, nk}, lengthOfFlattenDimension = ni*nj*nk, locationIndex[ni][nj][nk], ySortedlocationIndex[ni*nj*nk], location;
    // string timeFileOutput, filename;



    // mapping old and new coordinate from DomainManager coordinate index

    for (int count=0; count<lengthOfFlattenDimension; count++)
    {
        oldXCoordinate[count] = oldx[this->domainMgr_->coordinateIndex[count][0]];
        oldYCoordinate[count] = oldy[this->domainMgr_->coordinateIndex[count][1]];
        oldZCoordinate[count] = oldz[this->domainMgr_->coordinateIndex[count][2]];

        newXCoordinate[count] = x_[this->domainMgr_->coordinateIndex[count][0]];
        newYCoordinate[count] = y_[this->domainMgr_->coordinateIndex[count][1]];
        newZCoordinate[count] = z_[this->domainMgr_->coordinateIndex[count][2]];

    }



// // writing file to disc for further analysis

//     timeFileOutput = "time="+to_string(this->domainMgr_->timet_);
//     filename = "TempBeforeRemap"+timeFileOutput+".txt";

//     ofstream myoldfile (filename);
//     if (myoldfile.is_open())
//     {
//         for (int loc = 0; loc<lengthOfFlattenDimension; loc++)
//             myoldfile << temp_[loc] <<"\n";      
//     }
    
//     myoldfile.close();
    
//     filename = "oldMapX"+timeFileOutput+".txt";

//     ofstream myoldXfile (filename);
//     if (myoldXfile.is_open())
//     {
//         for (int i = 0; i<lengthOfFlattenDimension; i++)
//             myoldXfile << oldXCoordinate[i] <<"\n";      
//     }

//     myoldXfile.close();

//     filename = "oldMapY"+timeFileOutput+".txt";

//     ofstream myoldYfile (filename);
//     if (myoldYfile.is_open())
//     {
//         for (int j = 0; j<lengthOfFlattenDimension; j++)
//             myoldYfile << oldYCoordinate[j] <<"\n";      
//     }

//     myoldYfile.close();

//     filename = "oldMapZ"+timeFileOutput+".txt";

//     ofstream myoldZfile (filename);
//     if (myoldZfile.is_open())
//     {
//         for (int k = 0; k<lengthOfFlattenDimension; k++)
//             myoldZfile << oldZCoordinate[k] <<"\n";      
//     }

//     myoldZfile.close();



    // loop over the solution variables and map them to the new mesh
    // #pragma omp parallel for
    // for (int loc = 0; loc<ni*nj*nk; loc++)
    // {
    //     int i = loc % ni;

    //     // currently the solutions are remapped by a simple linear interpolation
    //     enthalpy_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], enthalpy_[loc], enthalpy_[loc+meshXMovementDirection], x_[i]);
    //     hnot_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], hnot_[loc], hnot_[loc+meshXMovementDirection], x_[i]);
    //     temp_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], temp_[loc], temp_[loc+meshXMovementDirection], x_[i]);
    //     fracl_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], fracl_[loc], fracl_[loc+meshXMovementDirection], x_[i]);
    //     uvel_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], uvel_[loc], uvel_[loc+meshXMovementDirection], x_[i]);
    //     vvel_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], vvel_[loc], vvel_[loc+meshXMovementDirection], x_[i]);
    //     wvel_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], wvel_[loc], wvel_[loc+meshXMovementDirection], x_[i]);
    //     unot_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], unot_[loc], unot_[loc+meshXMovementDirection], x_[i]);
    //     vnot_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], vnot_[loc], vnot_[loc+meshXMovementDirection], x_[i]);
    //     wnot_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], wnot_[loc], wnot_[loc+meshXMovementDirection], x_[i]);
    //     pressure_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], pressure_[loc], pressure_[loc+meshXMovementDirection], x_[i]);
    //     pp_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], pp_[loc], pp_[loc+meshXMovementDirection], x_[i]);

            
    //     if(meshObj_->species_)
    //     {                     
    //         avgconcentration_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], avgconcentration_[loc], avgconcentration_[loc+meshXMovementDirection], x_[i]);
    //         avgconcentrationnot_[loc] = interpolate(oldx[i], oldx[i+meshXMovementDirection], avgconcentrationnot_[loc], avgconcentrationnot_[loc+meshXMovementDirection], x_[i]);

    //     }//end if
    // }//end for(loc)

// interpolation scheme for mesh movement along X-direction
    int loc = 0;
    // #pragma omp parallel for
    for (int k = 0; k<nk; k++)
    {
        for (int j = 0; j<nj-1; j++)
        {
            for (int i = 1; i<ni-1; i++)    // i count is starting from 1 instead of 0 because if mesh is moving in -ve direction, it till try to access -1 th element of vector
            {
                loc = i + j*ni + k*ni*nj;
                if (abs(oldXCoordinate[loc] - newXCoordinate[loc]) > 1e-6)
                // along y direction
                {
                    enthalpy_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], enthalpy_[loc], enthalpy_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    hnot_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], hnot_[loc], hnot_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    temp_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], temp_[loc], temp_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    fracl_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], fracl_[loc], fracl_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    uvel_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], uvel_[loc], uvel_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    vvel_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], vvel_[loc], vvel_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    wvel_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], wvel_[loc], wvel_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    unot_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], unot_[loc], unot_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    vnot_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], vnot_[loc], vnot_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    wnot_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], wnot_[loc], wnot_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    pressure_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], pressure_[loc], pressure_[loc+meshXMovementDirection], newXCoordinate[loc]);
                    pp_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], pp_[loc], pp_[loc+meshXMovementDirection], newXCoordinate[loc]);
                
                    if(meshObj_->species_)
                    {                     
                        avgconcentration_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], avgconcentration_[loc], avgconcentration_[loc+meshXMovementDirection], newXCoordinate[loc]);
                        avgconcentrationnot_[loc] = interpolate(oldXCoordinate[loc], oldXCoordinate[loc+meshXMovementDirection], avgconcentrationnot_[loc], avgconcentrationnot_[loc+meshXMovementDirection], newXCoordinate[loc]);

                    }//end if
                }
            }//end for(loc)
        }

    }



// interpolation scheme for mesh movement along Y-direction
    // int loc = 0;
    // #pragma omp parallel for
    for (int k = 0; k<nk; k++)
    {
        for (int j = 0; j<nj-1; j++)
        {
            for (int i = 1; i<ni-1; i++)    // i count is starting from 1 instead of 0 because if mesh is moving in -ve direction, it till try to access -1 th element of vector
            {
                loc = i + j*ni + k*ni*nj;
                if (abs(oldYCoordinate[loc] - newYCoordinate[loc]) > 1e-6)
                // along y direction
                {
                    enthalpy_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], enthalpy_[loc], enthalpy_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    hnot_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], hnot_[loc], hnot_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    temp_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], temp_[loc], temp_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    fracl_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], fracl_[loc], fracl_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    uvel_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], uvel_[loc], uvel_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    vvel_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], vvel_[loc], vvel_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    wvel_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], wvel_[loc], wvel_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    unot_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], unot_[loc], unot_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    vnot_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], vnot_[loc], vnot_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    wnot_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], wnot_[loc], wnot_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    pressure_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], pressure_[loc], pressure_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                    pp_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], pp_[loc], pp_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                
                    if(meshObj_->species_)
                    {                     
                        avgconcentration_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], avgconcentration_[loc], avgconcentration_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);
                        avgconcentrationnot_[loc] = interpolate(oldYCoordinate[loc], oldYCoordinate[loc+meshYMovementDirection*ni], avgconcentrationnot_[loc], avgconcentrationnot_[loc+meshYMovementDirection*ni], newYCoordinate[loc]);

                    }//end if
                }
            }//end for(loc)
        }

    }
    









// // second set of file writing to the disc for analysis
//     filename = "TempAfterRemap"+timeFileOutput+".txt";
    
//     ofstream mynewfile (filename);
//     if (mynewfile.is_open())
//     {
//         for (int loc = 0; loc<lengthOfFlattenDimension; loc++)
//             mynewfile << temp_[loc] <<"\n";
//     }

//     mynewfile.close();

//     filename = "newMapX_"+timeFileOutput+".txt";

//     ofstream myX_file (filename);
//     if (myX_file.is_open())
//     {
//         for (int i = 0; i<lengthOfFlattenDimension; i++)
//             myX_file << newXCoordinate[i] <<"\n";      
//     }

//     myX_file.close();

//     filename = "newMapY_"+timeFileOutput+".txt";

//     ofstream myY_file (filename);
//     if (myY_file.is_open())
//     {
//         for (int j = 0; j<lengthOfFlattenDimension; j++)
//             myY_file << newYCoordinate[j] <<"\n";      
//     }

//     myY_file.close();


//     filename = "newMapZ_"+timeFileOutput+".txt";

//     ofstream myZ_file (filename);
//     if (myZ_file.is_open())
//     {
//         for (int k = 0; k<lengthOfFlattenDimension; k++)
//             myZ_file << newZCoordinate[k] <<"\n";      
//     }

//     myZ_file.close();    

// end of file writing for the new values


}//end remapSolutions

//////////////////////////////////////////////////////
//	        linear interpolation  	                //
//////////////////////////////////////////////////////

double 
CFDSolverManager::interpolate(double x0, double x1, double y0, double y1, double x)
{
    double interpolatedValue, yDifference, xDifference;

    yDifference = (y1 - y0);
    xDifference = (x1 - x0);
    interpolatedValue = y0 + (x - x0)*(yDifference/xDifference);
    return interpolatedValue;
}

//////////////////////////////////////////////////////
//	        discretPower         	            //
//////////////////////////////////////////////////////
inline double 
CFDSolverManager::discretPower(double pe) 
{
    // faster alternative to the pow() function 
    double base = 1.0 - 0.1*pe;
    double out = base;
    for(int i=1; i<5; i++)
        out *= base;

   return max(0.0, out);
}//end discretPower

//////////////////////////////////////////////////////
//	        discretUpwind         	            //
//////////////////////////////////////////////////////
inline double 
CFDSolverManager::discretUpwind(double pe) 
{
    // upwind only uses upstream data 
   return 1.0;
}//end discretUpwind

//////////////////////////////////////////////////////
//	        discretExponential     	            //
//////////////////////////////////////////////////////
inline double 
CFDSolverManager::discretExponential(double pe) 
{
    // exponential is closer to "exact solution" in <numerical heat transfer>
   return pe/(exp(pe) - 1.0);
}//end discretExponential
