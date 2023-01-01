// standard headers
#include <vector>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <omp.h>

// user-defined headers
#include <FreeSurfaceManager.h>

///////////////////////////////////////////
//		Constructor		 //
///////////////////////////////////////////
FreeSurfaceManager::FreeSurfaceManager(Mesh *meshObj,
                                       DomainManager *domainMgr,
                                       SolutionVariableManager *solVarMgr,
                                       BoundCondManager *bcMgr)
                                       : meshObj_(meshObj),
                                         domainMgr_(domainMgr),
                                         solVarMgr_(solVarMgr),
                                         bcMgr_(bcMgr)
{
}

//////////////////////////////////////////////////////
//		initializeFreeSurface		    //
//////////////////////////////////////////////////////
void 
FreeSurfaceManager::initializeFreeSurface()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;

    // set up appropipate pointers to cfd solver manager
    temp_ = solVarMgr_->mapSolutionVars_["temp"]; 
    areaij_ = solVarMgr_->mapGeometryVars_["areaij"];       
    xu_ = solVarMgr_->mapGeometryVars_["xu"];       
    yv_ = solVarMgr_->mapGeometryVars_["yv"];       
    x_ = solVarMgr_->mapGeometryVars_["x"];       
    y_ = solVarMgr_->mapGeometryVars_["y"];       
    z_ = solVarMgr_->mapGeometryVars_["z"];       
    
    // initialize local arrays
    for (int loc = 0; loc<ni*nj*nk; loc++)
    {
        zr_[loc] = z_[loc / (ni*nj)]; 
        tempr_[loc] = temp_[loc];

        if(loc < ni*nj)
        {
            fi_[loc] = 0.0; 
        }//end if
    }//end for(loc)

    // map adjusted quanities to the viewing window
    solVarMgr_->mapSolutionVars_["tempr"] = &tempr_[0];     
    solVarMgr_->mapSolutionVars_["zr"] = &zr_[0];     
}//end initializeFreeSurface

//////////////////////////////////////////////////////
//		execute		               	    //
//////////////////////////////////////////////////////
void
FreeSurfaceManager::execute()
{
    // grab some stuff from the domain and mesh
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    double tsolid = domainMgr_->tsolid_;
    double timet = domainMgr_->timet_;
    double delt = meshObj_->delt_;

    // initialize bi-section solver parameters based on input file
    int nitermax = (int)domainMgr_->maxiterbisect_;
    double alam1 = domainMgr_->initlam1_;
    double alam2 = domainMgr_->initlam2_;
    double tol = domainMgr_->tolbisect_;
    
    // define some local variables
    int i, j, k, niter;
    double alamopt, dx, xmid, fmid, f;
    double factor, vd;
    
    // calculate the intial time to start tracking free surface
    double startdeltn = 0.0;
    double melttime = delt*(1.0 + startdeltn);
    
    // calculate pressure
    fill(&apdis_[0], &apdis_[0]+(nj*ni), 0.0);
    for(j=0; j<nj; j++)
    {
        for(i=0; i<ni; i++)
        {
            //apdis_(i,j)=0.5e4*exp(-((beam_pos-x_[i])**2+(beam_posy-y_[j])**2)/(2.0*alasrb**2))
            apdis_[(j*ni)+i] = 0.0;
        }//end for(i)
    }//end for(j)

    // calculate optimum lagrangian parameter and free surface
    double vdp = 2e-9;
    if(vdp == 0.0 && *max_element(&apdis_[0], &apdis_[0]+(ni*nj)) == 0.0) 
        return;        
    
    // melttime is the inital time the substract is melted
    // vd is volume of the additional droplet
    if(bcMgr_->tpeak_ <= tsolid)
    {
        vd = 0.0;
        melttime = 0.0;
    }
    else
    {
        // volume addition here is adjustiable  
        vd = vdp*(timet - melttime - startdeltn*delt);

        // re-adjust volume addition if the laser is off
        if(bcMgr_->powerindicator_ < 1)   
            vd = vdp*(domainMgr_->lasertime_ - melttime);
    }//end if
    
    // caclulate f and fmid using alam1 and alam2
    f = getDeltaVol(alam1) + vd;
    fmid = getDeltaVol(alam2) + vd;
    
    if(f*fmid >= 0.0 && f != vd) 
        cout << "WARNING: lambda should be a bracked interval \n";
    
    // choose initial conditions for Bisection Method
    if(f < 0.0)
    {
        alamopt = alam1;
        dx = alam2 - alam1;
    }
    else
    {
        alamopt = alam2;
        dx = alam1 - alam2;
    }//end if

    // search for satisfactory lambda by using Bisection Method
    for (niter=0; niter<nitermax; niter++)
    {
        dx *= 0.5;
        xmid = alamopt + dx;
        fmid = getDeltaVol(xmid) + vd;
    
        if(fmid < 0.0)
            alamopt = xmid;
    
        if(abs(dx) < 1.0e-5 || abs(fmid) < tol && niter != 0)
            break;
    }//end for(niter)
    
    if (niter == nitermax-1)
        cout << "too many bisections\n";

    // calculate free surface profile
    for (k=0; k<nk; k++)
    {
        for (j=0; j<nj; j++)
        {
            for (i=0; i<ni; i++)
            {
                factor = fi_[(j*ni)+i]/z_[nkm1];
                zr_[(k*nj*ni)+(j*ni)+i] = z_[k]*(1.0 - factor);
            
                if(zr_[(k*nj*ni)+(j*ni)+i] > z_[nkm1])
                    continue;

                if(zr_[(k*nj*ni)+(j*ni)+i] == z_[nkm1]) 
                    tempr_[(k*nj*ni)+(j*ni)+i] = temp_[(nkm1*nj*ni)+(j*ni)+i];

                for (int l=nkm1-1; l>=0; l--)
                {
                    if(z_[l] <= zr_[(k*nj*ni)+(j*ni)+i])
                    {
                        tempr_[(k*nj*ni)+(j*ni)+i] = (zr_[(k*nj*ni)+(j*ni)+i] - z_[l])/(z_[l+1] - z_[l]) * 
                                                     (temp_[((l+1)*nj*ni)+(j*ni)+i]-temp_[(l*nj*ni)+(j*ni)+i]) + 
                                                     temp_[(l*nj*ni)+(j*ni)+i];
                        break;    
                    }//end if
                }//end for(l)   
            }//end for(i)
        }//end for(j)
    }//end for(k)

    // output some solver data for the Bisection-Method
    freesurfFile_.open(meshObj_->outFileName_ + "_FreeSurfaceMonitor.txt", ios::out | ios::app);
    freesurfFile_ << "====================================================================================================\n";
    freesurfFile_ << "  time         lambda        res_bisect    min_fi         max_fi        droplet_volume   iter_bisect\n";
    freesurfFile_ << setprecision(4) << scientific << setw(12) << timet << setw(14) << xmid << setw(14) << fmid
                                     << setw(14) << *min_element(&fi_[0],&fi_[0]+ni*nj) 
                                     << setw(14) << *max_element(&fi_[0],&fi_[0]+ni*nj) << setw(14) << vd 
                                     << setw(9) << niter << "\n";
    freesurfFile_ << "====================================================================================================\n";
    freesurfFile_.close();
}//end execute

//////////////////////////////////////////////////////
//		getDeltaVol            		    //
//////////////////////////////////////////////////////
double
FreeSurfaceManager::getDeltaVol(double &alambda)
{
    // This function in conjuction with energy balance
    //  routine define the function (delta Volume) 
        
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
   
    // calculate energy balence first
    ebal(alambda);

    // calculate change in volume
    double delv = 0.0;
    for (int j=1; j<domainMgr_->njm1_; j++)
    {
        for (int i=1; i<domainMgr_->nim1_; i++)
        {
            delv += areaij_[(j*ni)+i]*fi_[(j*ni)+i];
        }//end for(i)
    }//end for(j)

    return delv;
}//end getDeltaVol

//////////////////////////////////////////////////////
//		ebal	               		    //
//////////////////////////////////////////////////////
void
FreeSurfaceManager::ebal(double &alambda)
{
    // This function defines the free surface for given:
    //   lambda - lagrangian parameter
    //   gamma - surface tension coefficient

    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
     
    // initialize energy minization parameters based on input file
    int nitermax = (int)domainMgr_->maxiterenergy_;
    double c1 = domainMgr_->dens_ * meshObj_->grav_;
    double gamma = domainMgr_->surftens_;
    double urfree = domainMgr_->urfree_;
    double tol = domainMgr_->tolenergy_;

    // define some local variables
    int i, j, k;
    double c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, 
           diff1, dely, dely1, delyp, delym, 
           delx, delx1, delxp, delxm, dzdx, dzdy, dzdyp1, dzdym1, 
           d2zdxdy,term1, term2;

    // set up parameters for the iterative solve
    int niter = 0;
    double difmax = 0.0;

    // zero out the old configuration parameter
    for (int loc=0; loc < ni*nj; loc++)
        zzm_[loc] = 0.0;

    // iteratively solve for configuration parameter
    while(niter < nitermax)
    {
        niter++;
        
        // calculate coefficients: c1 ... c11
        #pragma omp parallel for private(dely,dely1,diff1,delyp,delym,delx,delx1,delxp,delxm,c2,dzdx,dzdy,dzdyp1,dzdym1,d2zdxdy,c3,c4,c5,c6,c7,c8,c9,c10,c11,term1,term2) \
                                        reduction(max:difmax) collapse(2)
        for (j=1; j<njm1; j++)
        {
            for (i=1; i<nim1; i++)
            {
                // if temp is less than solidus, profile is not calculated
                if(temp_[(nkm1*nj*ni)+(j*ni)+i] < domainMgr_->tsolid_) 
                    continue;   
                 
                // the differences in y are calculated here to parallelize loops with OpenMP
                dely = yv_[j+1] - yv_[j-1];
                dely1 = y_[j+1] - y_[j-1];
                delyp = y_[j+1] - y_[j];
                delym = y_[j] - y_[j-1];

                // calculate differences in x
                delx = xu_[i+1] - xu_[i-1];
                delx1 = x_[i+1] - x_[i-1];
                delxp = x_[i+1] - x_[i];
                delxm = x_[i] - x_[i-1];

                // calculte the coeffeicents for the energy balence  
                c2 = apdis_[(j*ni)+i] + alambda;
                dzdx = (zzm_[(j*ni)+i+1] - zzm_[(j*ni)+i-1])/delx1;
                dzdy = (zzm_[((j+1)*ni)+i] - zzm_[((j-1)*ni)+i])/dely1;
                dzdyp1 = (zzm_[((j+1)*ni)+i+1] - zzm_[((j-1)*ni)+i+1])/dely1;
                dzdym1 = (zzm_[((j+1)*ni)+i-1] - zzm_[((j-1)*ni)+i-1])/dely1;
                d2zdxdy = (dzdyp1 - dzdym1)/delx1;
                c3 = pow(1.0 + dzdx*dzdx + dzdy*dzdy, 1.5);
                c4 = 1.0 + dzdy*dzdy;
                c5 = 1.0 + dzdx*dzdx;
                c6 = 2.0*dzdx*dzdy*d2zdxdy;
                c7 = zzm_[(j*ni)+i+1]/(delxp*delx) + zzm_[(j*ni)+i-1]/(delxm*delx);
                c8 = zzm_[((j+1)*ni)+i]/(delyp*dely) + zzm_[((j-1)*ni)+i]/(delym*dely);
                c9 = 1.0/(delxp*delx) + 1.0/(delxm*delx);
                c10 = 1.0/(delyp*dely) + 1.0/(delym*dely);
                c11 = gamma/c3;
                term1 = c2 + c11*(c4*c7 + c5*c8 - c6);
                term2 = c1 + c11*(c4*c9 + c5*c10);

                // update the free surface configuration parameter 
                fi_[(j*ni)+i] = urfree*term1/term2 + (1.0 - urfree)*zzm_[(j*ni)+i];
                
                // calculate error
                diff1 = abs((fi_[(j*ni)+i] - zzm_[(j*ni)+i])/zzm_[(j*ni)+i]);
                difmax = max(diff1,difmax);
            }//end for(i)
        }//end for(j)
        
        // check for convergence
        if(difmax < tol) 
            break;
        
        // push back configuration parameter
        for (int loc=0; loc<ni*nj; loc++)
            zzm_[loc] = fi_[loc];
    }//end while(niter)
}//end ebal
