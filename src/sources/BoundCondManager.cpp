// standard headers
#include <algorithm>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <numeric>

// user-defined headers
#include <BoundCondManager.h>
#include <PrintManager.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
BoundCondManager::BoundCondManager(Mesh *meshObj,
                                   DomainManager *domainMgr,
                                   SolutionVariableManager *solVarMgr)
                                   : meshObj_(meshObj),
                                     domainMgr_(domainMgr),
                                     solVarMgr_(solVarMgr)
{
}

//////////////////////////////////////////////////////
//		initializeBoundaries           	    //
//////////////////////////////////////////////////////
void
BoundCondManager::initializeBoundaries()
{
    // initialize local arrays
    fill(&heatin_[0], &heatin_[0]+(domainMgr_->nj_*domainMgr_->ni_), 0.0);
    
    // set up appropipate pointers to main cfd solver (solutions)
    pp_ = solVarMgr_->mapSolutionVars_["pp"]; 
    vis_ = solVarMgr_->mapSolutionVars_["vis"]; 
    rho_ = solVarMgr_->mapSolutionVars_["rho"]; 
    uvel_ = solVarMgr_->mapSolutionVars_["uvel"]; 
    vvel_ = solVarMgr_->mapSolutionVars_["vvel"]; 
    wvel_ = solVarMgr_->mapSolutionVars_["wvel"]; 
    temp_ = solVarMgr_->mapSolutionVars_["temp"]; 
    diff_ = solVarMgr_->mapSolutionVars_["diff"]; 
    hnot_ = solVarMgr_->mapSolutionVars_["hnot"]; 
    fracl_ = solVarMgr_->mapSolutionVars_["fracl"]; 
    fraclnot_ = solVarMgr_->mapSolutionVars_["fraclnot"]; 
    enthalpy_ = solVarMgr_->mapSolutionVars_["enthalpy"]; 
    sourceinput_ = solVarMgr_->mapSolutionVars_["sourceinput"];
    if(meshObj_->species_)
        avgconcentration_ = solVarMgr_->mapSolutionVars_["avgconcentration"]; 

    // set up appropipate pointers to main cfd solver (geometry)
    volume_ = solVarMgr_->mapGeometryVars_["volume"];       
    areaij_ = solVarMgr_->mapGeometryVars_["areaij"];       
    areajk_ = solVarMgr_->mapGeometryVars_["areajk"];       
    areaik_ = solVarMgr_->mapGeometryVars_["areaik"];       
    dxpwinv_ = solVarMgr_->mapGeometryVars_["dxpwinv"];       
    dypsinv_ = solVarMgr_->mapGeometryVars_["dypsinv"];       
    dzpbinv_ = solVarMgr_->mapGeometryVars_["dzpbinv"];       
    fracx_ = solVarMgr_->mapGeometryVars_["fracx"];       
    fracy_ = solVarMgr_->mapGeometryVars_["fracy"];       
    x_ = solVarMgr_->mapGeometryVars_["x"];       
    y_ = solVarMgr_->mapGeometryVars_["y"];       
    z_ = solVarMgr_->mapGeometryVars_["z"];       
}//end initializeBoundaries

//////////////////////////////////////////////////////
//		updateLaser                  	    //
//////////////////////////////////////////////////////
void
BoundCondManager::updateLaser()
{
    // grab some stuff from the domain
    double ni = domainMgr_->ni_;
    double nj = domainMgr_->nj_;
    double alasfact = domainMgr_->alasfact_;

    // define some local variables
    int i, j, iout, jout, loc2d;
    double xloc, rb2, varlas, xdist, ydist, dist2;

    // determine position of the heat source
    if(domainMgr_->steady_) 
    {
        // for steady-state simulation, position of heat source is fixed
        beamposx_ = domainMgr_->xstart_;
    }
    else                   
    {        
        // for transient simulation, position of heat source is time-dependent
        getLaserCoords();
    }//end if

    iout = 0;            // reserve intial index of x[]
    jout = 0;            // reserve intial index of y[]

    xloc = beamposx_;    // inital x-coodinate value of laser beam center

    for (i=1; i<domainMgr_->nim1_; i++)
    {
        // x(i)>=xloc, exit loop, iout=i-1 at this time
        if (xloc <= x_[i])
            break;
        iout = i;
    }//end for(i)

    // judge which is closer to xloc, x[i] or x[i-1]
    if(abs(xloc - x_[iout+1]) < abs(xloc - x_[iout])) 
        iout++;      
    
    // reserve index of x(), which is nearest to xloc
    istart_ = iout;                                              

    // find the index j of y(), for which y(j)~beamposy
    for(j=1; j<domainMgr_->njm1_; j++)
    {
        // y(j)>=beamposy, exit loop, jout=j-1 at this time    
        if (beamposy_ <= y_[j])  
            break;
        jout = j;
    }//end for(j)

    // judge which is closer to yloc, y[j] or y[j-1]
    if(abs(beamposy_-y_[jout+1]) < abs(beamposy_-y_[jout])) 
        jout++;

    // reserve index of y[], which is nearest to beamposy
    jstart_ = jout;

    // energy deposition from the heat source 
    heatinlaser_ = 0.0;                                 // total heat deposition
    rb2 = domainMgr_->alasrb_*domainMgr_->alasrb_;      // square of beam radius
    varlas = domainMgr_->alaspow_*domainMgr_->alaseta_; // effective laser power
    peakhin_ = alasfact*varlas/(M_PI*rb2);              // peak power density

    for (j=0; j<nj; j++)
    {
        // distance between y-value and pre-position of y
        ydist = beamposy_ - y_[j];
        for (i=0; i<ni; i++)
        {
            // mapping from 1D flatten array to 2D matrix 
            loc2d = (ni*j)+i;
            
            // distance between x_-value and laser beam center
            xdist = beamposx_ - x_[i];

            // total distance from laser center 
            dist2 = (xdist*xdist) + (ydist*ydist);

            // gaussian distribution
            heatin_[loc2d] = peakhin_*exp(-alasfact*dist2/rb2)*(double)powerindicator_;

            // if laser factor is nearly zero, become uniform distribution
            if(alasfact < 1.0e-6 && dist2 < rb2) 
                heatin_[loc2d] = varlas/(M_PI*rb2)*(double)powerindicator_;

             // calculate total heat flux
            heatinlaser_ += areaij_[loc2d]*heatin_[loc2d];
        }//end for(i) 
    }//end for(j)

    //----initial dimension of the pool
    imin_ = istart_ - 2;
    imax_ = istart_ + 2;

    jmin_ = jstart_ - 2;
    jmax_ = jstart_ + 2;

    kmin_ = domainMgr_->nkm2_ - 3;
    kmax_ = domainMgr_->nkm2_;
}//end updateLaser

//////////////////////////////////////////////////////
//	        getLaserCoords                      //
//////////////////////////////////////////////////////
void
BoundCondManager::getLaserCoords()
{

    // determine if laser is single-track or not 
    if(!meshObj_->udtoolpath_)
    {
        // save off old toolpath
        oldbeamposx_ = beamposx_;

        // calculate x position
        if(domainMgr_->timet_ <= domainMgr_->lasertime_)
        {
            beamposx_ = domainMgr_->xstart_ + domainMgr_->timet_*domainMgr_->scanvel_;
            if(init_)
            {
                oldbeamposx_ = beamposx_;
                init_ = false;
            }//end if
        }
        else 
        {
            beamposx_ = domainMgr_->xstart_ + domainMgr_->lasertime_*domainMgr_->scanvel_;
            powerindicator_ = 0;
        }
        // define the position of ystart
        beamposy_ = domainMgr_->ystart_; 
    }
    else
    {
        // define some variables
        double tx, ty, tz, deltatx, deltaty, deltatz;
        double laserTimeNp1, laserTimeN;
        double small = 1.0e-8;

        // save off old toolpath
        oldbeamposx_ = beamposx_;
        oldbeamposy_ = beamposy_;
        oldbeamposz_ = beamposz_;
        int oldpowerindicator = powerindicator_;

        // loop over the file data
        for (int ii = 1; ii < domainMgr_->tooltxyz_.size(); ii++)
        {
            double *txyzNp1 = &domainMgr_->tooltxyz_[ii][0];
            laserTimeNp1 = txyzNp1[0];
            if (domainMgr_->timet_ <= laserTimeNp1 + small)
            {
                double *txyzN = &domainMgr_->tooltxyz_[ii-1][0];
                laserTimeN = txyzN[0];
                double numer = domainMgr_->timet_ - laserTimeN;
                double denom = laserTimeNp1 - laserTimeN;
                double ratio = numer/denom;
                deltatx = txyzNp1[1] - txyzN[1];
                deltaty = txyzNp1[2] - txyzN[2];
                deltatz = txyzNp1[3] - txyzN[3];

                tx = ratio * deltatx + txyzN[1];
                ty = ratio * deltaty + txyzN[2];
                tz = ratio * deltatz + txyzN[3];
                powerindicator_ = domainMgr_->laserOn_[ii];
                if(init_)
                {
                    oldbeamposx_ = domainMgr_->tooltxyz_[ii-1][1];
                    oldbeamposy_ = domainMgr_->tooltxyz_[ii-1][2];
                    oldbeamposz_ = domainMgr_->tooltxyz_[ii-1][3];
                    oldpowerindicator = domainMgr_->laserOn_[ii-1];
                    init_ = false;
                }//end if

                break;
            }//end if
        }//end for(ii)

        // store current laser location
        beamposx_ = tx;
        beamposy_ = ty;
        beamposz_ = tz;

        // calculate scanning velocity
        double disttraveled = std::sqrt(  
                        (beamposx_ - oldbeamposx_)*(beamposx_ - oldbeamposx_) +
                        (beamposy_ - oldbeamposy_)*(beamposy_ - oldbeamposy_) +
                        (beamposz_ - oldbeamposz_)*(beamposz_ - oldbeamposz_));
        domainMgr_->scanvel_ = disttraveled/meshObj_->delt_;

        // Kevochan: check which track we are on
        if(oldpowerindicator - powerindicator_ > 0 && deltaty > small)
            trackindx_++;

        // Kevochan: check which layer we are on 
        //if(oldpowerindicator - powerindicator_ > 0 && deltatz > small)
        if(powerindicator_ - oldpowerindicator > 0 && meshObj_->powderbed_)
        {
            trackindx_++;
            domainMgr_->nk_ += domainMgr_->ncvzlayer_;
            domainMgr_->nkm1_ = domainMgr_->nk_-1; 
            domainMgr_->nkm2_ = domainMgr_->nk_-2; 
        }//end if
    }//end if
}//end getLaserCoords

//////////////////////////////////////////////////////
//	        getMeltPoolShape                    //
//////////////////////////////////////////////////////
void
BoundCondManager::getMeltPoolShape()
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
    double tsolid = domainMgr_->tsolid_; 
    
    // define some local variables
    int i, j, k;
    double dtdxxinv, dtdzzinv, dtdyyinv, dep, wid;
    double xxmax, xxmin, yymax, yymin;
    
    // tpeak is the maximum temperature in the whole domain 
    tpeak_ = *max_element(&temp_[0], &temp_[0]+(ni*nj*nk)); 

    // determine if a meltpool exist
    if(tpeak_ <= tsolid) 
    {    
        alen_ = 0.0;
        depth_ = 0.0;
        width_ = 0.0;
        return;
    }//end if
    
    //----length-----------------------
    imax_ = istart_;
    imin_ = istart_;
    alen_ = 0.0;

    // get maximum position value in i direction
    for (i=istart_; i<nim1; i++)
    {
        imax_ = i;
        if(temp_[(nkm1*ni*nj)+(jstart_*ni)+i] <= tsolid) 
            break;
    }//end for(i)
    dtdxxinv = (x_[imax_] - x_[imax_+1])/(temp_[(nkm1*ni*nj)+(jstart_*ni)+imax_] - 
                                          temp_[(nkm1*ni*nj)+(jstart_*ni)+imax_+1]);
    xxmax = x_[imax_] + (tsolid - temp_[(nkm1*ni*nj)+(jstart_*ni)+imax_])*dtdxxinv;

    // get minimum position value in i direction
    for (i = istart_; i>0; i--)
    {
        imin_ = i;
        if (temp_[(nkm1*ni*nj)+(jstart_*ni)+i] < tsolid) 
            break;
    }//end for(i)
    dtdxxinv = (x_[imin_] - x_[imin_-1])/(temp_[(nkm1*ni*nj)+(jstart_*ni)+imin_] - 
                                          temp_[(nkm1*ni*nj)+(jstart_*ni)+imin_-1]);
    xxmin = x_[imin_] + (tsolid - temp_[(nkm1*ni*nj)+(jstart_*ni)+imin_])*dtdxxinv;

    // calulate melt pool length in i direction
    alen_ = xxmax - xxmin;

    //----depth-----------------------
    kmin_ = nkm2;
    depth_ = 0.0;
    
    // get minimum position value in k direction
    for (i=1; i<nim1; i++)
    {
        for (k=nkm2; k>0; k--) 
        {
            if (temp_[(k*ni*nj)+(jstart_*ni)+i] < tsolid)
                break;
            kmin_ = min(kmin_, k);
        }//end for(k)
    }//end for(i)
    kmin_--;

    // calulate melt pool depth in k direction
    if (kmin_ == 0)
    {
        depth_ = z_[nkm1] - z_[0];
    }
    else
    {
        for (i=1; i<nim1; i++)
        {
            if (temp_[((kmin_+1)*ni*nj)+(jstart_*ni)+i] < tsolid) 
                continue;
            dtdzzinv = (z_[kmin_] - z_[kmin_-1])/(temp_[(kmin_*ni*nj)+(jstart_*ni)+i] - 
                                                  temp_[((kmin_-1)*ni*nj)+(jstart_*ni)+i]);
            dep = z_[nkm1] - z_[kmin_] + (temp_[(kmin_*ni*nj)+(jstart_*ni)+i] - tsolid)*dtdzzinv;
            depth_ = max(dep, depth_);
        }//end for(i)
    }//end if
    kmax_ = nkm2;
        
    //----width-----------------------
    jmax_ = jstart_;
    jmin_ = jstart_;
    width_ = 0.0;

    // get maximum position value in j direction
    for (i=1; i<nim1; i++)
    {
        for (j=jstart_; j<njm1; j++)
        {
            if (temp_[(nkm1*nj*ni)+(j*ni)+i] < tsolid) 
                break;
            jmax_ = max(jmax_, j);
        }//end for(j)
    }//end for(i)
    jmax_++;

    yymax = y_[0];
    if(jmax_ != jstart_)
    {
        for (i=1; i<nim1; i++)
        {
            if (temp_[(nkm1*ni*nj)+((jmax_-1)*ni)+i] < tsolid) 
                continue;
            dtdyyinv = (y_[jmax_] - y_[jmax_+1])/(temp_[(nkm1*ni*nj)+(jmax_*ni)+i] - 
                                                  temp_[(nkm1*ni*nj)+((jmax_+1)*ni)+i]);
            wid = y_[jmax_] + (tsolid - temp_[(nkm1*ni*nj)+(jmax_*ni)+i])*dtdyyinv;
            yymax = max(wid, yymax);
        }//end for(i)
    }//end if

    // get minimum position value in j direction
    for (i=1; i<nim1; i++)
    {
        for (j=jstart_; j>0; j--)
        {
            if (temp_[(nkm1*ni*nj)+(j*ni)+i] < tsolid)
                break;
            jmin_ = min(jmin_, j);
        }//end for(j)
    }//end for(i)
    jmin_--;

    yymin = y_[jstart_];
    if(jmin_ != jstart_)
    {
        for (i=1; i<nim1; i++)
        {
            if (temp_[(nkm1*ni*nj)+((jmin_+1)*ni)+i] < tsolid)
                continue;
            dtdyyinv = (y_[jmin_]-y_[jmin_-1])/(temp_[(nkm1*ni*nj)+(jmin_*ni)+i] - 
                                                temp_[(nkm1*ni*nj)+((jmin_-1)*ni)+i]);
            wid = y_[jmin_] + (tsolid - temp_[(nkm1*ni*nj)+(jmin_*ni)+i])*dtdyyinv;
            yymin = min(wid,yymin);
        }//end for(i)
    }//end if

    // calulate melt pool width in j direction
    width_ = yymax - yymin;

    //----define solution domain for momentum equations
    istat_ = max(imin_-3, 1);
    iend_ = min(imax_+3, nim2);
    
    jstat_ = max(jmin_-3, 1);
    jend_ = min(jmax_+2, njm2);
    
    kstat_ = max(kmin_-2, 2);
    istatp1_ = istat_ + 1;
    iendm1_ = iend_ - 1;
}//end getMeltPoolShape

//////////////////////////////////////////////////////
//	        applyBCs                            //
//////////////////////////////////////////////////////
void
BoundCondManager::applyBCs(int &ivar)
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
    double a_oxygen_ = domainMgr_->a_oxygen_;

    // define some local variables
    int i, j, k;

    switch(ivar)
    {
        case 1:
        {
            // uvel variables
            double dtdx, fraclu, visu1, term1;

            //----k=nkm1
            for (j=jstat_; j<=jend_; j++)
            {
                for (i=istatp1_; i<=iendm1_; i++)
                {    
                    dtdx = (temp_[(nkm1*nj*ni)+(j*ni)+i] - temp_[(nkm1*nj*ni)+(j*ni)+i-1])*dxpwinv_[i];
                    fraclu = fracl_[(nkm1*nj*ni)+(j*ni)+i]*(1.0 - fracx_[i-1]) + fracl_[(nkm1*nj*ni)+(j*ni)+i-1]*fracx_[i-1];
                    visu1 = vis_[(nkm2*nj*ni)+(j*ni)+i]*vis_[(nkm2*nj*ni)+(j*ni)+i-1]/(vis_[(nkm2*nj*ni)+(j*ni)+i]*(1.0-fracx_[i-1]) + 
                            vis_[(nkm2*nj*ni)+(j*ni)+i-1]*fracx_[i-1]);
                    term1 = fraclu*getSurfaceTensionGradient(temp_[(nkm1*nj*ni)+(j*ni)+i],a_oxygen_)*dtdx/(visu1*dzpbinv_[nkm1]);

                    uvel_[(nkm1*nj*ni)+(j*ni)+i] = uvel_[(nkm2*nj*ni)+(j*ni)+i] + term1;
                }//end for(i)
            }//end for(j)

            // zero out u-velocities in solid
            for (k=kstat_; k<nkm1; k++)
            {
                for (j=jstat_; j<=jend_; j++)
                {
                    uvel_[(k*nj*ni)+(j*ni)+istat_] = 0.0;
                    uvel_[(k*nj*ni)+(j*ni)+iend_] = 0.0;
                }//end for(j)
            }//end for(k)

            break;
        }//end case(1)
        case 2:
        {
            // vvel variables
            double dtdy, fraclv, visv1, term1;
            
            //----k=nkm1 
            for (j=jstat_; j<=jend_; j++)
            {
                for (i=istatp1_; i<=iendm1_; i++)
                {
                    dtdy = (temp_[(nkm1*nj*ni)+(j*ni)+i]-temp_[(nkm1*nj*ni)+((j-1)*ni)+i])*dypsinv_[j];
                    //if(y_[j] >= 0.00035) dtdy=0;                 //avoid unreal temperature gradient at edge of layer	
                    fraclv = fracl_[(nkm1*nj*ni)+(j*ni)+i]*(1.0 - fracy_[j-1]) + fracl_[(nkm1*nj*ni)+((j-1)*ni)+i]*fracy_[j-1];
                    visv1 = vis_[(nkm2*nj*ni)+(j*ni)+i]*vis_[(nkm2*nj*ni)+((j-1)*ni)+i]/(vis_[(nkm2*nj*ni)+(j*ni)+i]*(1.0 - fracy_[j-1]) + 
                            vis_[(nkm2*nj*ni)+((j-1)*ni)+i]*fracy_[j-1]);
                    term1 = fraclv*getSurfaceTensionGradient(temp_[(nkm1*nj*ni)+(j*ni)+i],a_oxygen_)*dtdy/(visv1*dzpbinv_[nkm1]);

                    vvel_[(nkm1*nj*ni)+(j*ni)+i] = vvel_[(nkm2*nj*ni)+(j*ni)+i] + term1;
                }//end for(i)
            }//end for(j)
            
            // zero out v-velocities in solid
            for (k=kstat_; k<nkm1; k++)
            {
                for (j=jstat_; j<=jend_; j++)
                {
                    vvel_[(k*nj*ni)+(j*ni)+istat_] = 0.0;
                    vvel_[(k*nj*ni)+(j*ni)+iend_] = 0.0;
                }//end for(j)
            }//end for(k)

            break;
        }//end case(2)
        case 3:
        {
            // zero out w-velocities in solid
            for (k=kstat_; k<nkm1; k++)
            {
                for (j=jstat_; j<=jend_; j++)
                {
                    wvel_[(k*nj*ni)+(j*ni)+istat_] = 0.0;
                    wvel_[(k*nj*ni)+(j*ni)+iend_] = 0.0;
                }//end for(j)
            }//end for(k)

            break;
        }//end case(3)
        case 4:
        {
            // zero out pressures in solid
            for (k=kstat_; k<nkm1; k++)
            {
                for (j=jstat_; j<=jend_; j++)
                {
                    pp_[(k*nj*ni)+(j*ni)+istat_] = 0.0;
                    pp_[(k*nj*ni)+(j*ni)+iend_] = 0.0;
                }//end for(j)
            }//end for(k)

            break;
        }//end case(4)
        case 5:
        {
            // grab some stuff from the Domain
            double sigm = meshObj_->sigm_;
            double emiss = domainMgr_->emiss_;
            double htckn = domainMgr_->htckn_;
            double htck1 = domainMgr_->htck1_;
            double htci = domainMgr_->htci_;
            double htcj = domainMgr_->htcj_;

            // grab some stuff from the mesh
            double tempamb = meshObj_->tempamb_;

            // enthalpy variables
            double hlossradia, hlossconvec, ctmp1; 
            double mdot, psat, evaplatent, evaptemp, molarmass, evapfact, 
                   pressamb, gasconst;
            double hlossevap = 0.0;
            
            // if evporation is activated, grab necessary values
            if(meshObj_->isevap_)
            {
                evaplatent = domainMgr_->evaplatent_; 
                evaptemp = domainMgr_->evaptemp_; 
                molarmass = domainMgr_->molarmass_; 
                evapfact = domainMgr_->evapfact_;  
                pressamb = meshObj_->pressamb_;
                gasconst = meshObj_->gasconst_;
            }//end if(isevap)
             
            //----top (convection, radiation and maybe evporation)
            ahtoploss_ = 0.0;
            for (j=1; j<njm1; j++)
            {
                for (i=1; i<nim1; i++)
                {
                    // include heat loss due to evpaoration BC (if desired)
                    if(meshObj_->isevap_)
                    {
                        psat = pressamb*exp(-evaplatent*molarmass/gasconst*
                                (1.0/temp_[(nkm1*nj*ni)+(j*ni)+i] - 1.0/evaptemp));
                        mdot = evapfact*psat*sqrt(molarmass/(2.0*M_PI*gasconst*temp_[(nkm1*nj*ni)+(j*ni)+i]));
                        hlossevap = mdot*evaplatent;
                    }//end if(isevap)
                    
                    // calculate heat loss due to radiation and heat convection
                    hlossradia = emiss*sigm*(fastPow(temp_[(nkm1*nj*ni)+(j*ni)+i], 4) - fastPow(tempamb, 4));
                    hlossconvec = htckn*(temp_[(nkm1*nj*ni)+(j*ni)+i] - tempamb);

                    ctmp1 = diff_[(nkm2*nj*ni)+(j*ni)+i]*dzpbinv_[nkm1];
                    enthalpy_[(nkm1*nj*ni)+(j*ni)+i] = enthalpy_[(nkm2*nj*ni)+(j*ni)+i] + (heatin_[(j*ni)+i] - 
                                       hlossradia-hlossradia-hlossevap)/ctmp1;
                    
                    // heat loss from the top surface
                    ahtoploss_ += (hlossradia + hlossradia + hlossevap)*areaij_[(j*ni)+i]; 
                }//end for(i)
            }//end for(j)

            //----east and west (mainly convection) 
            for (k=1; k<nkm1; k++)
            {
                for (j=1; j<njm1; j++)
                {
                    hlossconvec = htci*(temp_[(k*nj*ni)+(j*ni)+0] - tempamb);
                    ctmp1 = diff_[(k*nj*ni)+(j*ni)+1]*dxpwinv_[1];
                    enthalpy_[(k*nj*ni)+(j*ni)+0] = enthalpy_[(k*nj*ni)+(j*ni)+1] - hlossconvec/ctmp1;

                    hlossconvec = htci*(temp_[(k*nj*ni)+(j*ni)+nim1] - tempamb);
                    ctmp1 = diff_[(k*nj*ni)+(j*ni)+nim2]*dxpwinv_[nim2];
                    enthalpy_[(k*nj*ni)+(j*ni)+nim1] = enthalpy_[(k*nj*ni)+(j*ni)+nim2] - hlossconvec/ctmp1;
                }//end for(j)
            }//end for(k)
            
            //----north and south (mainly convection)  
            for (k=1; k<nkm1; k++)
            {
                for (i=1; i<nim1; i++)
                {
                    hlossconvec = htcj*(temp_[(k*nj*ni)+(0*ni)+i] - tempamb);
                    ctmp1 = diff_[(k*nj*ni)+(1*ni)+i]*dypsinv_[1];
                    enthalpy_[(k*nj*ni)+(0*ni)+i] = enthalpy_[(k*nj*ni)+(1*ni)+i] - hlossconvec/ctmp1;

                    hlossconvec = htcj*(temp_[(k*nj*ni)+(njm1*ni)+i] - tempamb); 
                    ctmp1 = diff_[(k*nj*ni)+(njm2*ni)+i]*dypsinv_[njm2];
                    enthalpy_[(k*nj*ni)+(njm1*ni)+i] = enthalpy_[(k*nj*ni)+(njm2*ni)+i] - hlossconvec/ctmp1;
                }//end for(i)
            }//end for(k)

            break;
        }//end case(5)
        case 6:
        {
            double initval = domainMgr_->initconcentration_;
            double diffspeciesl = domainMgr_->dens_*domainMgr_->massdiff_;  
            //----top/bottom (k=nkm1/k=0)
            for (j=1; j<njm1; j++)
            {
                for (i=1; i<nim1; i++)
                {
                    if(temp_[(nkm1*nj*ni)+(j*ni)+i] > domainMgr_->tliquid_)
            	    {
                        avgconcentration_[(nkm1*nj*ni)+(j*ni)+i] = 2.0*(13.0-avgconcentration_[(nkm2*nj*ni)+(j*ni)+i])/dzpbinv_[nkm1]/diffspeciesl+avgconcentration_[(nkm2*nj*ni)+(j*ni)+i];
                        //avgconcentration_[(nkm1*nj*ni)+(j*ni)+i] = 1;
            	    }
                    else
            	    {
                        avgconcentration_[(nkm1*nj*ni)+(j*ni)+i] = avgconcentration_[(nkm2*nj*ni)+(j*ni)+i];
            	    }//end if

                   avgconcentration_[(0*nj*ni)+(j*ni)+i] = initval; 
                }//end for(i)
            }//end for(j)
            
            //----north/south (j=njm1/j=0)
            for (k=0; k<nk; k++)
            {
                for (i=0; i<ni; i++)
                {
                    avgconcentration_[(k*nj*ni)+(0*ni)+i] = avgconcentration_[(k*nj*ni)+(1*ni)+i];
                    avgconcentration_[(k*nj*ni)+(njm1*ni)+i] = avgconcentration_[(k*nj*ni)+((njm1-1)*ni)+i];
                }//end for(i)
            }//end for(k)
            
            //----east/west (i=nim1/i=0)
            for (k=0; k<nk; k++)
            {
                for (j=0; j<nj; j++)
                {
                    avgconcentration_[(k*nj*ni)+(j*ni)+0] = avgconcentration_[(k*nj*ni)+(j*ni)+1];
                    avgconcentration_[(k*nj*ni)+(j*ni)+nim1] = avgconcentration_[(k*nj*ni)+(j*ni)+nim1-1];
                }//end for(i)
            }//end for(k)

            break;
        }//end case(6)
    }//end switch
}//end applyBCs

//////////////////////////////////////////////////////
//	        getSurfaceTensionGradient           //
//////////////////////////////////////////////////////
double 
BoundCondManager::getSurfaceTensionGradient(double &t, 
                                            double &a0)
{
    // define some local variables
    double gamma0, Ts, k1, H0, K;
    double surface_tension_gradient;
    double A = 4.3e-4;
    double R = 8314;
    
    //----Fe-O system
    gamma0 = 1.943;
    Ts = 2.03e-8;
    k1 = 1.38e-2;
    H0 = -1.46e8;

    //----Fe-S system
    //gamma0 = 1.943;
    //Ts = 1.3e-8;
    //k1 = 3.18e-3;
    //H0 = -1.66e8;

    // calulate surfactant-affected surface tension gradient
    K = k1*exp(-H0/R/t);
    surface_tension_gradient = (-A-R*Ts*log(1+K*a0) - 
                                 (K*a0/(1 + K*a0))*(Ts*H0/t))*1;

    if(t >= 2500) 
        surface_tension_gradient = 0;

    // if there are no surfactants use user-given value
    if(a0 < 0) 
        surface_tension_gradient = domainMgr_->dgdtp_;

    return surface_tension_gradient;
}//end getSurfaceTensionGradient

//////////////////////////////////////////////////////
//	        calculateFluxes                     //
//////////////////////////////////////////////////////
void
BoundCondManager::calculateFluxes()
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
    double hlatnt = domainMgr_->hlatnt_;

    // grab some stuff from the mesh
    double delt = meshObj_->delt_;
      
    // initialize all fluxes to zero
    fluxwest_ = 0.0;
    fluxeast_ = 0.0;
    fluxbottom_ = 0.0;
    fluxtop_ = 0.0;
    fluxnorth_ = 0.0;
    fluxsouth_ = 0.0;
    accul_ = 0.0;
    heatvol_ = 0.0;
    
    // define some local variables
    int i, j, k; 
    double fluxi1, fluxl1, fluxk1, fluxn1, fluxj1, fluxm1, dh1;

    // calcaulate the fluxes in the system
    #pragma omp parallel private(i, j, k, fluxi1, fluxl1, fluxk1, fluxn1, fluxj1, fluxm1, dh1)
    {
        //----i=1 & i=nim1
        #pragma omp for reduction(+:fluxwest_,fluxeast_) collapse(2) nowait 
        for (k=1; k<nkm1; k++)
        {
            for (j=1; j<njm1; j++)
            {
                fluxi1 = diff_[(k*nj*ni)+(j*ni)+1]*(enthalpy_[(k*nj*ni)+(j*ni)+0] - 
                           enthalpy_[(k*nj*ni)+(j*ni)+1])*dxpwinv_[1];
                fluxl1 = diff_[(k*nj*ni)+(j*ni)+nim2]*(enthalpy_[(k*nj*ni)+(j*ni)+nim1] - 
                           enthalpy_[(k*nj*ni)+(j*ni)+nim2])*dxpwinv_[nim1];
                fluxwest_ += areajk_[(k*nj)+j]*fluxi1;
                fluxeast_ += areajk_[(k*nj)+j]*fluxl1;
            }//end for(j)
        }//end for(k)

        //----k=nkm1 and k=1
        #pragma omp for reduction(+:fluxbottom_,fluxtop_) collapse(2) nowait
        for (j=1; j<njm1; j++)
        {
            for (i=1; i<nim1; i++)
            {
                fluxk1 = diff_[(1*nj*ni)+(j*ni)+i]*(enthalpy_[(0*nj*ni)+(j*ni)+i] - 
                           enthalpy_[(1*nj*ni)+(j*ni)+i])*dzpbinv_[1];
                fluxn1 = diff_[(nkm2*nj*ni)+(j*ni)+i]*(enthalpy_[(nkm1*nj*ni)+(j*ni)+i] - 
                           enthalpy_[(nkm2*nj*ni)+(j*ni)+i])*dzpbinv_[nkm1];
                fluxbottom_ += areaij_[(j*ni)+i]*fluxk1;
                fluxtop_ += areaij_[(j*ni)+i]*fluxn1;
            }//end for(i)
        }//end for(j)

        //----j=1 and j=njm1
        #pragma omp for reduction(+:fluxsouth_,fluxnorth_) collapse(2) nowait
        for (k=1; k<nkm1; k++)
        {
            for (i=1; i<nim1; i++)
            {
                fluxj1 = diff_[(k*nj*ni)+(1*ni)+i]*(enthalpy_[(k*nj*ni)+(0*ni)+i] - 
                           enthalpy_[(k*nj*ni)+(1*ni)+i])*dypsinv_[1];         
                fluxm1 = diff_[(k*nj*ni)+(njm2*ni)+i]*(enthalpy_[(k*nj*ni)+(njm1*ni)+i] - 
                           enthalpy_[(k*nj*ni)+(njm2*ni)+i])*dypsinv_[njm1];
                fluxsouth_ += areaik_[(k*ni)+i]*fluxj1;
                fluxnorth_ += areaik_[(k*ni)+i]*fluxm1;
            }//end for(i)
        }//end for(k)

        // determine heat accumulation
        if(!domainMgr_->steady_)
        {
            #pragma omp for reduction(+:accul_,heatvol_) collapse(3) nowait 
            for (k=1; k<nkm1; k++)
            {
                for (j=1; j<njm1; j++)
                {
                    for (i=1; i<nim1; i++)
                    {
                        dh1 = enthalpy_[(k*nj*ni)+(j*ni)+i] - hnot_[(k*nj*ni)+(j*ni)+i] + (fracl_[(k*nj*ni)+(j*ni)+i] - 
                                fraclnot_[(k*nj*ni)+(j*ni)+i])*hlatnt;
                        accul_ += volume_[(k*nj*ni)+(j*ni)+i]*rho_[(k*nj*ni)+(j*ni)+i]*dh1/delt;
                        heatvol_ += sourceinput_[(k*nj*ni)+(j*ni)+i]*volume_[(k*nj*ni)+(j*ni)+i];
                    }//end for(i)
                }//end for(j)
            }//end for(k)
        }//end if    
    }//end omp parallel

    // calculate total heat loss
    heatout_ = (fluxnorth_ + fluxbottom_ + fluxwest_ + fluxeast_ + fluxsouth_) - ahtoploss_;

    // check energy balence (for convergence)
    if(domainMgr_->steady_)
    {
        ratio_ = -heatout_/heatinlaser_;
    }
    else
    {
        ratio_ = (heatvol_ + heatinlaser_)/(accul_ - heatout_);
    }//end if
}//end calculateFluxes

//////////////////////////////////////////////////////
//	        fastPow         	            //
//////////////////////////////////////////////////////
inline double 
BoundCondManager::fastPow(double base, int exponent) 
{
    // faster alternative to the pow() function 
    double out = base;
    for(int i=1; i<exponent; i++)
        out *= base;

   return out;
}//end fastPow
