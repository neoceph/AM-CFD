// standard headers
#include <vector>
#include <chrono>
#include <ctime>
#include <stdio.h>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <iomanip>
#include <omp.h>
#include "time.h"

// user-defined headers
#include <PrintManager.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
PrintManager::PrintManager(Mesh *meshObj,
                           DomainManager *domainMgr,
                           SolutionVariableManager *solVarMgr,
                           BoundCondManager *bcMgr,
                           CFDSolverManager *cfdSolvMgr,
                           FreeSurfaceManager *freesurMgr,
                           string &FileOut,
                           string &directoryName,
                           string &baseOutName)
                           : meshObj_(meshObj),
                             domainMgr_(domainMgr),
                             solVarMgr_(solVarMgr),
                             bcMgr_(bcMgr),
                             cfdSolvMgr_(cfdSolvMgr),
                             freesurMgr_(freesurMgr),
                             FileOut_(FileOut),
                             directoryName_(directoryName),
                             baseOutName_(baseOutName)
{
}

//////////////////////////////////////////////////////
//		initialization             	    //
//////////////////////////////////////////////////////
void 
PrintManager::initialization()
{
    outputInitHeader();
    outputAssignInputTimer();
    outputGridTime();
    outputTimestepStats();
    outputSolverMode();
    if(meshObj_->outputsection_)
        outputSectionInfo();
    outputInitFooter();

    cfdSolvMgr_->initializeSolver();
    bcMgr_->initializeBoundaries();
    freesurMgr_->initializeFreeSurface();
    initializePrinter();

    // Kevochan: For powder-bed simulations only activate first layer
    if(meshObj_->powderbed_)
    {
        domainMgr_->nk_ -= (domainMgr_->numlayers_ - 1)*domainMgr_->ncvzlayer_; 
        domainMgr_->nkm1_ = domainMgr_->nk_-1; 
        domainMgr_->nkm2_ = domainMgr_->nk_-2; 
        baseheightk_ = domainMgr_->nkm2_ - domainMgr_->ncvzlayer_;
    }//end if
}//end initialization

//////////////////////////////////////////////////////
//		outputInitHeader              	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputInitHeader()
{
    auto begin  = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    const int setSpace = 50;

    cout << "===================================================================\n";
    cout << "                   ** BEGAN SETTING UP DOMAIN **\n\n";
    cout << SecondEleWidth_ << "Timers for Creating Domain:\n";
}//end outputInitHeader

//////////////////////////////////////////////////////
//		outputAssignInputTimer        	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputAssignInputTimer()
{
    auto begin  = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    const int setSpace = 50;

    begin = std::chrono::high_resolution_clock::now();
    domainMgr_->assignInputValues();
    end = std::chrono::high_resolution_clock::now();

    cout << std::setw(setSpace) << std::right
         << "Time for assigning input vaules: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << std::endl;

    if(meshObj_->udtoolpath_)
    {
        begin = std::chrono::high_resolution_clock::now();
        domainMgr_->getToolpath();
        end = std::chrono::high_resolution_clock::now();

        cout << std::setw(setSpace) << std::right
             << "Time for reading in toolpath: "
             << std::setw(16) << std::setprecision(8) << std::right
             <<  1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << std::endl;
    }//end if
}//end outputAssignInputTimer

//////////////////////////////////////////////////////
//		outputGridTime              	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputGridTime()
{
    auto begin  = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    const int setSpace = 50;

    begin = std::chrono::high_resolution_clock::now();
    cfdSolvMgr_->generateGrid();
    end = std::chrono::high_resolution_clock::now();

    cout << std::setw(setSpace) << std::right
         << "Time for generating user-defined grid: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "\n";
    cout << "\n" << endl;
    cfdSolvMgr_->init_ = false;
}//end outputGridTime

//////////////////////////////////////////////////////
//		outputTimeStepStats            	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputTimestepStats()
{
    // grab some stuff from the mesh
    double outtime = meshObj_->outtime_;
    double finaltime = meshObj_->finaltime_;
    double delt = meshObj_->delt_;
     
    // define spacing
    const int setSpace = 50;

    // calculate outputs/approximate # of timesteps
    int numsteps = (int)ceil(finaltime/delt);
    int numout = (int)ceil(finaltime/outtime);
    int numoutsteps = (int)ceil(outtime/delt);

    // calulate number nodes with out boundaries
    int numboundnodes = 2;
    int sumcvx = domainMgr_->ni_-numboundnodes;
    int sumcvy = domainMgr_->nj_-numboundnodes;
    int sumcvz = domainMgr_->nk_-numboundnodes;

    cout << SecondEleWidth_ << "Mesh and Timestep Statistics:\n";
    cout << std::setw(setSpace) << std::right
         << "Total number of internal control volumes: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  sumcvx*sumcvy*sumcvz << "\n";
    cout << std::setw(setSpace) << std::right
         << "Given time step size: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  delt << "\n";
    cout << std::setw(setSpace) << std::right
         << "Approximate # of timesteps: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  numsteps << "\n";
    cout << std::setw(setSpace) << std::right
         << "Approximate # of outputs: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  numout << "\n";
    cout << std::setw(setSpace) << std::right
         << "Approximate # of steps between outputs: "
         << std::setw(16) << std::setprecision(8) << std::right
         <<  numoutsteps << "\n";
    cout << "\n" << endl;
}//end outputTimeStepStats

//////////////////////////////////////////////////////
//		outputSolverMode            	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputSolverMode()
{
    string base = "Solving the Following Equations (Discretized by ";
    string detail;
    if(meshObj_->upwind_)
    {
        detail = "Upwinding";
    }
    else if(meshObj_->powerlaw_)
    {
        detail = "Power-Law";
    }
    else if(meshObj_->exponential_)
    {
        detail = "Exponential-Law";
    }//end if

    string output = base + detail + "):";
    cout << SecondEleWidth_ << output <<"\n";

    if(meshObj_->navierstokes_)
    {
        cout << "\tConvervation of Mass" << "\n";
        cout << "\tConvervation of Momentum (3D)" << "\n";
    }//end if
    cout << "\tConvervation of Energy" << "\n";
    if(meshObj_->species_)
    {
        cout << std::right << "\tConvervation of Species (multi-component)" << "\n";
    }//end if
    if(meshObj_->energyfreesurface_)
    {
        cout << std::right << "\tBalence of Energy at the Free Surface" << "\n";
    }//end if
    cout << std::endl;
}//outputSolverMode

//////////////////////////////////////////////////////
//		outputSectionInfo            	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputSectionInfo()
{
    cout << "\n" << SecondEleWidth_ 
         << "Outputting an Extracted Section Defined by:" << "\n";
    cout << setprecision(5) << scientific 
         << "\t(xmin: " << domainMgr_->gxmin_ 
         << ", ymin: " << domainMgr_->gymin_
         << ", zmin: " << domainMgr_->gzmin_ << ")\n";
    cout << setprecision(5) << scientific 
         << "\t(xmax: " << domainMgr_->gxmax_ 
         << ", ymax: " << domainMgr_->gymax_
         << ", zmax: " << domainMgr_->gzmax_ << ")\n";
    cout << std::endl;
}//outputSectionInfo

//////////////////////////////////////////////////////
//		outputInitFooter              	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputInitFooter()
{
    cout << "                   ** ENDED SETTING UP DOMAIN **\n";
    cout << "===================================================================" 
         << endl;
}//end outputInitFooter

//////////////////////////////////////////////////////
//		initializePrinter		    //
//////////////////////////////////////////////////////
void 
PrintManager::initializePrinter()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    gridz_ = domainMgr_->nk_;

    // initialize local arrays
    fill(&auvel_[0], &auvel_[0]+nk*nj*ni, 0.0);
    fill(&avvel_[0], &avvel_[0]+nk*nj*ni, 0.0);
    fill(&awvel_[0], &awvel_[0]+nk*nj*ni, 0.0);

    // map interpolated velocites to viewing window
    solVarMgr_->mapSolutionVars_["auvel"] = &auvel_[0];     
    solVarMgr_->mapSolutionVars_["avvel"] = &avvel_[0];     
    solVarMgr_->mapSolutionVars_["awvel"] = &awvel_[0];     
    
    // set up appropipate pointers to cfd solver manager
    temp_ = solVarMgr_->mapSolutionVars_["temp"]; 
    uvel_ = solVarMgr_->mapSolutionVars_["uvel"]; 
    vvel_ = solVarMgr_->mapSolutionVars_["vvel"]; 
    wvel_ = solVarMgr_->mapSolutionVars_["wvel"]; 
    if(meshObj_->solidificationparams_)
    {
        rate_ = solVarMgr_->mapSolutionVars_["rate"]; 
        grad_ = solVarMgr_->mapSolutionVars_["grad"]; 
    }//end if

    x_ = solVarMgr_->mapGeometryVars_["x"];       
    y_ = solVarMgr_->mapGeometryVars_["y"];       
    z_ = solVarMgr_->mapGeometryVars_["z"];       

    // set up appropipate pointers to free surface manager (if necessary)
    if(meshObj_->energyfreesurface_)
    {
        tempr_ = solVarMgr_->mapSolutionVars_["tempr"]; 
        zr_ = solVarMgr_->mapSolutionVars_["zr"]; 
    }//end if
}//end initializePrinter

//////////////////////////////////////////////////////
//		scanOutputData			    //
//////////////////////////////////////////////////////
void 
PrintManager::scanOutputData()
{
    // Check for outputs
    vector<string> leftOverData;
    for (int i = 0; i < meshObj_->outputScalarNamesNODE_.size(); ++i)
    {
        string outName = meshObj_->outputScalarNamesNODE_[i];
        std::map<string, double*>::iterator found = 
            solVarMgr_->mapSolutionVars_.find(outName);
        if (found != solVarMgr_->mapSolutionVars_.end())
        {
            // make sure user is outputting consistent velocities
            // (we currently use a staggered grid in this model)
            if(outName == "uvel")
            {
                outName = "auvel";
            }
            else if(outName == "vvel")
            {
                outName = "avvel";
            }
            else if(outName == "wvel")
            {
                outName = "awvel";
            }//end if
            outScalarNames_.push_back(outName);
        }
        else
        {
            leftOverData.push_back(outName);
        }//end if
    }//end for(i)
  
    for (int i = 0; i < meshObj_->outputVectorNamesNODE_.size(); ++i)
    {
        string outName = meshObj_->outputVectorNamesNODE_[i];
        std::map<string, double*>::iterator found = 
        solVarMgr_->mapSolutionVars_.find(outName);
        if (found != solVarMgr_->mapSolutionVars_.end())
        {
            outVectorNames_.push_back(outName);
        }
        else
        {
            leftOverData.push_back(outName);
        }//end if
    }//end for(i)
  
    if(leftOverData.size() != 0)
    {
        cout << "===============================================================\n";
        cout << "**The following data outputs are not supported:   **\n";
        for (int i = 0; i < leftOverData.size(); ++i)
        {
            cout << "  " << leftOverData[i] << "\n";
        }//end for(i)
        cout << "===============================================================\n\n";
        cout << endl;
    }//end if

    leftOverData.shrink_to_fit();
}//end scanOutputData

//////////////////////////////////////////////////////
//		outputResiduals			    //
//////////////////////////////////////////////////////
void 
PrintManager::outputResiduals()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int nkm2 = domainMgr_->nkm2_;
    double tsolid = domainMgr_->tsolid_;
    double tliquid = domainMgr_->tliquid_;
    double scanvel = domainMgr_->scanvel_;

    // grab some stuff from the bc 
    int istatp1 = bcMgr_->istatp1_;
    int iendm1 = bcMgr_->iendm1_;
    int jstat = bcMgr_->jstat_;
    int jend = bcMgr_->jend_;
    double tpeak = bcMgr_->tpeak_;
    int jstart = bcMgr_->jstart_;

    // define local variables
    int i, j, k;
    double xso, tso, xsob, tsob, xliq, xliqm1, tliq, tliqm1, crs, crl;

    umax_ = 0.0;
    vmax_ = 0.0;
    wmax_ = 0.0;
    if(tpeak > tsolid)
    {
        for (j=jstat; j<=jend; ++j)
        {
           for (i=istatp1; i<=iendm1; ++i)
           {
               umax_ = max(umax_, abs(uvel_[(nkm2*nj*ni)+(j*ni)+i]));
               vmax_ = max(vmax_, abs(vvel_[(nkm2*nj*ni)+(j*ni)+i]));
               wmax_ = max(wmax_, abs(wvel_[(nkm2*nj*ni)+(j*ni)+i]));
           }//end for(i)
        }//end for(j)
    }//end if

    for (i=1; i<nim1; ++i)
    {
        if(temp_[(nkm2*nj*ni)+(jstart*ni)+i] > tsolid-290)
        {
            tsob = temp_[(nkm2*nj*ni)+(jstart*ni)+i-1];
            xsob = x_[i-1];
            break;
        }//end if
    }//end for(i)

    for (i=1; i<nim1; ++i)
    {
        if(temp_[(nkm2*nj*ni)+(jstart*ni)+i] > tsolid)
        {
            tso = temp_[(nkm2*nj*ni)+(jstart*ni)+i-1];
            xso = x_[i-1];
            break;
        }//end if
    }//end for(i)

    for (i=1; i<nim1; ++i)
    {
        if(temp_[(nkm2*nj*ni)+(jstart*ni)+i] > tliquid)
        {
            tliq = temp_[(nkm2*nj*ni)+(jstart*ni)+i];
            tliqm1 = temp_[(nkm2*nj*ni)+(jstart*ni)+i-1];
    
            xliq = x_[i];
            xliqm1 = x_[i-1];
            break;
        }//end if
    }//end for(i)
    
    crs = (tso-tsob)/(xso-xsob)*scanvel;
    //crl = (tliq-tliqm1)/(xliq-xliqm1)*scanvel;
    crl = (tliq-tso)/(xliq-xso)*scanvel;
    
    // grab some stuff from the cfd solver
    int niter = cfdSolvMgr_->niter_;
    double resoru = cfdSolvMgr_->resoru_;
    double resorv = cfdSolvMgr_->resorv_;
    double resorw = cfdSolvMgr_->resorw_;
    double resorh = cfdSolvMgr_->resorh_;
    double resorm = cfdSolvMgr_->resorm_;
    double resorc = cfdSolvMgr_->resorc_;

    // grab  stuff from the bc
    double flux_west = bcMgr_->fluxwest_; 
    double flux_east = bcMgr_->fluxeast_; 
    double flux_top = bcMgr_->fluxtop_;
    double flux_bottom = bcMgr_->fluxbottom_; 
    double flux_north = bcMgr_->fluxnorth_; 
    double flux_south = bcMgr_->fluxsouth_;
    double ratio = bcMgr_->ratio_;
    double accul = bcMgr_->accul_;
    double heatout = bcMgr_->heatout_;
    double heatinlaser = bcMgr_->heatinlaser_;
    double ahtoploss = bcMgr_->ahtoploss_;
    double alen = bcMgr_->alen_;
    double depth = bcMgr_->depth_;
    double width = bcMgr_->width_;

    // open the residual log file
    residFile_.open(baseOutName_ + "_Residuals.txt", ios::out | ios::app);
    
    // output the relevant solver data 
    if(domainMgr_->steady_)
    {
        residFile_ << "====================================================================================================\n";
        residFile_ << "  iter    res_enth   res_mass   res_u      res_v      res_w      res_concen\n";
        residFile_ << setprecision(2) << fixed << setw(8) << niter << scientific << setw(15) << resorh 
                                      << setw(11) << resorm << setw(11) << resoru << setw(11) << resorv 
                                      << setw(11) << resorw << setw(11) << resorc << " \n";

        residFile_ << "  Tmax       umax       vmax       wmax        length      depth       width\n";
        residFile_ << setprecision(2) << fixed << setw(10) << tpeak << scientific << setw(11) << umax_ << setw(11) << vmax_ 
                                      << setw(11) << wmax_ << setw(12) << alen << setw(12) << depth << setw(12) << width << "\n";

        residFile_ << "  north    south     top   toploss   bottom    west      east    hout    accu    hin      ratio\n";
        residFile_ << setprecision(1) << fixed << setw(6) << flux_north << setw(9) << flux_south << setw(9) << flux_top 
                                      << setw(9) << ahtoploss << setw(9) << flux_bottom << setw(9) << flux_west 
                                      << setw(9) << flux_east << setw(9) << heatout << setw(8) << accul << setw(8) << heatinlaser 
                                      << setw(8) << ratio << "\n";
        residFile_ << "====================================================================================================\n";
    }
    else
    {
        residFile_ << "====================================================================================================\n";
        residFile_ << "  time        iter  tot_iter    res_enth   res_mass   res_u      res_v      res_w      res_concen\n";
        residFile_ << setprecision(2) << scientific << setw(10) << domainMgr_->timet_ << fixed << setw(6) << niter 
                                      << setw(9) << cfdSolvMgr_->itertot_ << scientific << setw(15) << resorh 
                                      << setw(11) << resorm << setw(11) << resoru << setw(11) << resorv << setw(11) << resorw 
                                      << setw(11) << resorc << " \n";

        residFile_ << "  Tmax       umax       vmax       wmax        length      depth       width\n";
        residFile_ << setprecision(2) << fixed << setw(10) << tpeak << scientific << setw(11) << umax_ << setw(11) << vmax_ 
                                      << setw(11) << wmax_ << setw(12) << alen << setw(12) << depth << setw(12) << width << "\n";

        residFile_ << "  north    south     top   toploss   bottom    west      east    hout    accu    hin      ratio\n";
        residFile_ << setprecision(1) << fixed << setw(6) << flux_north << setw(9) << flux_south << setw(9) << flux_top 
                                      << setw(9) << ahtoploss << setw(9) << flux_bottom << setw(9) << flux_west 
                                      << setw(9) << flux_east << setw(9) << heatout << setw(8) << accul << setw(8) << heatinlaser 
                                      << setw(8) << ratio << "\n\n";

        residFile_ << setprecision(4) << scientific 
                   << "  Cooling rate from solidus to 290K below solidus:     " << setw(14) << crs << "\n"
                   << "  Cooling rate from liquidus to solidus:               " << setw(14) << crl << "\n";
        residFile_ << "====================================================================================================\n";
    }//end if
    
    // close the residual log file
    residFile_.close();
}//end outputResiduals

//////////////////////////////////////////////////////
//		outputSection			    //
//////////////////////////////////////////////////////
void 
PrintManager::outputSection()
{
    // exit the subroutine if necessary
    tracksect_ += meshObj_->delt_;
    if(tracksect_ < domainMgr_->toutsect_ && !domainMgr_->finish_)
        return;

    // begin timer
    auto begin = std::chrono::high_resolution_clock::now();

    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nkm1 = domainMgr_->nkm1_;
    double timet = domainMgr_->timet_;
    double gxmin = domainMgr_->gxmin_;  
    double gxmax = domainMgr_->gxmax_;
    double gymin = domainMgr_->gymin_;
    double gymax = domainMgr_->gymax_;
    double gzmin = domainMgr_->gzmin_;
    double gzmax = domainMgr_->gzmax_;

    // count number of nodes in each direction
    int outx = 0;
    int outy = 0;
    int outz = 0;

    for (int i=0; i<ni; i++)
    {
        if(x_[i] >= gxmin && x_[i] <= gxmax)
            outx++;
    }//end for(i)

    for (int j=0; j<nj; j++)
    {
        if(y_[j] >= gymin && y_[j] <= gymax)
            outy++;
    }//end for(j)

    for (int k=0; k<gridz_-1; k++)
    {
        if(z_[k] >= gzmin && z_[k] <= gzmax)
            outz++;
    }//end for(k)

    // determine the current layer of the simulation
    double tol = 1e-9;
    int layer = (int)((z_[nkm1] - z_[baseheightk_] + tol)/domainMgr_->layerheight_);
    
    // open the section output file
    sectionFile_.open(baseOutName_ + "_ExtractedSection" + to_string(layer) + ".txt", ios::out | ios::app);

    // Kevochan: output according to Lichao's file format (Feb., 2020)
    sectionFile_ << "t " <<  "ix " "jy " << "kz " << "\n"
                 << timet << " "<< outx << " " << outy << " " << outz << "\n";
    sectionFile_ << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << setw(15) << "tn" << "\n";
    for (int k=0; k<gridz_-1; ++k)
    {
        for (int j=0; j<nj; ++j)
        {
            for (int i=0; i<ni; ++i)
            {
                if(x_[i] >= gxmin && x_[i] <= gxmax && y_[j] >= gymin && y_[j] <= gymax && z_[k] >= gzmin && z_[k] <= gzmax && k<nkm1)
                {
                    sectionFile_ << setprecision(6) << scientific << setw(15) << x_[i] << setw(15) << y_[j] 
                                 << setw(15) << z_[k] << setw(15) << temp_[(k*nj*ni)+(j*ni)+i] << "\n";
                }
                else if(x_[i] >= gxmin && x_[i] <= gxmax && y_[j] >= gymin && y_[j] <= gymax && z_[k] >= gzmin && z_[k] <= gzmax && k>=nkm1)
                {
                    sectionFile_ << setprecision(6) << scientific << setw(15) << x_[i] << setw(15) << y_[j] 
                                 << setw(15) << z_[k] << setw(15) << 0.0 << "\n";
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)

    // close the section output file
    sectionFile_ << "\n";
    sectionFile_.close();
    tracksect_ = 0.0;

    // end timer and update total time spent outputting section data
    auto end = std::chrono::high_resolution_clock::now();
    secttime_ += 1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
    dumpsect_ += 1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count(); 
}//end outputSection

//////////////////////////////////////////////////////
//		outputFinalSteady       	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputFinalSteady()
{
    // open the residual log file
    residFile_.open(baseOutName_ + "_residuals.txt", ios::out | ios::app);
    residFile_ << setprecision(4) << fixed;

    // append final values
    residFile_ << "\n Some important calculated parameters at the end of heating cycle\n";
    residFile_ << "===================================================================\n";
    residFile_ << "  Length of the pool     (cm):        " << setw(14) << bcMgr_->alen_        << "\n";
    residFile_ << "  Depth of the pool      (cm):        " << setw(14) << bcMgr_->depth_       << "\n";
    residFile_ << "  Half-width of the pool (cm):        " << setw(14) << bcMgr_->width_       << "\n";
    residFile_ << "  Peak temperature       (K):         " << setw(14) << bcMgr_->tpeak_       << "\n";
    residFile_ << "  Maximum u-velocity     (cm/s):      " << setw(14) << umax_                << "\n";
    residFile_ << "  Maximum v-velocity     (cm/s):      " << setw(14) << vmax_                << "\n";
    residFile_ << "  Maximum w-velocity     (cm/s):      " << setw(14) << wmax_                << "\n";
    residFile_ << "  Rate of heat input     (cal/s):     " << setw(14) << bcMgr_->heatinlaser_ << "\n";
    residFile_ << "  Rate of heat output    (cal/s):     " << setw(14) << bcMgr_->heatout_     << "\n";
    residFile_ << "  Ratio of heat input to heat output: " << setw(14) << bcMgr_->ratio_       << "\n";
    residFile_ << "===================================================================\n";

    // close the residual log file
    residFile_.close();
}//end outputFinalSteady

//////////////////////////////////////////////////////
//		outputFinalGeometry       	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputFinalGeometry()
{

    // open the residual log file
    finalgeomFile_.open(meshObj_->geomOutName_, ios::out | ios::app);

    // output the final melt pool geometry 
    finalgeomFile_ << "\n This is the final melt pool geometry for the given process parameters\n";
    finalgeomFile_ << "===================================================================\n";
    finalgeomFile_ << "  Laser Power    Laser Speed      length      depth      width\n";
    finalgeomFile_ << setprecision(2) << fixed << setw(8) << domainMgr_->alaspow_ << setw(13) << domainMgr_->scanvel_ 
                   << scientific << setw(21) << bcMgr_->alen_ 
                                 << setw(12) << bcMgr_->depth_ 
                                 << setw(11) << bcMgr_->width_ << "\n";
    finalgeomFile_ << "===================================================================\n";

    // close the residual log file
    finalgeomFile_.close();
}//end outputFinalGeometry

//////////////////////////////////////////////////////
//		endTime                       	    //
//////////////////////////////////////////////////////
void 
PrintManager::endTime(double &comptime, double &dumptime)
{
    // output final simulation statistics
    outputFinalHeader();
    outputInfoCPU(comptime, dumptime);
    outputInfoMemory();
    outputFinalFooter();
}//end endTime

//////////////////////////////////////////////////////
//		outputFinalHeader       	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputFinalHeader()
{
    // get the final time and date
    auto end = std::chrono::high_resolution_clock::now();
    std::time_t endtime = std::chrono::system_clock::to_time_t(end);

    cout << "\n\n";
    cout << "\tCompleted simulation on " << std::ctime(&endtime) << std::endl;

    cout << "===================================================================\n";
    cout << "                ** BEGAN FINAL SIMULATION OUTPUT **\n\n";
}//end outputFinalHeader

//////////////////////////////////////////////////////
//		outputInfoCPU	        	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputInfoCPU(double &comptime, double &dumptime)
{
    // include dump time for outputting the extracted section
    dumptime += dumpsect_;

    // get information about the computing node being used
    cout << SecondEleWidth_ << "CPU time statistics:\n";
    cout << std::fixed;
    cout << std::setw(45) << std::right 
         << "Final simulation time: " 
         << std::setw(15) << setprecision(5) << domainMgr_->timet_ << " sec \n";
    cout << std::setw(45) << std::right
         << "Total compute time: "
	 << std::setw(15) << setprecision(5) << comptime << " sec \n";
    cout << std::setw(45) << std::right
         << "Total dump time: "
	 << std::setw(15) << setprecision(5) << dumptime << " sec \n";
    if(meshObj_->outputsection_)
    {
        cout << std::setw(45) << std::right
             << "Total section dump time: "
	     << std::setw(15) << setprecision(5) << dumpsect_ << " sec \n";
    }//end if
    cout << "\n\n";
}//end outputInfoCPU

//////////////////////////////////////////////////////
//		outputInfoMemory		    //
//////////////////////////////////////////////////////
void 
PrintManager::outputInfoMemory()
{
    // define some local variables
    double availMem = 0.0;
    double localMem = reportMemoryUsed();
    double localMaxMemUsed = reportMaxMemoryUsed();

    // determine amount of mermory available on the CPU
    getMemoryAvailable(availMem);

    // output information about CPU memory
    cout << std::fixed;
    cout << SecondEleWidth_ << "Memory statistics:\n";
    cout << std::setw(45) << std::right << "Total memory used: "
         << std::setw(15) << std::right << std::setprecision(5) << localMem << " MB\n";
    cout << std::setw(45) << std::right << "Max. memory used at one time: "
         << std::setw(15) << std::right << std::setprecision(5) << localMaxMemUsed << " MB\n";
    cout << std::setw(45) << std::right << "Total available memory on this node: "
         << std::setw(15) << std::right << std::setprecision(5) << availMem << " MB\n";
    cout << "\n";
}//end outputInfoMemory

//////////////////////////////////////////////////////
//		reportMemoryUsed	       	    //
//////////////////////////////////////////////////////
double
PrintManager::reportMemoryUsed()
{
    double memoryUsed = 0.0;
    #if defined(__APPLE__)
        struct task_basic_info t_info;
        mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
        
        if(KERN_SUCCESS != task_info(mach_task_self(),
            TASK_BASIC_INFO,
            reinterpret_cast<task_info_t>(&t_info),
            &t_info_count))
        {
            return 0;
        }//end if
        memoryUsed = t_info.resident_size/1.0e6;
    #elif defined(PROCFS)
        unsigned long rss_pages = 0;
        std::ifstream proc_stat("/proc/self/stat");
        if (proc_stat)
        {
            std::string buf;
            proc_stat
               >> buf >> buf >> buf >> buf >> buf
               >> buf >> buf >> buf >> buf >> buf
               >> buf >> buf >> buf >> buf >> buf
               >> buf >> buf >> buf >> buf >> buf
               >> buf >> buf >> buf
               >> rss_pages;
            proc_stat.close();
        }//end if
        memoryUsed = rss_pages*sysconf(_SC_PAGESIZE)/1.0e6;
    #endif
    return memoryUsed;
}//reportMemoryUsed

//////////////////////////////////////////////////////
//		reportMaxMemoryUsed	       	    //
//////////////////////////////////////////////////////
double
PrintManager::reportMaxMemoryUsed()
{
    // define some local variables
    double memoryUsed = 0.0;
    struct rusage usage;

    getrusage(RUSAGE_SELF, &usage);
    memoryUsed = usage.ru_maxrss;       // reported in kilobytes
    #if defined(__APPLE__)
        memoryUsed = memoryUsed/1.0e6;
    #elif defined(PROCFS)
        memoryUsed = memoryUsed/1000.0; // reported in megabytes
    #endif
return memoryUsed;
}//reportMaxMemoryUsed

//////////////////////////////////////////////////////
//		getMemoryAvailable             	    //
//////////////////////////////////////////////////////
void 
PrintManager::getMemoryAvailable(double &avail)
{
    avail = 0.0;
    #if defined(__APPLE__)
        int mib[2] = { CTL_HW, HW_MEMSIZE };
        u_int namelen = sizeof(mib) / sizeof(mib[0]);
        uint64_t size;
        size_t len = sizeof(size);
        sysctl(mib, namelen, &size, &len, NULL, 0);
        avail = size;
        const double bytes_in_MB = 1e6;
        avail /= bytes_in_MB; // output in MB
    #elif defined (PROCFS)
        std::string memtotal;
        std::string line(128,'\0');
        
        // Read available memory size data from /proc/meminfo 
        std::ifstream proc_meminfo("/proc/meminfo");
        if (!proc_meminfo) return;
        
        while (memtotal.empty())
        {
            if(!std::getline(proc_meminfo, line)) 
            {
              proc_meminfo.close();
              return;
            }
            
            // Find MemTotal
            else if (line.substr(0, 9) == "MemTotal:")
            {
              memtotal = line.substr(10);
              std::istringstream iss(memtotal);
              iss >> avail;
              avail *= 1000.0;
            }//end if
        }//end while
        const double bytes_in_MB = 1e6;
        avail /= bytes_in_MB; // output in MB
        proc_meminfo.close();
    #endif
}//end getMemoryAvailable

//////////////////////////////////////////////////////
//		outputFinalFooter       	    //
//////////////////////////////////////////////////////
void 
PrintManager::outputFinalFooter()
{
    cout << "                ** ENDED FINAL SIMULATION OUTPUT **\n";
    cout << "===================================================================\n";
}//end outputFinalFooter

