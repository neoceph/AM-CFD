// standard headers
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include "time.h"
#include <fstream>

// user defined headers
#include <Mesh.h>
#include <BoundCondManager.h>
#include <CFDSolverManager.h>
#include <DomainManager.h>
#include <FreeSurfaceManager.h>
#include <PrintManager.h>
#include <SolutionVariableManager.h>
#include <tecWriter.h>
#include <vtsBinWriter.h>

using namespace std;

int main(int arg, char *argv[])
{
    // make sure user provides an input file
    if (arg != 2)
    {
      cout << "Missing parameter: input file name" << endl;
      exit(1);
      return 0;
    }

    // define some variables for time-keeping
    auto begin  = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    string inputFile = argv[1];
    Mesh meshObj(inputFile);

    // get the mesh 
    meshObj.getDomainInfo();

    // get user parameters
    meshObj.assignParameters();
    meshObj.outputInputs(argv[1]);
     
    // create computational domain based on the input file 
    DomainManager *domainMgr = new DomainManager(&meshObj);

    // create viewing window for solution variables
    SolutionVariableManager *solVarMgr = new SolutionVariableManager();
     
    // create systems for solving and initialize them
    BoundCondManager *bcMgr = new BoundCondManager(&meshObj, domainMgr, solVarMgr);
    CFDSolverManager *cfdSolvMgr = new CFDSolverManager(&meshObj, domainMgr, 
                                                        solVarMgr, bcMgr);
     
    FreeSurfaceManager *freesurMgr = new FreeSurfaceManager(&meshObj, domainMgr, 
                                                            solVarMgr, bcMgr);
    // define some variables for outputting data
    PrintManager *printMgr;
    string extName, outFile;
    int outCt = 0;
    string directoryName = meshObj.directoryName_;
    string folderName = directoryName + "/";

    // set up output manager
    if(meshObj.istec_)
    {
        extName = ".tec";
        outFile = folderName + meshObj.outFileName_ + extName;
        printMgr = new tecWriter(&meshObj, domainMgr, solVarMgr,
                                 bcMgr, cfdSolvMgr, freesurMgr, outFile, 
                                 directoryName, meshObj.outFileName_);
    }
    else
    {
        extName = ".vts";
        outFile = folderName + meshObj.outFileName_ + "_" + 
                  to_string(outCt) + extName;
        printMgr = new vtsBinWriter(&meshObj, domainMgr, solVarMgr,
                                    bcMgr, cfdSolvMgr, freesurMgr, outFile, 
                                    directoryName, meshObj.outFileName_);
    }//end if
     
    // initialize the simulation
    printMgr->initialization();
     
    // check if the problem is steady-state
    if(domainMgr->steady_) 
        bcMgr->updateLaser();
    
    // set up the time trackers
    double dumptime = 0.0;
    double comptime = 0.0;
    double outtrack = 0.0;
    domainMgr->timet_ = 0.0;
    double simtime = meshObj.finaltime_;
    double *secttime = &printMgr->secttime_;

    // output data at initial time step
    printMgr->scanOutputData();
    printMgr->visualOut();
 
    // start the time loop and and simulation timer
    auto simstart = std::chrono::high_resolution_clock::now();
    while (domainMgr->timet_ < simtime)
    {
        // update timesteps
        domainMgr->timet_ += meshObj.delt_;
        outtrack += meshObj.delt_;
        
        // update the laser location for transient simulation
        if(!domainMgr->steady_) 
            bcMgr->updateLaser();
        
        // update the mesh according to the laser
        if(meshObj.movingmesh_)
            cfdSolvMgr->generateGrid();

        // iteratively solve conservation equations
        cfdSolvMgr->iterativeLoopSolve();

        // post-process the data to calculate G/R
        if(meshObj.solidificationparams_)
            cfdSolvMgr->getSolidificationParameters();

        // output an extracted section of the domain 
        if(meshObj.outputsection_)
            printMgr->outputSection();

        // determine the output for the simulation
        if(!domainMgr->steady_)
        {
            // push back the solution variables
            cfdSolvMgr->updateSolutions();

            // check for outputting data
            if(outtrack >= meshObj.outtime_)
            {
                // output for paraview
                outCt++;
                outFile = folderName + meshObj.outFileName_ + "_" + 
                                to_string(outCt) + extName;

                // begin outputting
                begin = std::chrono::high_resolution_clock::now();
                  // update the free surface via energy balence and bi-section method 
                if(meshObj.energyfreesurface_)
                    freesurMgr->execute();
                  // output solver information to file 
                printMgr->outputResiduals();
                  // print data for visualization
                printMgr->visualOut();
                end = std::chrono::high_resolution_clock::now();
                outtrack = 0.0;

                  // output progress to completion
                cout << "===============================================================\n";
                cout << left << 
                	"   Output at time: " << setw(43) << domainMgr->timet_ 
                			      << "|\n";
                cout << left << 
                	"   Percentage done: " << setw(42) << (domainMgr->timet_/simtime) * 100.0 
                			       << "|\n";
                cout << left << 
                	"   Timer for outputting files: " << setw(31) 
                     <<  1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()+(*secttime)
                			       << "|\n";
                cout << "===============================================================" << endl;
                dumptime += 1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
                *secttime = 0.0;
            }//end if
        }//end if(!steady)
    }//end while(timet)

    // set simulation to finished
    domainMgr->finish_ = true;

    // output velocity and temperature at surface and symmetry plane to file
    if(domainMgr->steady_)  
        printMgr->outputFinalSteady();
         
    // update final free surface
    if(meshObj.energyfreesurface_)
        freesurMgr->execute();
    
    // output final extracted section data
    if(meshObj.outputsection_)
        printMgr->outputSection();

    // output final visualization and solver residual information
    outCt++;
    outFile = folderName + meshObj.outFileName_ + "_" + to_string(outCt) + extName;
    printMgr->outputResiduals();
    printMgr->visualOut();

    // output the final melt pool geometry
    if(meshObj.outputgeom_)
        printMgr->outputFinalGeometry();

    // calculate total simulation time
    auto simend = std::chrono::high_resolution_clock::now();
    comptime = 1.0e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(simend-simstart).count();
    
    // output summary of simulation statistics
    printMgr->endTime(comptime, dumptime);

    // final clean-up (may change these to smart pointers later)
    delete bcMgr;
    delete cfdSolvMgr;
    delete freesurMgr;
    delete printMgr;
    delete solVarMgr;
    delete domainMgr;

    // end program
    return 0;
}//end main
