#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include "time.h"
#include <algorithm>
#include <iterator>
#include <fstream>
#include <omp.h>
#include <sys/stat.h>

// user-defined headers
#include <tecWriter.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
tecWriter::tecWriter(Mesh *meshObj,
                     DomainManager *domainMgr,
                     SolutionVariableManager *solVarMgr,
                     BoundCondManager *bcMgr,
                     CFDSolverManager *cfdSolvMgr,
                     FreeSurfaceManager *freesurMgr,
                     string &FileOut,
                     string &directoryName,
                     string &baseOutName)
                     : PrintManager(meshObj,
                                    domainMgr,
                                    solVarMgr,
                                    bcMgr,
                                    cfdSolvMgr,
                                    freesurMgr,
                                    FileOut,
                                    directoryName,
                                    baseOutName)

{
    mkdir(directoryName_.c_str(), 0777);
    isinit_ = true;
}

//////////////////////////////////////////////////////
//		visualOut                      	    //
//////////////////////////////////////////////////////
void 
tecWriter::visualOut()
{
    writeTec_openFile();
    writeTec_data();
    writeTec_closeFile();
}//end visualOut

//////////////////////////////////////////////////////
//		writeTec_openFile              	    //
//////////////////////////////////////////////////////
void 
tecWriter::writeTec_openFile()
{
        outputFile_.open(FileOut_, ios::out | ios::app);
        
        if(isinit_)
        {
            outputFile_ << "TITLE = \"Thermo-Capillary Flow in Laser-Generated Melt Pool\"\n";
            outputFile_ << "VARIABLES = ";
                for (int i = 0; i<outScalarNames_.size(); ++i)
                {
                       outputFile_ << "\"" << outScalarNames_[i]  << "\", ";
                }
            outputFile_ << "VARIABLES = \"X\", \"Y\", \"Z\"\n";

            isinit_ = false;
        }//end if
}//end writeTec_openFile

//////////////////////////////////////////////////////
//		writeTec_data              	    //
//////////////////////////////////////////////////////
void 
tecWriter::writeTec_data()
{
    if(meshObj_->tectype_ == "Tec")
    {
        outputTecFile();
    }
    else if(meshObj_->tectype_ == "Custom")
    {
        outputCustomFile();
    }
    else if(meshObj_->tectype_ == "UserDefined")
    {
        outputUserDefinedFile();
    }
    else if(meshObj_->tectype_ == "Surface")
    {
        outputSurfaceFile();
    }
}//end writeTec_data

//////////////////////////////////////////////////////
//		writeTec_closeFile             	    //
//////////////////////////////////////////////////////
void 
tecWriter::writeTec_closeFile()
{
    outputFile_.close();
}//end writeTec_closeFile

//////////////////////////////////////////////////////
//		outputTecFile			    //
//////////////////////////////////////////////////////
void 
tecWriter::outputTecFile()
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
    double xstart = domainMgr_->xstart_;

    // define local variables
    int i, j, k;
    
    // interpolate velocities at pressure nodes
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            for (i=1; i<nim1; ++i)
            {
                auvel_[(k*nj*ni)+(j*ni)+i] = (uvel_[(k*nj*ni)+(j*ni)+i] + uvel_[(k*nj*ni)+(j*ni)+i+1])*0.5;
                avvel_[(k*nj*ni)+(j*ni)+i] = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j+1)*ni)+i])*0.5;
                awvel_[(k*nj*ni)+(j*ni)+i] = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k+1)*nj*ni)+(j*ni)+i])*0.5;
                if(temp_[(k*nj*ni)+(j*ni)+i] <= tsolid);
                {
                    auvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    avvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    awvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
 
    // output Data
    outputFile_.width(4);
    outputFile_ << "ZONE  I= "<< nim2 << "J= " << njm2 << "K= " << nkm2 << "F=POINT \n";
    outputFile_.width(15);
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            for (i=1; i<nim1; ++i)
            {
                outputFile_ << setprecision(4) << scientific 
                            << 1000.0*(x_[i]-xstart) << 1000.0*y_[j] << 1000.0*z_[k] 
                            << auvel_[(k*nj*ni)+(j*ni)+i] << avvel_[(k*nj*ni)+(j*ni)+i] << awvel_[(k*nj*ni)+(j*ni)+i] 
                            << temp_[(k*nj*ni)+(j*ni)+i] << "\n";
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end outputTecFile

//////////////////////////////////////////////////////
//		outputCustomFile               	    //
//////////////////////////////////////////////////////
void 
tecWriter::outputCustomFile()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    double tsolid = domainMgr_->tsolid_;
    double timet = domainMgr_->timet_;
    double a_oxygen = domainMgr_->a_oxygen_;
    
    // grab some stuff from the mesh
    double delt = meshObj_->delt_;
    
    // define local variables
    int i, j, k;

    gridx = 0;
    gridy = 0;
    gridz = 0;

    for (i=1; i<nim1; ++i)
    {
        //if(x_[i] >= 0.005 && x_[i] <= 0.010)
        //{
            gridx++;
        //}
    }//end for(i)

    for (j=1; j<njm1; ++j)
    {

        //if(y_[j] >= 0.005 && y_[j] <= 0.010)
        //{
            gridy++;
        //}
    }//end for(j)

    for (k=1; k<nkm1; ++k)
    {

        //if(z_[k] >= 0.0 && z_[k] <= 0.003)
        //{
            gridz++;
        //}
    }//end for(k)

    // interpolate velocities at pressure nodes
    #pragma omp parallel for collapse(3)
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            for (i=1; i<nim1; ++i)
            {
                auvel_[(k*nj*ni)+(j*ni)+i] = (uvel_[(k*nj*ni)+(j*ni)+i] + uvel_[(k*nj*ni)+(j*ni)+i+1])*0.5;
                avvel_[(k*nj*ni)+(j*ni)+i] = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j+1)*ni)+i])*0.5;
                awvel_[(k*nj*ni)+(j*ni)+i] = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k+1)*nj*ni)+(j*ni)+i])*0.5;
                if(temp_[(k*nj*ni)+(j*ni)+i] <= tsolid) 
                {
                    auvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    avvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    awvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
 
    // output Data
    outputFile_.width(5);
    outputFile_ << "ZONE T = \"XYZ " << timet << "\"" << "I= "<< gridx << "J= " << gridy << "K= " << gridz << "F=POINT \n";
    outputFile_.width(15);
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            for (i=1; i<nim1; ++i)
            {
                //outputtemp_ = min(temp_[(k*nj*ni)+(j*ni)+i],3100.0);
                outputtemp_ = temp_[(k*nj*ni)+(j*ni)+i];
                //if(x_[i] >= 0.005 && x_[i] <= 0.010 && y_[j] >= 0.005 && y_[j] <= 0.010 && z_[k] >= 0.00 && z_[k] <= 0.003)
                //{
                outputFile_ << setprecision(4) << scientific 
                            << x_[i] << y_[j] << z_[k] 
                            << auvel_[(k*nj*ni)+(j*ni)+i] << avvel_[(k*nj*ni)+(j*ni)+i] << awvel_[(k*nj*ni)+(j*ni)+i] 
                            << outputtemp_ << grad_[(k*nj*ni)+(j*ni)+i] << rate_[(k*nj*ni)+(j*ni)+i] 
                            << bcMgr_->getSurfaceTensionGradient(temp_[(k*nj*ni)+(j*ni)+i], a_oxygen) << "\n";
                //}
            }//end for(i)
        }//end for(j)
    }//end for(k)

}//end outputCustomFile

//////////////////////////////////////////////////////
//		outputOtherFile               	    //
//////////////////////////////////////////////////////
void 
tecWriter::outputOtherFile()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nkm1 = domainMgr_->nkm1_;

    // define some local variables
    double area = bcMgr_->width_*bcMgr_->alen_;

    //speOut = fopen("otherFile.txt", "w");
    outputFile_ << area << "\n";
    outputFile_ << ni << " " << nj << "\n";

    outputFile_.width(9);
    for (int j=0; j<nj; ++j)
    {
        for (int i=0; i<ni; ++i)
        {
            outputFile_ << setprecision(3) << x_[i] << y_[j] << temp_[(nkm1*nj*ni)+(j*ni)+i] << "\n";
        }//end for(i)
    }//end for(j)
}//end outputOtherFile

//////////////////////////////////////////////////////
//		outputUserDefinedFile          	    //
//////////////////////////////////////////////////////
void 
tecWriter::outputUserDefinedFile()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    double tsolid = domainMgr_->tsolid_;
    double timet = domainMgr_->timet_;

    // grab some stuff from the mesh
    double delt = meshObj_->delt_;

    // define some local variables
    int i, j, k;
    double gxmin, gxmax, gymin, gymax, gzmin, gzmax;

    // define region of output
    gxmin = 2.5e-3;  
    gxmax = 3.5e-3;
    gymin = 0.5e-3;
    gymax = 0.75e-3;
    gzmin = 0.5e-3;
    gzmax = 0.6e-3;

    if(bcMgr_->powerindicator_ == 1)
    {
        // reserve transient datas every 10 time steps(1e-4)
        if((int)(timet/delt) % 2 != 0) 
            return;
    }
    else
    {
        // reserve transient datas every 40 time steps (0.01)
        if((int)(timet/delt) % 2 != 0) 
            return;
    }//end if

    gridx = 0;
    gridy = 0;
    gridz = 0;

    for (i=0; i<ni; ++i)
    {
        if(x_[i] >= gxmin && x_[i] <= gxmax)
            gridx = gridx + 1;
    }//end for(i)

    for (j=0; j<nj; ++j)
    {
        if(y_[j] >= gymin && y_[j] <= gymax)
            gridy = gridy + 1;
    }//end for(j)

    for (k=0; k<nkm1; ++k)
    {
        if(z_[k] >= gzmin && z_[k] <= gzmax)
            gridz = gridz + 1;
    }//end for(k)

    #pragma omp parallel for collapse(3)
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            for (i=1; i<nim1; ++i)
            {
                auvel_[(k*nj*ni)+(j*ni)+i] = (uvel_[(k*nj*ni)+(j*ni)+i] + uvel_[(k*nj*ni)+(j*ni)+i+1])*0.5;
                avvel_[(k*nj*ni)+(j*ni)+i] = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j+1)*ni)+i])*0.5;
                awvel_[(k*nj*ni)+(j*ni)+i] = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k+1)*nj*ni)+(j*ni)+i])*0.5;
                if(temp_[(k*nj*ni)+(j*ni)+i] <= tsolid) 
                {
                    auvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    avvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    awvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
 
    //----top plane
    for (j=1; j<njm1; ++j)
    {
        for (i=1; i<nim1; ++i)
        {
            auvel_[(nkm1*nj*ni)+(j*ni)+i] = (uvel_[(nkm1*nj*ni)+(j*ni)+i] + uvel_[(nkm1*nj*ni)+(j*ni)+i+1])*0.5;
            avvel_[(nkm1*nj*ni)+(j*ni)+i] = (vvel_[(nkm1*nj*ni)+(j*ni)+i] + vvel_[(nkm1*nj*ni)+((j+1)*ni)+i])*0.5;
            if(temp_[(nkm1*nj*ni)+(j*ni)+i] <= tsolid) 
            {
                auvel_[(nkm1*nj*ni)+(j*ni)+i] = 0.0;
                avvel_[(nkm1*nj*ni)+(j*ni)+i] = 0.0;
            }//end if
        }//end for(i)
    }//end for(j)

    //----symmetry plane
    for (k=1; k<nkm1; ++k)
    {
        for (i=1; i<nim1; ++i)
        {
            auvel_[(k*nj*ni)+(0*ni)+i] = (uvel_[(k*nj*ni)+(0*ni)+i] + uvel_[(k*nj*ni)+(0*ni)+i+1])*0.5;
            awvel_[(k*nj*ni)+(0*ni)+i] = (wvel_[(k*nj*ni)+(0*ni)+i] + wvel_[((k+1)*nj*ni)+(0*ni)+i])*0.5;
            if(temp_[(k*nj*ni)+(0*ni)+i] <= tsolid) 
            {
                auvel_[(k*nj*ni)+(0*ni)+i] = 0.0;
                awvel_[(k*nj*ni)+(0*ni)+i] = 0.0;
            }//end if
        }//end for(i)
    }//end for(k)

    //----left plane
    for (k=1; k<nkm1; ++k)
    {
        for (j=1; j<njm1; ++j)
        {
            avvel_[(k*nj*ni)+(j*ni)+0] = (vvel_[(k*nj*ni)+(j*ni)+0] + vvel_[(k*nj*ni)+(j*ni)+i+1])*0.5;
            awvel_[(k*nj*ni)+(j*ni)+0] = (wvel_[(k*nj*ni)+(j*ni)+0] + wvel_[((k+1)*nj*ni)+(j*ni)+0])*0.5;
            if(temp_[(k*nj*ni)+(j*ni)+0] <= tsolid) 
            {
                avvel_[(k*nj*ni)+(j*ni)+0] = 0.0;
                awvel_[(k*nj*ni)+(j*ni)+0] = 0.0;
            }//end if
        }//end for(j)
    }//end for(k)

    //for (k = nk,1,-1
    //{
        //for (j=1; j<nj; ++j)
        //{
            //if(y_[j] <= alasrb*1.1 && z_[nkm1] - z_[k] <= 0.00035)
            //{
                //for (i=ni-1; i>=0 && i<ni; --i)
                //{
                    //if(temp_[(nkm1*nj*ni)+(j*ni)+i] <= tsolid && bcMgr_->ipoweroff == 1) 
                    //{
                        //active[(k*nj*ni)+(j*ni)+i] = 0;
                    //}
                    //else
                    //{
                        //break
                    //}
                //}
            //}
        //}
    //}
    
      //yanping type
    outputFile_ << "ZONE T =," << timet << "," << " I= "<< gridx << " J= " << gridy << " K= " << gridz << " F=POINT \n";
      //tecplot type
    //outputFile_ << "ZONE T = \"XYZ " << timet << "\"" << "I= "<< gridx << "J= " << gridy << "K= " << gridz << "F=POINT \n";
    outputFile_.width(15);
    for (k=0; k<nkm1; ++k)
    {
        for (j=0; j<nj; ++j)
        {
            for (i=0; i<ni; ++i)
            {
                //if(active[(k*nj*ni)+(j*ni)+i] == 1) 
                //{
                    //showtemp[(k*nj*ni)+(j*ni)+i] = temp_[(k*nj*ni)+(j*ni)+i];
                //}
                //else
                //{
                    //showtemp[(k*nj*ni)+(j*ni)+i] = -1;
                //}

                if(x_[i] >= gxmin && x_[i] <= gxmax && y_[j] >= gymin && y_[j] <= gymax && z_[k] >= gzmin && z_[k] <= gzmax)
                {
                    outputtemp_ = min(temp_[(k*nj*ni)+(j*ni)+i],3100.0);

                    outputFile_ << setprecision(4) << scientific 
                                << x_[i] << y_[j] << z_[k] 
                                << auvel_[(k*nj*ni)+(j*ni)+i] << avvel_[(k*nj*ni)+(j*ni)+i] << awvel_[(k*nj*ni)+(j*ni)+i] 
                                << outputtemp_ << "\n";
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end outputUserDefinedFile

//////////////////////////////////////////////////////
//		outputSurfaceFile          	    //
//////////////////////////////////////////////////////
void 
tecWriter::outputSurfaceFile()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = domainMgr_->nk_;
    int nim1 = domainMgr_->nim1_;
    int njm1 = domainMgr_->njm1_;
    int nkm1 = domainMgr_->nkm1_;
    double a_oxygen = domainMgr_->a_oxygen_;
    double tsolid = domainMgr_->tsolid_;
    double timet = domainMgr_->timet_;

    // grab some stuff from the mesh
    double delt = meshObj_->delt_;

    // define local variables
    int i, j, k;

    if(bcMgr_->powerindicator_ == 0)
    {
        // reserve transient datas every 1000 time steps(1e-4)
        if((int)(timet/delt) % 10 != 0)        
            return;
    }
    else
    {
        // reserve transient datas every 40 time steps (0.01)
        if((int)(timet/delt) % 40 != 0)        
            return;
    }//end if

    //----top plane
    for (j=1; j<njm1; ++j)
    {
        for (i=1; i<nim1; ++i)
        {
            auvel_[(nkm1*nj*ni)+(j*ni)+i] = (uvel_[(nkm1*nj*ni)+(j*ni)+i] + uvel_[(nkm1*nj*ni)+(j*ni)+i+1])*0.5;
            avvel_[(nkm1*nj*ni)+(j*ni)+i] = (vvel_[(nkm1*nj*ni)+(j*ni)+i] + vvel_[(nkm1*nj*ni)+((j+1)*ni)+i])*0.5;
            if(temp_[(nkm1*nj*ni)+(j*ni)+i] <= tsolid) 
            {
                auvel_[(nkm1*nj*ni)+(j*ni)+i] = 0.0;
                avvel_[(nkm1*nj*ni)+(j*ni)+i] = 0.0;
            }//end if
        }//end for(i)
    }//end for(j)

    k = nkm1;

    // output data
    outputFile_ << "ZONE T = \"XYZ " << timet << "\"" << "I= "<< ni << "J= " << nj << "K= " << 1 << "F=POINT \n";
    outputFile_.width(15);
    for (j=0; j<nj; ++j)
    {
        for (i=0; i<ni; ++i)
        {
            outputFile_ << setprecision(4) << scientific 
                        << x_[i] << y_[j] << z_[k] 
                        << auvel_[(k*nj*ni)+(j*ni)+i] << avvel_[(k*nj*ni)+(j*ni)+i] << awvel_[(k*nj*ni)+(j*ni)+i] 
                        << temp_ << grad_[(k*nj*ni)+(j*ni)+i] << rate_[(k*nj*ni)+(j*ni)+i] 
                        << bcMgr_->getSurfaceTensionGradient(temp_[(k*nj*ni)+(j*ni)+i], a_oxygen) << "\n";
        }//end for(i)
    }//end for(j)
}//end outputSurfaceFile
