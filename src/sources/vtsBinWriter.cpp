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
#include <vtsBinWriter.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
vtsBinWriter::vtsBinWriter(Mesh *meshObj,
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
    pvdName_ = baseOutName_ + ".pvd";
    nDim_ = 3;
}

//////////////////////////////////////////////////////
//		visualOut			    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::visualOut()
{
    // set up some initial stuff 
    offSetCtr_ = 0;
    byteCtr_ = 0;
    outputFile_.open(FileOut_, ios::out | ios::binary);
 
    // define number of nodes here and set up initial headers
    nPoints_ = domainMgr_->ni_*domainMgr_->nj_*(gridz_-1);

    writeVTS_header();
    
    writeVTS_coordsHeader();
    
    writeVTS_coordsEnd();
    
    writeVTS_pointDataHeader();
    
    // set up nodal data outputs (scalar)
    for (int iout = 0; iout < outScalarNames_.size(); iout++)
    {
        string outname = outScalarNames_[iout];
        writeVTS_pointOutHeader(outname);
    }//end for(iout)

    // this provides the header for layer addition
    if(meshObj_->powderbed_)
    {
        string outname = "vtkGhostType";
        writeVTS_pointOutHeader(outname);
    }//end if
    
    // set up nodal data outputs (vector)
    for (int iout = 0; iout < outVectorNames_.size(); iout++)
    {
        string outname = outVectorNames_[iout];
        writeVTS_pointOutVecHeader(outname);
    }//end for(iout)
    
    writeVTS_pointDataEnd();
    
    // Kevochan: currently not used
    writeVTS_cellDataHeader();
    
    // Kevochan: currently not used
    //writeVTS_cellOutHeader();
    
    writeVTS_cellDataEnd();
    
    writeVTS_end();
    
    writeVTS_appendHeader();
    
    transferVelNodes();
    
    // append data
    writeVTS_coordinates();
    
    for (int iout = 0; iout < outScalarNames_.size(); iout++)
    {
        string outname = outScalarNames_[iout];
        double *pointDataOut = solVarMgr_->mapSolutionVars_[outname];
        writeVTS_pointData(pointDataOut);
    }//end for(iout)

    // this handles visualization present layers
    if(meshObj_->powderbed_)
        writeVTS_GhostPoints();
    
    for (int iout = 0; iout < outVectorNames_.size(); iout++)
    {
        string outname = outVectorNames_[iout];
        double *pointDataOut = solVarMgr_->mapSolutionVars_[outname];
        writeVTS_pointVecData(pointDataOut);
    }//end for(iout)
    
    // Kevochan: currently not used
    //writeVTS_cellData();
    
    // finish appending data
    writeVTS_appendEnd();
    
    // write new pvd file
    initializePVD();
    appendPVD();
    closePVD();

    // close the output file
    outputFile_.close();
}//end visualOut

//////////////////////////////////////////////////////
//		writeVTS_coordinates		    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_coordinates()
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;

    // x, y, z coordintes are the intial appended data
    byteCtr_ = sizeof(double) * nPoints_ * 3;
    outputFile_.write((char*) &byteCtr_, sizeof(int));

    // loop over all the coordinates, except top boundary
    for (int k=0; k<gridz_-1; k++)
    {
        for (int j=0; j<domainMgr_->nj_; j++)
        {
            for (int i=0; i<domainMgr_->ni_; i++)
            {
                outputFile_.write((char*) &x_[i], sizeof(double));
                outputFile_.write((char*) &y_[j], sizeof(double));
                if(!meshObj_->energyfreesurface_)
                {
                    outputFile_.write((char*) &z_[k], sizeof(double));
                }
                else
                {
                    outputFile_.write((char*) &zr_[(k*nj*ni)+(j*ni)+i], sizeof(double));
                }
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end writeVTS_coordinates

//////////////////////////////////////////////////////
//		writeVTS_coordsHeader		    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_coordsHeader()
{
    outputFile_ << "<Points> ";
    
    outputFile_ << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset = \"" 
                << offSetCtr_ << "\" /> ";
    
    offSetCtr_ += sizeof(double) * nPoints_ * 3 + sizeof(offSetCtr_);
}//end writeVTS_coordsHeader

//////////////////////////////////////////////////////
//		writeVTS_coordsEnd      	    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_coordsEnd()
{
    outputFile_  << "</Points> ";
}//end writeVTS_coordsEnd

//////////////////////////////////////////////////////
//		writeVTS_pointData		    //
//////////////////////////////////////////////////////
void
vtsBinWriter::writeVTS_pointData(double *pointVec)
{
    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nk = gridz_;

    byteCtr_ = sizeof(double) * nPoints_;
    outputFile_.write((char*) &byteCtr_, sizeof(int));

    // write out point data
    for (int k=0; k<nk-1; k++)
    {
        for (int j=0; j<nj; j++)
        {
            for (int i=0; i<ni; i++)
            {
                outputFile_.write((char*) &pointVec[(k*nj*ni)+(j*ni)+i], sizeof(double));
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end writeVTS_pointData

//////////////////////////////////////////////////////
//		writeVTS_GhostPoints		    //
//////////////////////////////////////////////////////
void
vtsBinWriter::writeVTS_GhostPoints()
{
    // write out point data
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;
    int nkm1 = domainMgr_->nkm1_;

    byteCtr_ = sizeof(unsigned char) * nPoints_;
    outputFile_.write((char*) &byteCtr_, sizeof(int));

    unsigned char visible = 0;
    unsigned char hidden = 2;

    for (int k=0; k<gridz_-1; k++)
    {
        for (int j=0; j<nj; j++)
        {
            for (int i=0; i<ni; i++)
            {
                if(k < nkm1)
                    outputFile_.write((char*) &visible, sizeof(unsigned char));
                else
                    outputFile_.write((char*) &hidden, sizeof(unsigned char));
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end writeVTS_GhostPoints

//////////////////////////////////////////////////////
//		writeVTS_header			    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_header()
{
    // define internal node indecies
    int x1 = 0; 
    int y1 = 0; 
    int z1 = 0; 
    int x2 = domainMgr_->nim1_; 
    int y2 = domainMgr_->njm1_; 
    //int z2 = domainMgr_->nkm2_; 
    int z2 = gridz_-2; 

    // output necessary data
    outputFile_ << "<?xml version=\"1.0\"?> ";
    outputFile_ << "<VTKFile type=\"StructuredGrid\" version=\"2.1\" byte_order=\"LittleEndian\"> ";
    
    outputFile_ << "<StructuredGrid WholeExtent=\""
                                 << x1 << " " << x2 << " " 
                                 << y1 << " " << y2 << " " 
                                 << z1 << " " << z2 << "\"> ";
    
    // local values are just min and maxes on the thread
    outputFile_ << "<Piece Extent=\""
                                 << x1 << " " << x2 << " " 
                                 << y1 << " " << y2 << " " 
                                 << z1 << " " << z2 << "\"> ";
}//end writeVTS_header

//////////////////////////////////////////////////////
//		writeVTS_end   			    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_end()
{
    outputFile_ << "</Piece> ";
    outputFile_ << "</StructuredGrid> ";
}//end writeVTS_end

//////////////////////////////////////////////////////
//		writeVTS_pointDataHeader	    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_pointDataHeader()
{
    outputFile_ <<"<PointData> ";
}//end writeVTS_pointDataHeader

//////////////////////////////////////////////////////
//		writeVTS_pointDataEnd	   	    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_pointDataEnd()
{
    outputFile_ <<"</PointData> ";
}//end writeVTS_pointDataEnd

//////////////////////////////////////////////////////
//		writeVTS_pointOutHeader	    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_pointOutHeader(string &outScalarName)
{
    if(outScalarName == "vtkGhostType")
    {
        outputFile_ << "<DataArray type=\"UInt8\" Name=\"" 
                    << outScalarName <<"\" format=\"appended\" offset=\""
                    << offSetCtr_ << "\" />  ";
        offSetCtr_ += sizeof(unsigned char) * nPoints_ + sizeof(offSetCtr_);
    }
    else
    {
        outputFile_ << "<DataArray type=\"Float64\" Name=\"" 
                    << outScalarName <<"\" format=\"appended\" offset=\""
                    << offSetCtr_ << "\" />  ";
        offSetCtr_ += sizeof(double) * nPoints_ + sizeof(offSetCtr_);
    }//end if
}//end writeVTS_pointOutHeader

//////////////////////////////////////////////////////
//		writeVTS_pointOutEnd	    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_pointOutEnd()
{
    string FifthEleWidth  = "        ";
    
    outputFile_ << FifthEleWidth <<"</DataArray>\n";
}//end writeVTS_pointOutEnd

//////////////////////////////////////////////////////
//		writeVTS_appendHeader	    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_appendHeader()
{
    outputFile_ <<"<AppendedData encoding=\"raw\">_";
}//end writeVTS_appendHeader

//////////////////////////////////////////////////////
//		writeVTS_appendEnd	    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_appendEnd()
{
    outputFile_ << "\n";
    outputFile_ << "</AppendedData> </VTKFile>";
}//end writeVTS_appendEnd

//////////////////////////////////////////////////////
//		initializePVD 		    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::initializePVD()
{
    outputPVD_.open(pvdName_, ios::out | ios::binary);
    string SecondEleWidth = "  ";

    // Write Headers
    outputPVD_ << "<?xml version=\"1.0\"?>\n";
    outputPVD_ << "<VTKFile type=\"Collection\" version=\"2.1\" byte_order=\"LittleEndian\">\n";
    outputPVD_ << SecondEleWidth << "<Collection>\n";
}//end initializePVD

//////////////////////////////////////////////////////
//		appendPVD 	    	            //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::appendPVD()
{
    string ThirdEleWidth = "   ";
    string timeStamp = to_string(domainMgr_->timet_);
    
    timeStampList_.push_back(domainMgr_->timet_); 
    pvtsOutList_.push_back(FileOut_);
    
    for (int i = 0; i < timeStampList_.size(); i++)
    {
        outputPVD_ << ThirdEleWidth << "<DataSet timestep=\"" << 
        timeStampList_[i] << "\" group=\"\" part=\"0\" file=\"" <<
        pvtsOutList_[i] <<"\"/>\n";
    }//end for(i)
}//end appendPVD

//////////////////////////////////////////////////////
//		closePVD 	    	            //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::closePVD()
{
    string SecondEleWidth = "  ";
    
    outputPVD_ << SecondEleWidth << "</Collection>\n";
    outputPVD_ << "</VTKFile>\n";
    
    outputPVD_.close();
}//end closePVD

//////////////////////////////////////////////////////
//		writeVTS_cellDataHeader		    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_cellDataHeader()
{
    outputFile_ <<"<CellData> ";
}//end writeVTS_cellDataHeader

//////////////////////////////////////////////////////
//		writeVTS_cellOutHeader	    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_cellOutHeader()
{
    outputFile_ << "<DataArray type=\"Float64\" Name=\"consolidFrac\" format=\"appended\" offset=\""
                << offSetCtr_ << "\" />  ";
    offSetCtr_ += sizeof(double) * nCells_ + sizeof(offSetCtr_);
}//end writeVTS_cellOutHeader

//////////////////////////////////////////////////////
//		writeVTS_cellDataEnd	   	    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_cellDataEnd()
{
    outputFile_ <<"</CellData> ";
}//end writeVTS_cellDataEnd

//////////////////////////////////////////////////////
//		writeVTS_cellData	   	    //
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_cellData()
{
    // write out point data
    byteCtr_ = sizeof(double) * nCells_;
    outputFile_.write((char*) &byteCtr_, sizeof(int));
    
    // loop over the contol volumes for appropiate output data

}//end writeVTS_cellData

//////////////////////////////////////////////////////
//		writeVTS_pointOutVecHeader    	    //	
//////////////////////////////////////////////////////
void 
vtsBinWriter::writeVTS_pointOutVecHeader(string &outVectorName)
{
    outputFile_ << "<DataArray type=\"Float64\" Name=\"" 
                << outVectorName <<"\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
                << offSetCtr_ << "\" />  ";
    offSetCtr_ += sizeof(double) * nPoints_*nDim_ + sizeof(offSetCtr_);
}//end writeVTS_pointOutVecHeader

//////////////////////////////////////////////////////
//		writeVTS_pointVecData		    //
//////////////////////////////////////////////////////
void
vtsBinWriter::writeVTS_pointVecData(double *pointVec)
{
    // write out point data
    byteCtr_ = sizeof(double) * nPoints_ * nDim_;
    outputFile_.write((char*) &byteCtr_, sizeof(int));
    
    for (int i = 0; i < nPoints_; i++)
    {
        int nodeOffSetG = i * nDim_;
        for (int j = 0; j < nDim_; j++)
        {
            outputFile_.write((char*) &pointVec[nodeOffSetG + j], sizeof(double));
        }//end for(j)
    }//end for(i)
}//end writeVTS_pointVecData

//////////////////////////////////////////////////////
//		transferVelNodes		    //
//////////////////////////////////////////////////////
void
vtsBinWriter::transferVelNodes()
{
    // define some local variables
    int i, j, k;

    // grab some stuff from the domain
    int ni = domainMgr_->ni_;
    int nj = domainMgr_->nj_;

    // shift velocity nodes from cell faces to cell centers 
    for (k=1; k<domainMgr_->nkm1_; k++)
    {
        for (j=1; j<domainMgr_->njm1_; j++)
        {
            for (i=1; i<domainMgr_->nim1_; i++)
            {
                auvel_[(k*nj*ni)+(j*ni)+i] = (uvel_[(k*nj*ni)+(j*ni)+i] + uvel_[(k*nj*ni)+(j*ni)+i+1])*0.5;
                avvel_[(k*nj*ni)+(j*ni)+i] = (vvel_[(k*nj*ni)+(j*ni)+i] + vvel_[(k*nj*ni)+((j+1)*ni)+i])*0.5;
                awvel_[(k*nj*ni)+(j*ni)+i] = (wvel_[(k*nj*ni)+(j*ni)+i] + wvel_[((k+1)*nj*ni)+(j*ni)+i])*0.5;
                if(temp_[(k*nj*ni)+(j*ni)+i] <= domainMgr_->tsolid_) 
                {
                    auvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    avvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                    awvel_[(k*nj*ni)+(j*ni)+i] = 0.0;
                }//end if
            }//end for(i)
        }//end for(j)
    }//end for(k)
}//end transferVelNodes
