// standard headers
#include <numeric>
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

// user-defined headers
#include <DomainManager.h>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
DomainManager::DomainManager(Mesh *meshObj)
                             : meshObj_(meshObj)
                               
{
}

//////////////////////////////////////////////////////
//		assignInputValues	            //
//////////////////////////////////////////////////////
void
DomainManager::assignInputValues()
{
    // geometric parameters to create the mesh
    nzx_ = (int)meshObj_->x_grid_[0];
    xzone_ = &meshObj_->x_grid_[1];
    ncvx_ = &meshObj_->x_grid_[1 + nzx_];
    powrx_ = &meshObj_->x_grid_[1 + 2*nzx_];

    nzy_ = (int)meshObj_->y_grid_[0];
    yzone_ = &meshObj_->y_grid_[1];
    ncvy_ = &meshObj_->y_grid_[1 + nzy_];
    powry_ = &meshObj_->y_grid_[1 + 2*nzy_];

    nzz_ = (int)meshObj_->z_grid_[0];
    zzone_ = &meshObj_->z_grid_[1];
    ncvz_ = &meshObj_->z_grid_[1 + nzz_];
    powrz_ = &meshObj_->z_grid_[1 + 2*nzz_];

    // process parameters (single track)
    if(!meshObj_->udtoolpath_)
    {
        lasertime_ = meshObj_->proc_params_[0];
        xstart_ = meshObj_->proc_params_[1];
        ystart_ = meshObj_->proc_params_[2];
        scanvel_ = meshObj_->proc_params_[3];
        vector<double>().swap(meshObj_->proc_params_);
    }//end if

    // material properties
    dens_ = meshObj_->mat_props_[0];
    denl_ = meshObj_->mat_props_[1];
    viscos_ = meshObj_->mat_props_[2];
    tsolid_ = meshObj_->mat_props_[3];
    tliquid_ = meshObj_->mat_props_[4];
    hsmelt_ = meshObj_->mat_props_[5];
    hlfriz_ = meshObj_->mat_props_[6];
    acpa_ = meshObj_->mat_props_[7];
    acpb_ = meshObj_->mat_props_[8];
    acpl_= meshObj_->mat_props_[9];
    thconsa_ = meshObj_->mat_props_[10];
    thconsb_ = meshObj_->mat_props_[11];
    thconsc_ = meshObj_->mat_props_[12];
    thconla_ = meshObj_->mat_props_[13];
    thconlb_ = meshObj_->mat_props_[14];

    // under_relaxation parameters
    urfu_ = meshObj_->num_relax_[0];
    urfv_ = meshObj_->num_relax_[1];
    urfw_ = meshObj_->num_relax_[2];
    urfp_ = meshObj_->num_relax_[3];
    urfh_ = meshObj_->num_relax_[4];

    // boundary conditions: energy
    htci_ = meshObj_->bound_cond_energy_[0];
    htcj_ = meshObj_->bound_cond_energy_[1];
    htck1_ = meshObj_->bound_cond_energy_[2];
    htckn_ = meshObj_->bound_cond_energy_[3];
    emiss_ = meshObj_->bound_cond_energy_[4];
     
    // boundary conditions: pre-heat temperatures 
    tempwest_ = meshObj_->bound_cond_temp_[0]; 
    tempeast_ = meshObj_->bound_cond_temp_[1];
    tempnorth_ = meshObj_->bound_cond_temp_[2];
    tempbottom_ = meshObj_->bound_cond_temp_[3]; 
    temppreheat_ = meshObj_->bound_cond_temp_[4];
     
    // boundary conditions: momentum 
    beta_ = meshObj_->bound_cond_momentum_[0];
    dgdtp_ = meshObj_->bound_cond_momentum_[1];
    darcyKo_ = meshObj_->bound_cond_momentum_[2];
    a_oxygen_ = meshObj_->bound_cond_momentum_[3];

    // gaussian heat source parameters 
    if(meshObj_->gaussheatsource_)
    {
        alaspow_ = meshObj_->gauss_[0];
        alaseta_ = meshObj_->gauss_[1]; 
        alasrb_ = meshObj_->gauss_[2];
        alasfact_ = meshObj_->gauss_[3];
        vector<double>().swap(meshObj_->gauss_); 
    }//end if

    // evaporation parameters
    if(meshObj_->isevap_)
    {
        evaplatent_ = meshObj_->evap_params_[0]; 
        evaptemp_ = meshObj_->evap_params_[1]; 
        evapfact_ = meshObj_->evap_params_[2]; 
        molarmass_ = meshObj_->evap_params_[3]; 
        vector<double>().swap(meshObj_->evap_params_); 
    }//end if 
    
    // species solver parameters
    if(meshObj_->species_)
    {
        initconcentration_ = meshObj_->species_params_[0];
        massdiff_ = meshObj_->species_params_[1];
        kp_ = meshObj_->species_params_[2];
        urfc_ = meshObj_->species_params_[3];
        vector<double>().swap(meshObj_->species_params_); 
    }//end if
     
    // free-surface solver parameters
    if(meshObj_->energyfreesurface_)
    {
        initlam1_ = meshObj_->energyfree_params_[0];
        initlam2_ = meshObj_->energyfree_params_[1];
        surftens_ = meshObj_->energyfree_params_[2];
        maxiterenergy_ = meshObj_->energyfree_params_[3];
        tolenergy_ = meshObj_->energyfree_params_[4];
        maxiterbisect_ = meshObj_->energyfree_params_[5];
        tolbisect_ = meshObj_->energyfree_params_[6];
        urfree_ = meshObj_->energyfree_params_[7];
        vector<double>().swap(meshObj_->energyfree_params_); 
    }//end if
     
    // powder-bed parameters 
    if(meshObj_->powderbed_)
    {
        numlayers_ = (int)meshObj_->powderbed_params_[0];
        ncvzlayer_ = (int)meshObj_->powderbed_params_[1];
        layerheight_ = meshObj_->powderbed_params_[2];
        powderden_ = meshObj_->powderbed_params_[3];
        powdercpa_ = meshObj_->powderbed_params_[4];
        powdercpb_ = meshObj_->powderbed_params_[5];
        powderthcona_ = meshObj_->powderbed_params_[6];
        powderthconb_ = meshObj_->powderbed_params_[7];
        vector<double>().swap(meshObj_->powderbed_params_); 

        // add in the build mesh
        zzone_[nzz_-1] += numlayers_*layerheight_;
        ncvz_[nzz_-1] += numlayers_*ncvzlayer_;
    }//end if

    // volumetric heat source parameters 
    if(meshObj_->volheatsource_)
    {
        if (!meshObj_->volheatsourcetable_)
        {
            avolpow_ = meshObj_->volheatsource_params_[0];
            apowseta_ = meshObj_->volheatsource_params_[1];
            heatthick_ = meshObj_->volheatsource_params_[2];
            heatrb_ = meshObj_->volheatsource_params_[3];
            avolfact_ = meshObj_->volheatsource_params_[4];
            vector<double>().swap(meshObj_->volheatsource_params_); 
        }
        else 
        {
            getHeatSourceParameters();
            avolpow_ = heatSourceParameterValues[0][1];
            apowseta_ = heatSourceParameterValues[0][2];
            heatthick_ = heatSourceParameterValues[0][3];
            heatrb_ = heatSourceParameterValues[0][4];
            avolfact_ = heatSourceParameterValues[0][5];
            vector<double>().swap(meshObj_->volheatsource_params_); 
        }
    }//end if

    // output section parameters
    if(meshObj_->outputsection_)
    {
        gxmin_ = meshObj_->outputsection_params_[0];
        gymin_ = meshObj_->outputsection_params_[1];
        gzmin_ = meshObj_->outputsection_params_[2];
        gxmax_ = meshObj_->outputsection_params_[3];
        gymax_ = meshObj_->outputsection_params_[4];
        gzmax_ = meshObj_->outputsection_params_[5];
        toutsect_ = meshObj_->outputsection_params_[6];
        vector<double>().swap(meshObj_->outputsection_params_); 
    }//end if

    // determine if the analysis is transient or steady-state 
    // if time step less than 100, then transient-state simulation
    if(meshObj_->delt_ < 100.0) 
    {        
        steady_ = false;
    }
    else
    {
        steady_ = true;     
    }

    // calculate total number of control volumes
    int sumcvx = std::accumulate(&ncvx_[0], &ncvx_[0]+nzx_, 0);
    int sumcvy = std::accumulate(&ncvy_[0], &ncvy_[0]+nzy_, 0);
    int sumcvz = std::accumulate(&ncvz_[0], &ncvz_[0]+nzz_, 0);
     
    // assign the numer of nodes for entire domain including boundaries
    int numboundnodes = 2;
    ni_ = sumcvx + numboundnodes; 
    nj_ = sumcvy + numboundnodes; 
    nk_ = sumcvz + numboundnodes; 

    // identifying cell center coordinate array index 
    int count = 0;

    for (int k=0; k<nk_; k++)
    {
        for (int j=0; j<nj_; j++)
        {
            for (int i=0; i<ni_; i++)
            {
                coordinateIndex[count][0] = i;
                coordinateIndex[count][1] = j;
                coordinateIndex[count][2] = k;
                count+=1;
            }
        }
        
    }
    
    // define the boundary nodes
    nim1_ = ni_-1;
    njm1_ = nj_-1;
    nkm1_ = nk_-1;

    // define the last internal nodes
    nim2_ = ni_-2;
    njm2_ = nj_-2;
    nkm2_ = nk_-2;

    // caulate latent heat
    double cpavg = (acpa_*tsolid_ + acpb_ + acpl_)*0.5;
    hlcal_ = hsmelt_ + cpavg*(tliquid_ - tsolid_);
    hlatnt_ = hlfriz_ - hlcal_; 

    // clean up the mesh of some inputfile vectors
    vector<double>().swap(meshObj_->mat_props_);
    vector<double>().swap(meshObj_->num_relax_);
    vector<double>().swap(meshObj_->bound_cond_temp_);
    vector<double>().swap(meshObj_->bound_cond_energy_);
    vector<double>().swap(meshObj_->bound_cond_momentum_);
    vector<double>().swap(meshObj_->gauss_);
}//end assignInputValues

//////////////////////////////////////////////////////
//		getToolpath()			    //
//////////////////////////////////////////////////////
void
DomainManager::getToolpath()
{
    ifstream file;
    file.open(meshObj_->toolFileName_.c_str());
    string line;
    if (file.is_open())
    {
        while(getline(file, line))
        {
            istringstream lines(line);
            vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
            vector<double> txyz(coords.begin(), coords.begin() + 4);
     
            int state = (int)coords[4];
            laserOn_.push_back(state);
            tooltxyz_.push_back(txyz);
        }//end while
    }//end if
}//end getToolpath


//////////////////////////////////////////////////////
//		getHeatSourceParameters()			    //
//////////////////////////////////////////////////////
void
DomainManager::getHeatSourceParameters()
{
    ifstream file;
    file.open(meshObj_->heatSourceFileName_.c_str());
    string line;
    if (file.is_open())
    {
        while(getline(file, line))
        {
            istringstream lines(line);
            vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
            vector<double> parameterValues(coords.begin(), coords.end());

            heatSourceParameterValues.push_back(parameterValues);
        }//end while
    }//end if
}//end getToolpath