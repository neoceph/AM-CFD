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
#include <stdio.h>
#include <ctype.h>
#include "Mesh.h"
#include <iomanip>
#include <chrono>
#include <ctime>

//////////////////////////////////////////////////////
//		Constructor			    //
//////////////////////////////////////////////////////
Mesh::Mesh(string fileName) 
{
    file.open(fileName.c_str());
}

//////////////////////////////////////////////////////
//		destructor			    //
//////////////////////////////////////////////////////
Mesh::~Mesh() 
{
    file.close();
}

//////////////////////////////////////////////////////
//		getDomainInfo 			    //
//////////////////////////////////////////////////////
void Mesh::getDomainInfo()
{
    // define some local parameters
    string line;
    int option = -1;

    // loop over the input file, line-by-line
    if(file.is_open()) 
    {
        while(getline(file, line)) 
        {
            if(line.at(0) == '$')
                continue;

            if(line == "*FLUID_FLOW_OFF")
            {
                navierstokes_ = false;
                continue;
            }
            else if(line == "*CALCULATE_SOLIDIFICATION_PARAMETERS")
            {
                solidificationparams_ = true;
                continue;
            }
            else if(line == "*SINGLE_TRACK_TOOL")
            {
                option = 3;
                continue;
            }
            else if(line == "*MATERIAL_PROPERTIES")
            {
                option = 4;
                continue;
            }
            else if(line == "*UNDER_RELAXATION")
            {
                option = 5;
                continue;
            }
            else if(line == "*BOUNDARY_TEMPERATURES")
            {
                option = 6;
                continue;
            }
            else if(line == "*MESH_FILE")
            {
                option = 7;
                continue;
            }
            else if(line == "*KEYWORD_ID")
            {
                option = 8;
                continue;
            }
            else if(line == "*PARAMETERS")
            {
                option = 9;
                continue;
            }
            else if(line == "*BOUNDARY_ENERGY_EQUATION")
            {
                option = 10;
                continue;
            }
            else if(line == "*CONTROL_TIMESTEP")
            {
                option = 11;
                continue;
            }
            else if(line == "*CONTROL_TERMINATION")
            {
                option = 12;
                continue;
            }
            else if(line == "*GAUSS_LASER_PARAMETERS")
            {
                option = 13;
                gaussheatsource_ = true; 
                continue;
            }
            else if(line == "*BOUNDARY_MOMENTUM_EQUATION")
            {
                option = 14;
                continue;
            }
            else if(line == "*VECTOR_OUT")
            {
                option = 15;
                continue;
            }
            else if(line == "*SCALAR_OUT")
            {
                option = 16;
                continue;
            }
            else if(line == "*CELL_OUT")
            {
                option = 17;
                continue;
            }
            else if(line == "*DIRECTORY")
            {
                option = 18;
                continue;
            }
            else if(line == "*OUTPUT_TEC")
            {
                option = 19;
                istec_ = true;
                continue;
            }
            else if(line == "*BOUNDARY_EVAPORATION")
            {
                option = 20;
                isevap_ = true;
                continue;
            }
            else if(line == "*SOLVE_SPECIES")
            {
                option = 21;
                species_ = true;
                continue;
            }
            else if(line == "*TOOLPATH_FILE")
            {
                option = 22;
                udtoolpath_ = true;
                continue;
            }
            else if(line == "*ACTIVATE_MOVING_MESH")
            {
                option = 23;
                movingmesh_ = true;
                continue;
            }
            else if(line == "*TRACK_ENERGY_BASED_FREE_SURFACE")
            {
                option = 24;
                energyfreesurface_ = true;
                continue;
            }
            else if(line == "*SELECT_DISCRETIZATION")
            {
                option = 25;
                continue;
            }
            else if(line == "*POWDER_BED_MODE")
            {
                option = 26;
                powderbed_ = true;
                continue;
            }
            else if(line == "*VOLUMETRIC_HEAT_SOURCE")
            {
                option = 27;
                volheatsource_ = true;
                continue;
            }
            else if(line == "*OUTPUT_FINAL_GEOMETRY")
            {
                option = 28;
                outputgeom_ = true;
                continue;
            }
            else if(line == "*OUTPUT_SECTION")
            {
                option = 29;
                outputsection_ = true;
                continue;
            }
            else if(line.front() == '*') 
            {
                option = -1;
                continue;
            }//end if

            // read in the vaules
            if(option == 3)	// *SINGLE_TRACK_TOOL
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    proc_params_.push_back(coords[ii]);
            }// option 3
            else if(option == 4)	// *MATERIAL_PROPERTIES
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    mat_props_.push_back(coords[ii]);
            }// option 4
            else if(option == 5)	// *UNDER_RELAXATION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    num_relax_.push_back(coords[ii]);
            }// option 5
            else if(option == 6)	// *BOUNDARY_TEMPERATURES
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    bound_cond_temp_.push_back(coords[ii]);
            }// option 6
            else if(option == 7)	// *MESH_FILE
            {
	        istringstream lines(line);
	        vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                meshFileName_ = coords[0];
                foundMesh_ = true;
            }// option 7
            else if (option == 8)	// *KEYWORD_ID
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                outFileName_ = coords[0];
            }// option 8
            else if (option == 9)	// *PARAMETERS
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                string paramName = coords[0]; 
                double value = stod(coords[1]);
                paramValues_[paramName] = value;
            }// option 9
            else if(option == 10)	// *BOUNDARY_ENERGY_EQUATION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    bound_cond_energy_.push_back(coords[ii]);
            }// option 10
            else if(option == 11)	// *CONTROL_TIMESTEP
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                delt_ = coords[0];
                outtime_ = coords[1];
            }// option 11
            else if(option == 12)	// *CONTROL_TERMINATION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                finaltime_ = coords[0];
                maxit_ = coords[1];
                nonlintol_ = coords[2];
            }// option 12
            else if(option == 13)	// *GAUSS_LASER_PARAMETERS
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    gauss_.push_back(coords[ii]);
            }// option 13
            else if(option == 14)	// *BOUNDARY_MOMENTUM_EQUATION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    bound_cond_momentum_.push_back(coords[ii]);
            }// option 14
            else if (option == 15)	// *VECTOR_OUT
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                outputVectorNamesNODE_.push_back(coords[0]);
            }// option 15
            else if (option == 16)	// *SCALAR_OUT
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                outputScalarNamesNODE_.push_back(coords[0]);
            }// option 16
            else if (option == 17)	// *CELL_OUT
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                outputScalarNamesCELL_.push_back(coords[0]);
            }// option 17
            else if (option == 18)	// *DIRECTORY
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                directoryName_ = coords[0];
            }// option 18
            else if (option == 19)	// *OUTPUT_TEC
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                tectype_ = coords[0];
            }// option 19
            else if (option == 20)	// *BOUNDARY_EVAPORATION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    evap_params_.push_back(coords[ii]);
            }// option 20
            else if (option == 21)	// *SOLVE_SPECIES
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    species_params_.push_back(coords[ii]);
            }// option 21
            else if (option == 22)	// *TOOLPATH_FILE
            {
	        istringstream lines(line);
	        vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                toolFileName_ = coords[0]; 
            }// option 22
            else if(option == 23)	// *ACTIVATE_MOVING_MESH
            {
                istringstream lines(line);
                vector<int> coords((istream_iterator<int>(lines)), istream_iterator<int>());
                xmoveindex_ = coords[0]-1;
                ymoveindex_ = coords[1]-1;
            }// option 23
            else if(option == 24)	// *TRACK_ENERGY_BASED_FREE_SURFACE
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    energyfree_params_.push_back(coords[ii]);
            }// option 24
            else if (option == 25)	// *SELECT_DISCRETIZATION
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                if (coords[0] == "power_law")
                {
                    powerlaw_ = true;
                    upwind_ = false;
                }
                else if (coords[0] == "exponential_law")
                {
                    exponential_ = true;
                    upwind_ = false;
                }
            }// option 25
            else if(option == 26) 	// *POWDER_BED_MODE
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    powderbed_params_.push_back(coords[ii]);
            }// option 26
            else if(option == 27) 	// *VOLUMETRIC_HEAT_SOURCE
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    volheatsource_params_.push_back(coords[ii]);
            }// option 27
            else if (option == 28)	// *OUTPUT_FINAL_GEOMETRY
            {
                istringstream lines(line);
                vector<string> coords((istream_iterator<string>(lines)), istream_iterator<string>());
                geomOutName_ = coords[0];
            }// option 28
            else if(option == 29)	// *OUTPUT_SECTION
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    outputsection_params_.push_back(coords[ii]);
            }// option 29
        }//end option conditional loops
    }//end if 
}//end getDomainInfo

//////////////////////////////////////////////////////
//		assignParameters		    //
//////////////////////////////////////////////////////
void 
Mesh::assignParameters()
{
    // Assign parameter values
    for (map<string,double>::iterator it = paramValues_.begin();
         it != paramValues_.end(); it++)
    {
        string paramName = it->first;
        double value = it->second;
        if (paramName == "accelGravity")
        {
            grav_ = value;
        }
        else if (paramName == "radBoltzmann")
        {
            sigm_ = value;
        }
        else if (paramName == "ambientTemp")
        {
            tempamb_ = value;
        }
        else if (paramName == "ambientPres")
        {
            pressamb_ = value;
        }
        else if (paramName == "gasConstant")
        {
            gasconst_ = value;
        }
    }//end for(it)
}//assignParameters

//////////////////////////////////////////////////////
//		readMeshFile    		    //
//////////////////////////////////////////////////////
void 
Mesh::readMeshFile()
{
    // define some local variables
    ifstream meshFile;
    meshFile.open(meshFileName_);
    string line;
    int option = -1;

    // make sure the correct mesh was given
    ifstream meshCheck(meshFileName_);
    if (!meshCheck)
    {
        cout << "****************************************************************\n";
        cout << " ERROR: SPECIFICED MESH WAS NOT FOUND, SEARCHED IN:\n" ;
        cout << meshFileName_ << "\n";
        cout << "****************************************************************\n";
        exit(EXIT_FAILURE);
    }//end if

    // loop throgh the mesh file to get grid informtion
    if(meshFile.is_open()) 
    {
        while(getline(meshFile, line)) 
        {
            if(line.at(0) == '$')
                continue;

            if(line == "*X_GRID_PARAMETERS")
            {
                option = 1;
                continue;
            }
            else if(line == "*Y_GRID_PARAMETERS")
            {
                option = 2;
                continue;
            }
            else if(line == "*Z_GRID_PARAMETERS")
            {
                option = 3;
                continue;
            }
            else if(line.at(0) == '*') 
            {
              option = -1;
              continue;
            }

            if(option == 1)	// *X_GRID_PARAMETERS
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    x_grid_.push_back(coords[ii]);
            }// option 1
            else if(option == 2)	// *Y_GRID_PARAMETERS
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    y_grid_.push_back(coords[ii]);
            }// option 2
            else if(option == 3)	// *Z_GRID_PARAMETERS
            {
                istringstream lines(line);
                vector<double> coords((istream_iterator<double>(lines)), istream_iterator<double>());
                for (int ii = 0; ii < coords.size(); ++ii)
                    z_grid_.push_back(coords[ii]);
            }// option 3
        }//end while
    }//end if

    // close the mesh file
    meshFile.close();
}//readMeshFile

//////////////////////////////////////////////////////
//		outputInputs		    	    //
//////////////////////////////////////////////////////
void 
Mesh::outputInputs(char* &inputName)
{
    // define some local variables
    const std::string SecondEleWidth  = "  ";
    const std::string ThirdEleWidth   = "    ";
    auto begin = std::chrono::high_resolution_clock::now();

    // close the input file
    file.close();
 
    // title output
    cout << "|==================================================================|\n";
    cout << "|       Northwestern University AM-CFD Software                    |\n";
    cout << "|                                                                  |\n";
    cout << "|                        Developer:                                |\n";
    cout << "|                      Abdullah A Amin                             |\n";
    cout << "|==================================================================|\n";
    cout << "\n\n";

    // read in mesh file
    if(foundMesh_) 
    {
        readMeshFile();
    }
    else 
    {
      cout << "****************************************************************\n";
      cout << " ERROR: Need to specify mesh file name in: "<< inputName << " with keyword \"*MESH_FILE\" \n";
      cout << "****************************************************************\n";
      exit(EXIT_FAILURE);
    }//end if

    // output user given information about the simulation
    std::time_t starttime = std::chrono::system_clock::to_time_t(begin);
    cout << "\tStarted simulation on " << std::ctime(&starttime) << "\n";

    cout << "===================================================================\n";
    cout << "\tUser Input Information\n\n";
    cout << SecondEleWidth << "Reading in Input File: " << inputName << "\n";
    cout << SecondEleWidth << "Reading in Mesh File: " << meshFileName_ << "\n";
    if(udtoolpath_)
    {
        cout << SecondEleWidth << "Reading in Toolpath File: " << toolFileName_ << "\n";
    }
    cout << SecondEleWidth << "Writing out to output database name: " << outFileName_ << "\n";
    cout << SecondEleWidth << "Input parameter list: " << "\n";
    for (map<string,double>::iterator it = paramValues_.begin();
         it != paramValues_.end(); it++)
    {
        string paramName = it->first;
        double value = it->second;
        cout << ThirdEleWidth << paramName << ":  \t  " << value << "\n";
    }//end for(it)
    cout << SecondEleWidth << "User output scalar nodal variables: \n";
    for (vector<string>::iterator it = outputScalarNamesNODE_.begin();
         it != outputScalarNamesNODE_.end(); it++)
    {
      cout << ThirdEleWidth << *it << "\n";
    }//end for(it)
    cout << SecondEleWidth << "User output vector nodal variables: \n";
    for (vector<string>::iterator it = outputVectorNamesNODE_.begin();
         it != outputVectorNamesNODE_.end(); it++)
    {
      cout << ThirdEleWidth << *it << "\n";
    }//end for(it)
    cout << SecondEleWidth << "User output scalar cell variables: \n";
    for (vector<string>::iterator it = outputScalarNamesCELL_.begin();
         it != outputScalarNamesCELL_.end(); it++)
    {
      cout << ThirdEleWidth << *it << "\n";
    }//end for(it)

    cout << "===================================================================\n\n";
    cout << endl;
}//end outputInputs
