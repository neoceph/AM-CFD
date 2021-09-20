#ifndef MESH_H_
#define MESH_H_

// standard headers
#include <vector>
#include <fstream>
#include <map>

using namespace std;

class Mesh 
{
    public:
        // Constructor
        Mesh(string fileName);
         
        // Destructor
        virtual ~Mesh();
  
        // Local variables
          // storage for the input file values
        vector<double> proc_params_;
        vector<double> mat_props_;
        vector<double> num_relax_;
        vector<double> bound_cond_temp_;
        vector<double> bound_cond_energy_;
        vector<double> bound_cond_momentum_;
        vector<double> x_grid_;
        vector<double> y_grid_;
        vector<double> z_grid_;
        vector<double> gauss_;
        vector<double> evap_params_;
        vector<double> species_params_;
        vector<double> energyfree_params_;
        vector<double> powderbed_params_;
        vector<double> volheatsource_params_;
        vector<double> outputsection_params_;
        vector<string>  outputScalarNamesNODE_, outputScalarNamesCELL_, outputVectorNamesNODE_;
        double finaltime_, maxit_, nonlintol_, delt_, outtime_;
        std::string toolFileName_;
        std::string outFileName_;
        std::string meshFileName_;
        std::string directoryName_ = "Results";
        std::string tectype_;
        std::string geomOutName_;
        bool istec_ = false;
        bool navierstokes_ = true;
        bool species_ = false;
        bool isevap_ = false;
        bool udtoolpath_ = false;
        bool solidificationparams_ = false;
        bool movingmesh_ = false;
        bool energyfreesurface_ = false;
        int xmoveindex_, ymoveindex_;
        bool upwind_ = true;
        bool powerlaw_ = false;
        bool exponential_= false;
        bool foundMesh_ = false;
        bool powderbed_ = false;
        bool volheatsource_ = false;
        bool gaussheatsource_ = false;
        bool outputgeom_ = false;
        bool outputsection_ = false;
          // parameters from input file
        map<string, double> paramValues_;
        double grav_, sigm_, tempamb_, pressamb_, gasconst_; 
  
        // Local functions/subroutines
        void getDomainInfo();

        void openNew(string fileName);

        void assignParameters();

        void readMeshFile();

        void outputInputs(char* &inputName);

    private:
        // Private variables
        ifstream file;

        // Private functions/subroutines
        string getParam();
};

#endif /* MESH_H_ */
