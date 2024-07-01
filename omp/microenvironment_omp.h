#ifndef __microenvironment_omp_h__
#define __microenvironment_omp_h__
#include<iostream>
#include<vector>
#include<chrono>
using namespace std;
class microenvironment_omp {
    public:
        std::vector< std::vector<double> > aux; 
        std::vector< std::vector<double> >* p_density_vectors; 
        int  x_size, y_size, z_size,number_of_densities;
        double cube_side, delta;
        int i_jump, j_jump, k_jump;

        //Timing experiments variables
        std::string timing_csv; //Path to csv used to store timing results

        //Substrates values
        vector<double> diffusion_coefficients;
        vector<double> decay_rates;
        vector<bool> dirichlet_activation_vector;
        vector<double> dirichlet_value_vectors;
        //Tridiagonal variables
        vector<double> thomas_constant1, thomas_constant1a, thomas_constant2, thomas_constant3, thomas_constant3a;
        vector<vector<double>> thomas_denomx,thomas_denomy, thomas_denomz;
        vector<vector<double>> thomas_cx, thomas_cy, thomas_cz;
        //Dirichlet conditions
        vector<int> dirichlet_indices;

        void init_densities();

        int voxel_index(int i, int j, int k);

        void init_diffusion_coeficients(double cube_side, double delta, int d_num);

        void set_dirichlets();

        void apply_dirichlet_conditions(); 

        void diffusion_decay_3D_solver();

        void print_voxels_densities(std::string *file_name);

        microenvironment_omp();
};
#endif