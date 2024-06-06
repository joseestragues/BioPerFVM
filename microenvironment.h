#ifndef __microenvironment_h__
#define __microenvironment_h__
#include<iostream>
#include<vector>
#include<chrono>
#include<mpi.h>
using namespace std;
class microenvironment {
    public:
        MPI_Comm mpi_comm;
        vector<double> densities;
        int  x_size, y_size, z_size,number_of_densities;
        double cube_side, delta;
        int i_jump, j_jump, k_jump;

        //MPI parameters
        int mpi_rank, mpi_size;
        int granurality;
        int factor;
        int snd_data_size, snd_data_size_last;
        int rcv_data_size, rcv_data_size_last;
        
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
        //Tridiagonal variables for vectorization
        vector<double> gthomas_constant1, gthomas_constant1a, gthomas_constant2, gthomas_constant3, gthomas_constant3a;
        vector<vector<double>> gthomas_denomx,gthomas_denomy, gthomas_denomz;
        vector<vector<double>> gthomas_cx, gthomas_cy, gthomas_cz;
        //Dirichlet conditions
        vector<int> dirichlet_indices;
        int gvec_size;

        void init_densities();

        int p_index (int i, int j, int k);

        void init_diffusion_coeficients(double cube_side, double delta, int d_num, int rank, int size, int vl, MPI_Comm mpi_world);
        void init_diffvec_coefficients();

        void set_dirichlets();

        void apply_dirichlet_conditions(); 
        void apply_dirichlet_conditions_boundaries();

        void diffusion_decay_3D_solver();
        void diffusion_decay_3D_solver_256D();
        void diffusion_decay_3D_solver_512D(); 

        microenvironment();
};
#endif