#include<vector>
#include <immintrin.h>
#include<chrono>
#include<fstream>
#include "./microenvironment.h"



void microenvironment::apply_dirichlet_conditions_boundaries()
{
    #pragma omp parallel
    {
    //X boundaries
    if (mpi_rank == 0) {
        #pragma omp for collapse(2)
        for (int j = 0 ; j < y_size; j++){
            for (int k = 0; k < z_size; k++){
                int index = j*j_jump +
                            k*k_jump;
                for (int d = 0; d < number_of_densities; d++)
                {
                    if (dirichlet_activation_vector[d] == true)
                    {
                        densities[index+d] = dirichlet_value_vectors[d];
                    }
                }
            }
        }
    }

    if (mpi_rank == mpi_size -1) {
        #pragma omp for collapse(2)
        for (int j = 0 ; j < y_size; j++){
            for (int k = 0; k < z_size; k++){
                int index = (x_size - 1)*i_jump + 
                        j*j_jump +
                        k*k_jump;
                for (int d = 0; d < number_of_densities; d++)
                {
                    if (dirichlet_activation_vector[d] == true)
                    {
                        densities[index+d] = dirichlet_value_vectors[d];
                    }
                }
            }
        }
    }

    //Y boundaries
    #pragma omp for collapse(2)
    for (int i = 0 ; i < x_size; i++){
        for (int k = 0; k < z_size; k++){
            int index = i*i_jump +
                        k*k_jump;
            for (int d = 0; d < number_of_densities; d++)
            {
                if (dirichlet_activation_vector[d] == true)
                {
                    densities[index+d] = dirichlet_value_vectors[d];
                }
            }
            index = i*i_jump + 
                    (y_size -1)*j_jump +
                    k*k_jump;
            for (int d = 0; d < number_of_densities; d++)
            {
                if (dirichlet_activation_vector[d] == true)
                {
                    densities[index+d] = dirichlet_value_vectors[d];
                }
            }
        }
    }

    //Z boundaries
    #pragma omp for collapse(2)
    for (int i = 0 ; i < x_size; i++){
        for (int j = 0; j < y_size; j++){
            int index = i*i_jump +
                        j*j_jump;
            for (int d = 0; d < number_of_densities; d++)
            {
                if (dirichlet_activation_vector[d] == true)
                {
                    densities[index+d] = dirichlet_value_vectors[d];
                }
            }
            index = i*i_jump + 
                    j*j_jump +
                    (z_size-1)*k_jump;
            for (int d = 0; d < number_of_densities; d++)
            {
                if (dirichlet_activation_vector[d] == true)
                {
                    densities[index+d] = dirichlet_value_vectors[d];
                }
            }
        }
    }
    }
}
