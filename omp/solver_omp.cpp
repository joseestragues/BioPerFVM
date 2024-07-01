#include<vector>
#include <immintrin.h>
#include<chrono>
#include<fstream>
#include "./microenvironment_omp.h"
#include "../utils.h"



void microenvironment_omp::diffusion_decay_3D_solver() {

    std::ofstream file(timing_csv, std::ios::app);
    file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto apply_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for( int k=0; k < z_size ; k++ )
    {
        for( int j=0; j < y_size ; j++ )
        {
            // Thomas solver, x-direction

            // remaining part of forward elimination, using pre-computed quantities 
            int n = voxel_index(0,j,k);
            (*p_density_vectors)[n] /= thomas_denomx[0]; 

            for( int i=1; i < x_size ; i++ )
            {
                n = voxel_index(i,j,k); 
                axpy( &(*p_density_vectors)[n] , thomas_constant1 , (*p_density_vectors)[n-i_jump] ); 
                (*p_density_vectors)[n] /= thomas_denomx[i]; 
            }

            for( int i = x_size-2 ; i >= 0 ; i-- )
            {
                n = voxel_index(i,j,k); 
                naxpy( &(*p_density_vectors)[n] , thomas_cx[i] , (*p_density_vectors)[n+i_jump] ); 
            }

        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    file << duration_us << ",";

    // y-diffusion 
    start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for( int k=0; k < z_size ; k++ )
    {
        for( int i=0; i < x_size ; i++ )
        {
        // Thomas solver, y-direction

        // remaining part of forward elimination, using pre-computed quantities 

        int n = voxel_index(i,0,k);
        (*p_density_vectors)[n] /= thomas_denomy[0]; 

        for( int j=1; j < y_size ; j++ )
        {
        n = voxel_index(i,j,k); 
        axpy( &(*p_density_vectors)[n] , thomas_constant1 , (*p_density_vectors)[n-j_jump] ); 
        (*p_density_vectors)[n] /= thomas_denomy[j]; 
        }

        // back substitution 
        // n = voxel_index( mesh.x_coordinates.size()-2 ,j,k); 

        for( int j = y_size -2 ; j >= 0 ; j-- )
        {
        n = voxel_index(i,j,k); 
        naxpy( &(*p_density_vectors)[n] , thomas_cy[j] , (*p_density_vectors)[n+j_jump] ); 
        }

        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    file << duration_us << ",";

    // z-diffusion 

    start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for( int j=0; j < y_size ; j++ )
    {
        for( int i=0; i < x_size ; i++ )
        {
            // Thomas solver, y-direction

            // remaining part of forward elimination, using pre-computed quantities 

            int n = voxel_index(i,j,0);
            (*p_density_vectors)[n] /= thomas_denomz[0]; 

            // should be an empty loop if mesh.z_coordinates.size() < 2  
            for( int k=1; k < z_size ; k++ )
            {
            n = voxel_index(i,j,k); 
            axpy( &(*p_density_vectors)[n] , thomas_constant1 , (*p_density_vectors)[n-k_jump] ); 
            (*p_density_vectors)[n] /= thomas_denomz[k]; 
            }

            // back substitution 

            // should be an empty loop if mesh.z_coordinates.size() < 2 
            for( int k = z_size-2 ; k >= 0 ; k-- )
            {
            n = voxel_index(i,j,k); 
            naxpy( &(*p_density_vectors)[n] , thomas_cz[k] , (*p_density_vectors)[n+k_jump] ); 
            // n -= i_jump; 
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    file << duration_us << ",";

    start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    file << apply_us << std::endl;

    return; 
}