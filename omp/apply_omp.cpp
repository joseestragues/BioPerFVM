#include<vector>
#include "./microenvironment_omp.h"

void microenvironment_omp::apply_dirichlet_conditions( void )
{
	#pragma omp parallel for 
	for( int i=0 ; i < dirichlet_indices.size() ; i++ )
	{ 
		(*p_density_vectors)[dirichlet_indices[i]] = dirichlet_value_vectors; 
    }
}

void microenvironment_omp::set_dirichlets(){

    vector<double> densities(number_of_densities, 100.0);
    dirichlet_value_vectors = densities;
    
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j){
            int voxel = voxel_index(i,j,0);
            dirichlet_indices.push_back(voxel);
            voxel = voxel_index(i,j,z_size-1);
            dirichlet_indices.push_back(voxel);
        }
    }
    
    for (int k = 0; k < z_size; ++k) {
        for (int j = 0; j < y_size; ++j){
           int voxel = voxel_index(0,j,k);
           dirichlet_indices.push_back(voxel);
           voxel = voxel_index(x_size-1,j,k);
           dirichlet_indices.push_back(voxel);
        }
    }
    
    for (int k = 0; k < z_size; ++k) {
        for (int i = 0; i < x_size; ++i){
            int voxel = voxel_index(i,0,k);
            dirichlet_indices.push_back(voxel);
            voxel = voxel_index(i,y_size-1,k);
            dirichlet_indices.push_back(voxel);
        }
    }
}
