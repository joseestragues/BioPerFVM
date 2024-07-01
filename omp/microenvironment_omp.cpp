#include <fstream>
#include "./microenvironment_omp.h"
#include "../utils.h"

microenvironment_omp::microenvironment_omp(){
    x_size = 1;
    y_size = 1;
    z_size = 1;
    number_of_densities = 1;
    cube_side = 1.0;
    delta = 1.0;
    i_jump = 1;
    j_jump = 1;
    k_jump = 1;
}

int microenvironment_omp::voxel_index(int i, int j, int k) {
    return i*i_jump + j*j_jump + k*k_jump;
}

void microenvironment_omp::init_densities(){

    set_dirichlets();

    vector<double> densities(number_of_densities, 100.0);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j){
            int voxel = voxel_index(i,j,0);
            (*p_density_vectors)[voxel]=densities;
            voxel = voxel_index(i,j,z_size-1);
            (*p_density_vectors)[voxel]=densities;
        }
    }
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < z_size; ++k) {
        for (int j = 0; j < y_size; ++j){
           int voxel = voxel_index(0,j,k);
            (*p_density_vectors)[voxel]=densities;
            voxel = voxel_index(x_size-1,j,k);
            (*p_density_vectors)[voxel]=densities;
        }
    }
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < z_size; ++k) {
        for (int i = 0; i < x_size; ++i){
            int voxel = voxel_index(i,0,k);
            (*p_density_vectors)[voxel]=densities;
            voxel = voxel_index(i,y_size-1,k);
            (*p_density_vectors)[voxel]=densities;
        }
    }

}

void microenvironment_omp::init_diffusion_coeficients(double cube_side, double delta, int d_num){

    int voxels = cube_side/delta;

    x_size =voxels;
    y_size = voxels;
    z_size = voxels;
    i_jump = 1;
    j_jump = x_size;
    k_jump = x_size*y_size;
    number_of_densities = d_num;
    vector<double> zero(number_of_densities, 0.0);
    vector<double> one(number_of_densities, 1.0);
    double dt = 0.01;

    //Simulation
    aux.resize(x_size*y_size*z_size, zero);
    
    p_density_vectors = &aux;
    //*p_density_vectors = vector<vector<double>> (x_size*y_size*z_size, zero);
    /*(*p_density_vectors).resize(x_size*y_size*z_size); 
    for (int i = 0; i < x_size*y_size*z_size; ++i) {
        (*p_density_vectors)[i].resize(number_of_densities, 0.0);
    }*/
        
    diffusion_coefficients.resize(number_of_densities);
    decay_rates.resize(number_of_densities);
    dirichlet_activation_vector.resize(number_of_densities);
    dirichlet_value_vectors.resize(number_of_densities);

    for (int i = 0; i < number_of_densities; i++)
    {
        diffusion_coefficients[i] = 75000.0 + (i*1000.0);
        decay_rates[i] = 0.05 + (i*0.005);
        dirichlet_activation_vector[i] = true;
        dirichlet_value_vectors[i] = 100.0;
    }

    //Thomas initialization
    thomas_denomx.resize(x_size, zero); // sizeof(x_coordinates) = local_x_nodes, denomx is the main diagonal elements
    thomas_cx.resize(x_size, zero);     // Both b and c of tridiagonal matrix are equal, hence just one array needed

    /*-------------------------------------------------------------*/
    /* y_coordinates are of size of local_y_nodes.                 */
    /* Each line of Voxels going                                   */
    /* from bottom to top forms a tridiagonal system of Equations  */
    /*-------------------------------------------------------------*/

    thomas_denomy.resize(y_size, zero);
    thomas_cy.resize(y_size, zero);

    /*-------------------------------------------------------------*/
    /* z_coordinates are of size of local_z_nodes.                 */
    /* Each line of Voxels going                                   */
    /* from front to back forms a tridiagonal system of Equations  */
    /*-------------------------------------------------------------*/

    thomas_denomz.resize(z_size, zero);
    thomas_cz.resize(z_size, zero);

    /*-------------------------------------------------------------*/
    /* For X-decomposition thomas_i_jump - 1 can be in the previous*/
    /* process and thomas_i_jump+1 can be in the next processs     */
    /* hence we can use thomas_j_jump and thomas_k_jump safely     */
    /* but we CANNOT use thomas_i_jump safely                      */
    /*-------------------------------------------------------------*/

    //i_jump = y_size*z_size*number_of_densities;
    //j_jump = z_size*number_of_densities;
    //k_jump = number_of_densities; // M.thomas_j_jump * M.mesh.y_coordinates.size();

    /*-------------------------------------------------------------*/
    /* This part below of defining constants SHOULD typically      */
    /* not change during parallelization.                          */
    /*-------------------------------------------------------------*/

    thomas_constant1 = diffusion_coefficients; // dt*D/dx^2
    thomas_constant1a = zero;                  // -dt*D/dx^2;
    thomas_constant2 = decay_rates;            // (1/3)* dt*lambda
    thomas_constant3 = one;                    // 1 + 2*constant1 + constant2;
    thomas_constant3a = one;                   // 1 + constant1 + constant2;

    thomas_constant1 *= dt;
    thomas_constant1 /= delta; //dx
    thomas_constant1 /= delta; //dx

    thomas_constant1a = thomas_constant1;
    thomas_constant1a *= -1.0;

    thomas_constant2 *= dt;
    thomas_constant2 /= 3.0; // for the LOD splitting of the source, division by 3 is for 3-D

    thomas_constant3 += thomas_constant1;
    thomas_constant3 += thomas_constant1;
    thomas_constant3 += thomas_constant2;

    thomas_constant3a += thomas_constant1;
    thomas_constant3a += thomas_constant2;

    // Thomas solver coefficients 

    thomas_cx.assign( x_size , thomas_constant1a ); 
    thomas_denomx.assign( x_size  , thomas_constant3 ); 
    thomas_denomx[0] = thomas_constant3a; 
    thomas_denomx[ x_size-1 ] = thomas_constant3a; 
    if( x_size == 1 )
    { thomas_denomx[0] = one; thomas_denomx[0] += thomas_constant2; } 

    thomas_cx[0] /= thomas_denomx[0]; 
    for( int i=1 ; i <= x_size-1 ; i++ )
    { 
        axpy( &thomas_denomx[i] , thomas_constant1 , thomas_cx[i-1] ); 
        thomas_cx[i] /= thomas_denomx[i]; // the value at  size-1 is not actually used  
    }

    thomas_cy.assign( y_size , thomas_constant1a ); 
    thomas_denomy.assign( y_size  , thomas_constant3 ); 
    thomas_denomy[0] = thomas_constant3a; 
    thomas_denomy[ y_size-1 ] = thomas_constant3a; 
    if( y_size == 1 )
    { thomas_denomy[0] = one; thomas_denomy[0] += thomas_constant2; } 

    thomas_cy[0] /= thomas_denomy[0]; 
    for( int i=1 ; i <= y_size-1 ; i++ )
    { 
        axpy( &thomas_denomy[i] , thomas_constant1 , thomas_cy[i-1] ); 
        thomas_cy[i] /= thomas_denomy[i]; // the value at  size-1 is not actually used  
    }

    thomas_cz.assign( z_size , thomas_constant1a ); 
    thomas_denomz.assign( z_size  , thomas_constant3 ); 
    thomas_denomz[0] = thomas_constant3a; 
    thomas_denomz[ z_size-1 ] = thomas_constant3a; 
    if( z_size == 1 )
    { thomas_denomz[0] = one; thomas_denomz[0] += thomas_constant2; } 

    thomas_cz[0] /= thomas_denomz[0]; 
    for( int i=1 ; i <= z_size-1 ; i++ )
    { 
        axpy( &thomas_denomz[i] , thomas_constant1 , thomas_cz[i-1] ); 
        thomas_cz[i] /= thomas_denomz[i]; // the value at  size-1 is not actually used  
    }	
    
    std::ofstream file(timing_csv, std::ios::app);
    file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
     
}

void microenvironment_omp::print_voxels_densities(std::string *file_name)
{
    // std::string filename = std::to_string(rank) + "_" + *file_name;
    //std::cout << "Rank " << rank << " esta imprimiendo densidades" << std::endl;
    std::string filename = *file_name;
    std::ofstream outputFile(filename, std::ios::app);
    
    for (int i = 0; i < x_size; i++)
    {
        for (int j = 0; j < y_size; j++)
        {
            for (int k = 0; k < z_size; k++)
            {
                int voxel = voxel_index(i,j,k);
                outputFile << + i << " " << j << " " << k << " : ";
                for (int d = 0; d < number_of_densities; ++d)
                {
                    outputFile << (*p_density_vectors)[voxel][d] << " ";
                    //std::cout << densities[index] << " ";
                }
                //std::cout << endl;
                outputFile << std::endl;
            }
        }
    }
}
