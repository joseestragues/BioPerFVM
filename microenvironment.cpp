#include<mpi.h>
#include <fstream>
#include "./microenvironment.h"
#include "./utils.h"

microenvironment::microenvironment(){
    x_size = 1;
    y_size = 1;
    z_size = 1;
    number_of_densities = 1;
    cube_side = 1.0;
    delta = 1.0;
    i_jump = 1;
    j_jump = 1;
    k_jump = 1;
    mpi_rank = 0;
    mpi_size = 0;
}

int microenvironment::p_index (int i, int j, int k) {
    return i*i_jump + j*j_jump + k*k_jump;
}

void microenvironment::init_densities(){

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j){
            int voxel = p_index(i,j,0);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
            voxel = p_index(i,j,z_size-1);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
        }
    }
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < z_size; ++k) {
        for (int j = 0; j < y_size; ++j){
            int voxel = p_index(0,j,k);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
            voxel = p_index(x_size -1,j,k);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
        }
    }
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < z_size; ++k) {
        for (int i = 0; i < x_size; ++i){
            int voxel = p_index(i,0,k);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
            voxel = p_index(i,y_size-1,k);
            for (int d = 0; d < number_of_densities; ++d)
                densities[voxel+d]=100;
        }
    }

}

int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// Function to compute the least common multiple (LCM)
int lcm(int a, int b) {
    return (a * b) / gcd(a, b);
}

void microenvironment::apply_dirichlet_conditions( void )
{ /*	#pragma omp parallel for 
	for( int i=0 ; i < dirichlet_indices.size() ; i++ )
	{ 
		density_vector( dirichlet_indices[i] ) = dirichlet_value_vectors[i]; }
    }*/
}
/*
void microenvironment::init_dirichlet(vector<double> dirichlet_values){
    dirichlet_value_vectors = dirichlet_values;
     //X boundaries
    if (mpi_rank == 0) {
        #pragma parallel for collapse(2)
        for (int j = 0 ; j < y_size; j++){
            for (int k = 0; k < z_size; k++){
                int index = j*j_jump +
                            k*k_jump;
                
            }
        }
    }

    if (mpi_rank == mpi_size -1) {
        #pragma parallel for collapse(2)
        for (int j = 0 ; j < y_size; j++){
            for (int k = 0; k < z_size; k++){
                int index = (x_size - 1)*i_jump + 
                        j*j_jump +
                        k*k_jump;
                
                
            }
        }
    }

    //Y boundaries
    #pragma parallel for collapse(2)
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
*/


void microenvironment::init_diffusion_coeficients(double cube_side, double delta, int d_num, int rank, int size, int vl, MPI_Comm mpi_world){

    mpi_comm = mpi_world;
    mpi_rank = rank;
    mpi_size = size;
    int voxels = cube_side/delta;
    MPI_Request send_req[mpi_size];
    MPI_Request recv_req[mpi_size];

    if (size == 0){
        std::cout << "Error: Number of MPI processes is equal to 0!" << std::endl;
        exit(1);
    }
    //Domain partition
    
    if (voxels%mpi_size != 0){
        if (mpi_rank == 0){
            std::cout << "Warning: X-dimension is not perfectly divisible between MPI processes!" << std::endl;
            std::cout << "  Microenvironment mesh not equally distributed" << std::endl;
        }
        int residual = voxels%mpi_size;
        if (mpi_rank < residual) x_size = (voxels/mpi_size) + 1; 
        else x_size = (voxels/mpi_size);
    }
    else{
        x_size =voxels/mpi_size;
    }
    
    y_size= voxels;
    z_size = voxels;
    number_of_densities = d_num;
    vector<double> zero(number_of_densities, 0.0);
    vector<double> one(number_of_densities, 1.0);
    double dt = 0.01;

    int step_size = (z_size * y_size) / granurality;

    snd_data_size = step_size * number_of_densities; // Number of data elements to be sent
    rcv_data_size = step_size * number_of_densities; // All p_density_vectors elements have same size, use anyone

    snd_data_size_last = ((z_size * y_size) % granurality) * number_of_densities; // Number of data elements to be sent
    rcv_data_size_last = ((z_size * y_size) % granurality) * number_of_densities;

    //Simulation
    densities.resize(x_size*y_size*z_size*number_of_densities, 0.0); 
    densities = vector<double> (x_size*y_size*z_size*number_of_densities);
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

    i_jump = y_size*z_size*number_of_densities;
    j_jump = z_size*number_of_densities;
    k_jump = number_of_densities; // M.thomas_j_jump * M.mesh.y_coordinates.size();

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

    /*--------------------------------------------------------------------*/
    /* In 1-D X decomposition, y and z-lines are contiguous and typically */
    /* the assignments below for y,z should not be changed                */
    /*--------------------------------------------------------------------*/

    thomas_cx.assign(x_size, thomas_constant1a);    // Fill b and c elements with -D * dt/dx^2
    thomas_denomx.assign(x_size, thomas_constant3); // Fill diagonal elements with (1 + 1/3 * lambda * dt + 2*D*dt/dx^2)

    if (rank == 0)
        thomas_denomx[0] = thomas_constant3a; // First diagonal element is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

    if (rank == (size - 1))
        thomas_denomx[x_size - 1] = thomas_constant3a; // Last diagonal element  is   (1 + 1/3 * lambda * dt + 1*D*dt/dx^2)

    if (rank == 0)
        if (x_size == 1) // This is an extreme case, won't exist, still if it does
        {                                     // then this must be at rank 0
            thomas_denomx[0] = one;
            thomas_denomx[0] += thomas_constant2;
        }
    if (rank == 0)
        thomas_cx[0] /= thomas_denomx[0]; // The first c element of tridiagonal matrix is div by first diagonal el.

    // axpy(1st, 2nd, 3rd) => 1st = 1st + 2nd * 3rd
    // the value at  size-1 is not actually used
    // Since value of size-1 is not used, it means it is the value after the last Diagonal element
    // cout << "Rank " << rank << endl;
    for (int ser_ctr = 0; ser_ctr <= size - 1; ser_ctr++)
    {
        if (rank == ser_ctr)
        {
            if (rank == 0 && rank <= size - 1) // If size=1, then this process does not send data
            {

                for (int i = 1; i <= x_size - 1; i++)
                {
                    axpy(&thomas_denomx[i], thomas_constant1, thomas_cx[i - 1]);
                    thomas_cx[i] /= thomas_denomx[i]; // the value at  size-1 is not actually used
                }
            }
            else
            {
                for (int i = 1; i <= x_size - 1; i++)
                {
                    axpy(&thomas_denomx[i], thomas_constant1, thomas_cx[i - 1]);
                    thomas_cx[i] /= thomas_denomx[i]; // the value at  size-1 is not actually used
                }
            }

            if (rank < (size - 1))
            {
                MPI_Isend(&(thomas_cx[x_size - 1][0]), thomas_cx[x_size - 1].size(), MPI_DOUBLE, ser_ctr + 1, 1111, mpi_comm, &send_req[0]);
            }
        }

        if (rank == (ser_ctr + 1) && (ser_ctr + 1) <= (size - 1))
        {

            std::vector<double> temp_cx(thomas_cx[0].size());

            MPI_Irecv(&temp_cx[0], temp_cx.size(), MPI_DOUBLE, ser_ctr, 1111, mpi_comm, &recv_req[0]);
            MPI_Wait(&recv_req[0], MPI_STATUS_IGNORE);

            axpy(&thomas_denomx[0], thomas_constant1, temp_cx); // CHECK IF &temp_cz[0] is OK, axpy() in BioFVM_vector.cpp
            thomas_cx[0] /= thomas_denomx[0];                   // the value at  size-1 is not actually used
        }

        MPI_Barrier(mpi_comm);
    }
    
    /*--------------------------------------------------------------------*/
    /* In 1-D X decomposition, z and y-lines are contiguous adn typically */
    /* the assignments below for z,y should not be changed                */
    /* Both the first voxel i.e. index 0 and last voxel i.e. index=       */
    /* y_coordinates.size()-1 are on the same process                     */
    /*--------------------------------------------------------------------*/

    thomas_cy.assign(y_size, thomas_constant1a);
    thomas_denomy.assign(y_size, thomas_constant3);
    thomas_denomy[0] = thomas_constant3a;
    thomas_denomy[y_size - 1] = thomas_constant3a;
    if (y_size == 1)
    {
        thomas_denomy[0] = one;
        thomas_denomy[0] += thomas_constant2;
    }
    thomas_cy[0] /= thomas_denomy[0];
    for (int i = 1; i <= y_size - 1; i++)
    {
        axpy(&thomas_denomy[i], thomas_constant1, thomas_cy[i - 1]);
        thomas_cy[i] /= thomas_denomy[i]; // the value at  size-1 is not actually used
    }

    thomas_cz.assign(z_size, thomas_constant1a);
    thomas_denomz.assign(z_size, thomas_constant3);
    thomas_denomz[0] = thomas_constant3a;
    thomas_denomz[z_size - 1] = thomas_constant3a;
    if (z_size == 1)
    {
        thomas_denomz[0] = one;
        thomas_denomz[0] += thomas_constant2;
    }
    thomas_cz[0] /= thomas_denomz[0];
    for (int i = 1; i <= z_size - 1; i++)
    {
        axpy(&thomas_denomz[i], thomas_constant1, thomas_cz[i - 1]);
        thomas_cz[i] /= thomas_denomz[i]; // the value at  size-1 is not actually used
    }
    //if (rank == 0) file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;

    /*
    string path  = "./timing/voxels_" + std::to_string((int)cube_side/2.0) + "/substrates_" + std::to_string(number_of_densities) 
    + "/factor_" + std::to_string(factor) +  "/" + std::to_string(mpi_size) + "_node.csv";*/
    std::ofstream file(timing_csv, std::ios::app);
    if (mpi_rank == 0) {
        file << "X-diffusion,Y-diffusion,Z-diffusion,Apply Dirichlet" << std::endl;
    } 
}

void microenvironment::init_vec_coeficients( int vl) {
    //Vectorization initialization
    gvec_size = lcm(number_of_densities, vl);
    
    //X-diffusion
    gthomas_constant1.resize(gvec_size, 0.0);
    auto dest_iter =  gthomas_constant1.begin();
    for (int j = 0; j < gvec_size; j+=number_of_densities){
        copy(thomas_constant1.begin(), thomas_constant1.end(), dest_iter);
        dest_iter+=number_of_densities;
    }

    gthomas_denomx.resize(x_size);
    gthomas_cx.resize(x_size);
    for (int i = 0; i < x_size; ++i){
        gthomas_denomx[i].resize(gvec_size, 0.0);
        gthomas_cx[i].resize(gvec_size, 0.0);
        auto dest_denomx = gthomas_denomx[i].begin();
        auto dest_cx = gthomas_cx[i].begin();
        for (int d = 0; d < gvec_size; d+=number_of_densities){
            copy(thomas_denomx[i].begin(), thomas_denomx[i].end(), dest_denomx);
            copy(thomas_cx[i].begin(), thomas_cx[i].end(), dest_cx);
            dest_denomx+=number_of_densities;
            dest_cx+=number_of_densities;
        }
    }
    //Y-diffusion
   
    gthomas_denomy.resize(y_size);
    gthomas_cy.resize(y_size);
    for (int j = 0; j < y_size; ++j){
        gthomas_denomy[j].resize(gvec_size, 0.0);
        gthomas_cy[j].resize(gvec_size, 0.0);
        auto dest_denomy = gthomas_denomy[j].begin();
        auto dest_cy = gthomas_cy[j].begin();
        for (int d = 0; d < gvec_size; d+=number_of_densities){
            copy(thomas_denomy[j].begin(), thomas_denomy[j].end(), dest_denomy);
            copy(thomas_cy[j].begin(), thomas_cy[j].end(), dest_cy);
            dest_denomy+=number_of_densities;
            dest_cy+=number_of_densities;
        }
    }
    //Z - diffusion

    gthomas_denomz.resize(z_size);
    gthomas_cz.resize(z_size);
    for (int k = 0; k < z_size; ++k){
        gthomas_denomz[k].resize(gvec_size, 0.0);
        gthomas_cz[k].resize(gvec_size, 0.0);
        auto dest_denomz = gthomas_denomz[k].begin();
        auto dest_cz = gthomas_cz[k].begin();
        for (int d = 0; d < gvec_size; d+=number_of_densities){
            copy(thomas_denomz[k].begin(), thomas_denomz[k].end(), dest_denomz);
            copy(thomas_cz[k].begin(), thomas_cz[k].end(), dest_cz);
            dest_denomz+=number_of_densities;
            dest_cz+=number_of_densities;
        }
    }
    /*
    string path  = "./timing/voxels_" + std::to_string((int)cube_side/2.0) + "/substrates_" + std::to_string(number_of_densities) 
    + "/factor_" + std::to_string(factor) +  "/" + std::to_string(mpi_size) + "_node.csv";*/
    
}
