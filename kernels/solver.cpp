#include<vector>
#include <immintrin.h>
#include<chrono>
#include<fstream>
#include "./microenvironment.h"


void microenvironment::diffusion_decay_3D_solver()
{
    MPI_Request send_req[granurality], recv_req[granurality];
    double block3d[i_jump]; //Aux structure of the size: Y*Z*Substrates
    /*
    string path  = "./timing/voxels_" + std::to_string((int)cube_side/2.0) + "/substrates_" + std::to_string(number_of_densities) 
    + "/factor_" + std::to_string(factor) +  "/" + std::to_string(mpi_size) + "_node.csv";*/
    std::ofstream file(timing_csv, std::ios::app);
    int vl = 4;
    auto start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions_boundaries();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto apply_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    //cout << "Peto aqui!" << endl;
    start_time = std::chrono::high_resolution_clock::now();
    //X-diffusion vector_x_v2
    //start_time = std::chrono::high_resolution_clock::now();
     if (mpi_rank == 0)
        {
            // cout << "Rank " << rank << "is computing" << endl;
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = step * snd_data_size;
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size; index += number_of_densities)
                {
                    int index_dec = index; 
                    for (int d = 0; d < number_of_densities; d++)
                    {
                        densities[index + d] /= thomas_denomx[0][d];
                    }

                    for (int i = 1; i < x_size; i++)
                    {
                        
                        int index_inc = index_dec + i_jump;
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index_inc + d] += thomas_constant1[d] * densities[index_dec + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index_inc + d] /= thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (mpi_size > 1) {
                    int x_end = x_size - 1;
                    int offset = step * snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&(densities[x_end * i_jump + offset]), snd_data_size, MPI_DOUBLE, mpi_rank + 1, step, MPI_COMM_WORLD, &send_req[step]);
                }
            }
            //Last iteration
            if (snd_data_size_last != 0) {
                int initial_index = granurality * snd_data_size;
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size_last; index += number_of_densities)
                {
                    int index_dec = index; 
                    for (int d = 0; d < number_of_densities; d++)
                    {
                        densities[index + d] /= thomas_denomx[0][d];
                    }

                    for (int i = 1; i < x_size; i++)
                    {
                        int index_inc = index_dec + i_jump;
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index_inc + d] += thomas_constant1[d] * densities[index_dec + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index_inc + d] /= thomas_denomx[i][d];
                        }
                        index_dec = index_inc;
                    }
                }

                if (mpi_size > 1) {
                    int x_end = x_size - 1;
                    int offset = granurality * snd_data_size;
                    MPI_Status status;
                    MPI_Isend(&(densities[x_end * i_jump + offset]), snd_data_size_last, MPI_DOUBLE, mpi_rank + 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                    
                }
            }
        }
        else
        {
            if (mpi_rank >= 1 && mpi_rank <= (mpi_size - 1))
            {
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * snd_data_size;
                    MPI_Irecv(&(block3d[initial_index]), rcv_data_size, MPI_DOUBLE, mpi_rank-1, step, MPI_COMM_WORLD, &(recv_req[step]));
                }
                if (snd_data_size_last != 0)
                    MPI_Irecv(&(block3d[granurality*snd_data_size]), rcv_data_size_last, MPI_DOUBLE, mpi_rank-1, granurality, MPI_COMM_WORLD, &(recv_req[granurality]));
                for (int step = 0; step < granurality; ++step)
                {
                    int initial_index = step * snd_data_size;
                    MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                    #pragma omp parallel for
                    for (int index = initial_index; index < initial_index + snd_data_size; index += number_of_densities)
                    {
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index + d] += thomas_constant1[d] * block3d[index + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[0];
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index + d] /= thomas_denomx[0][d];
                        }

                        for (int i = 1; i < x_size; i++)
                        {

                            int index_inc = index_dec + i_jump;
                            // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                            for (int d = 0; d < k_jump; d++)
                            {
                                densities[index_inc + d] += thomas_constant1[d] * densities[index_dec + d];
                            }
                            //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                            for (int d = 0; d < number_of_densities; d++)
                            {
                                densities[index_inc + d] /= thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    if (mpi_rank < (mpi_size - 1))
                    {
                        int x_end = x_size - 1;
                        MPI_Isend(&(densities[x_end * i_jump + initial_index]), snd_data_size, MPI_DOUBLE, mpi_rank + 1, step, MPI_COMM_WORLD, &send_req[step]);
                    }
                }
                if (snd_data_size_last != 0)
                {
                    int initial_index = granurality * snd_data_size;
                    MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE); //Need to change
                    #pragma omp parallel for
                    for (int index = initial_index; index < initial_index + snd_data_size_last; index += number_of_densities)
                    {
                        // axpy(&(*M.microenvironment)[n], M.thomas_constant1, block3d[k][j]);
                        int index_dec = index;
                        for (int d = 0; d < k_jump; d++)
                        {
                            densities[index + d] += thomas_constant1[d] * block3d[index + d];
                        }

                        //(*M.microenvironment)[n] /= M.thomas_denomx[0];
                        for (int d = 0; d < number_of_densities; d++)
                        {
                            densities[index + d] /= thomas_denomx[0][d];
                        }

                        for (int i = 1; i < x_size; i++)
                        {

                            int index_inc = index_dec + i_jump;
                            // axpy(&(*M.microenvironment)[n], M.thomas_constant1, (*M.microenvironment)[n - M.i_jump]);
                            for (int d = 0; d < number_of_densities; d++)
                            {
                                densities[index_inc + d] += thomas_constant1[d] * densities[index_dec + d];
                            }
                            //(*M.microenvironment)[n] /= M.thomas_denomx[i];
                            for (int d = 0; d < number_of_densities; d++)
                            {
                                densities[index_inc + d] /= thomas_denomx[i][d];
                            }

                            index_dec = index_inc;
                        }
                        
                    }
                    // End of computation region
                    if (mpi_rank < (mpi_size - 1))
                    {
                        int x_end = x_size - 1;
                        MPI_Request aux;
                        MPI_Isend(&(densities[x_end * i_jump + initial_index]), snd_data_size_last, MPI_DOUBLE, mpi_rank + 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                      
                    }
                }
            }
        }
        
        /*-----------------------------------------------------------------------------------*/
        /*                         CODE FOR BACK SUBSITUTION                                 */
        /*-----------------------------------------------------------------------------------*/
        
        //cout << "Rank " << rank << " starting backward substitution" << endl;

        if (mpi_rank == (mpi_size - 1))
        {
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = ((x_size - 1)*i_jump) + (step * snd_data_size);
                #pragma omp parallel for

                for (int index = initial_index; index < initial_index + snd_data_size; index += number_of_densities)
                {
                    int index_aux = index;
                    //int index = j * M.j_jump + k * M.k_jump + (M.mesh.x_coordinates.mpi_size() - 1) * M.i_jump;
                    for (int i = x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < k_jump; d++)
                        {
                            densities[index_dec + d] -= thomas_cx[i][d] * densities[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (mpi_size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&(densities[step * snd_data_size]), snd_data_size, MPI_DOUBLE, mpi_rank - 1, step, MPI_COMM_WORLD, &send_req[step]);
                }
            }

            //Last iteration
            if (snd_data_size_last != 0) {
                int initial_index = ((x_size - 1)*i_jump) + (granurality * snd_data_size);
                #pragma omp parallel for
                for (int index = initial_index; index < initial_index + snd_data_size_last; index += number_of_densities)
                {
                    int index_aux = index;
                    for (int i = x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < k_jump; d++)
                        {
                            densities[index_dec + d] -= thomas_cx[i][d] * densities[index_aux + d];
                        }
                        index_aux = index_dec;
                    }
                }
                if (mpi_size > 1) {
                    MPI_Request aux;
                    MPI_Isend(&(densities[granurality * snd_data_size]), snd_data_size_last, MPI_DOUBLE, mpi_rank - 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                    //cout << "Rank " << rank << " has send" << endl;
                }
            
            }
        }
        else
        {
            MPI_Request status[granurality]; 
            MPI_Request status_last;
            for (int step = 0; step < granurality; ++step)
            {
                MPI_Irecv(&(block3d[step*snd_data_size]), rcv_data_size, MPI_DOUBLE, mpi_rank+1, step, MPI_COMM_WORLD, &recv_req[step]);
            }
            if (snd_data_size_last != 0)
                MPI_Irecv(&(block3d[granurality*snd_data_size]), rcv_data_size_last, MPI_DOUBLE, mpi_rank+1, granurality, MPI_COMM_WORLD, &recv_req[granurality]);

            
            for (int step = 0; step < granurality; ++step)
            {
                int initial_index = ((x_size - 1)*i_jump) + (step * snd_data_size);
                int index_3d_initial = (step * snd_data_size);
                MPI_Wait(&recv_req[step], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int offset = 0; offset < snd_data_size; offset += number_of_densities)
                {
                    int index_aux = initial_index + offset;
                    int index_3d = index_3d_initial + offset;
                    for (int d = 0; d < k_jump; d++)
                    {
                        densities[index_aux + d] -= thomas_cx[x_size - 1][d] * block3d[index_3d + d];
                    }

                    for (int i = x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < k_jump; d++)
                        {
                            densities[index_dec + d] -= thomas_cx[i][d] * densities[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (mpi_rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&(densities[step * snd_data_size]), snd_data_size, MPI_DOUBLE, mpi_rank - 1, step, MPI_COMM_WORLD, &send_req[step]);
                    // cout << "Rank " << rank << " has send" << endl;
                }
            }
            if (snd_data_size_last != 0)
            {
                int initial_index = ((x_size - 1)*i_jump) + (granurality * snd_data_size);
                int index_3d_initial = (granurality * snd_data_size);
                MPI_Wait(&recv_req[granurality], MPI_STATUS_IGNORE);
                #pragma omp parallel for
                for (int offset = 0; offset < snd_data_size_last; offset += number_of_densities)
                {
                    int index_aux = initial_index + offset;
                    //int index = j * M.j_jump + k * M.k_jump + (M.mesh.x_coordinates.size() - 1) * M.i_jump;
                    int index_3d = index_3d_initial + offset;
                    // naxpy(&(*M.microenvironment)[n], M.thomas_cx[M.mesh.x_coordinates.size() - 1], block3d[k][j]);
                    for (int d = 0; d < k_jump; d++)
                    {
                        densities[index_aux + d] -= thomas_cx[x_size - 1][d] * block3d[index_3d + d];
                    }

                    for (int i = x_size - 2; i >= 0; i--)
                    {

                        int index_dec = index_aux - i_jump;
                        // naxpy(&(*M.microenvironment)[n], M.thomas_cx[i], (*M.microenvironment)[n + M.i_jump]);
                        for (int d = 0; d < k_jump; d++)
                        {
                            densities[index_dec + d] -= thomas_cx[i][d] * densities[index_aux + d];
                        }
                        index_aux = index_dec;
                        
                    }
                }
                if (mpi_rank > 0)
                {
                    MPI_Request aux;
                    MPI_Isend(&(densities[granurality * snd_data_size]), snd_data_size_last, MPI_DOUBLE, mpi_rank - 1, granurality, MPI_COMM_WORLD, &send_req[granurality]);
                }
            }
        }
    end_time = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();

    if (mpi_rank == 0)
        file << duration_us << ",";
    start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions_boundaries();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    //Y-diffusion
    start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for (int k = 0; k < z_size; k++)
    {
        for (int i = 0; i < x_size; i++)
        {

            int index = i * i_jump + k * k_jump;
            //(*M.p_density_vectors)[n] /= M.thomas_denomy[0];
            for (int d = 0; d < number_of_densities; d++)
            {
                densities[index + d] /= thomas_denomy[0][d];
            }

            for (int j = 1; j < y_size; j++)
            {

                int index_inc = index + j_jump;
                // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_j_jump]);
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_inc + d] += thomas_constant1[d] * densities[index + d];
                }
                //(*M.p_density_vectors)[n] /= M.thomas_denomy[j];
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_inc + d] /= thomas_denomy[j][d];
                }
                index = index_inc;
            }

            // back substitution

            index = i * i_jump + k * k_jump + (j_jump * (y_size - 1));
            for (int j = y_size - 2; j >= 0; j--)
            {

                int index_dec = index - j_jump;
                // naxpy(&(*M.p_density_vectors)[n], M.thomas_cy[j], (*M.p_density_vectors)[n + M.thomas_j_jump]);
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_dec + d] -= thomas_cy[j][d] * densities[index + d];
                }
                index = index_dec;
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    if (mpi_rank == 0)
        file << duration_us << ",";
    apply_dirichlet_conditions_boundaries();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    //Z-diffusion
    //cout << "Peto aqui 3!" << endl;

    start_time = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < x_size; i++)
        {
        for (int j = 0; j < y_size; j++)
            {
            int index = i * i_jump + j * j_jump;
            //(*M.p_density_vectors)[n] /= M.thomas_denomz[0];
            for (int d = 0; d < number_of_densities; d++)
            {
                densities[index + d] /= thomas_denomz[0][d];
            }

            // should be an empty loop if mesh.z_coordinates.size() < 2
            for (int k = 1; k < z_size; k++)
            {

                int index_inc = index + k_jump;
                // axpy(&(*M.p_density_vectors)[n], M.thomas_constant1, (*M.p_density_vectors)[n - M.thomas_k_jump]);
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_inc + d] += densities[d] * densities[index + d];
                }
                //(*M.p_density_vectors)[n] /= M.thomas_denomz[k];
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_inc + d] /= thomas_denomz[k][d];
                }

                index = index_inc;
            }

            index = i * i_jump + j * j_jump + (k_jump * (z_size - 1));
            for (int k = z_size - 2; k >= 0; k--)
            {

                int index_dec = index - k_jump;
                // naxpy(&(*M.p_density_vectors)[n], M.thomas_cz[k], (*M.p_density_vectors)[n + M.thomas_k_jump]);
                for (int d = 0; d < number_of_densities; d++)
                {
                    densities[index_dec + d] -= thomas_cz[k][d] * densities[index + d];
                }
                index = index_dec;
            }
        }
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
    if (mpi_rank == 0)
        file << duration_us << ",";
    
    start_time = std::chrono::high_resolution_clock::now();
    apply_dirichlet_conditions_boundaries();
    end_time = std::chrono::high_resolution_clock::now();
    apply_us += std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count();
     if (mpi_rank == 0)
            file << apply_us/4.0 << std::endl;

}