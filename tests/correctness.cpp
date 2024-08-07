#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <omp.h>
#include <mpi.h>
#include "../microenvironment.h"
#include "../omp/microenvironment_omp.h"
#include "../utils.h"



using namespace std;


double error = 0.00001;

int main(int argc, char* argv[]) {
  
    #ifdef EXTRAE_LINK
        Extrae_shutdown();
    #endif

    double side, delta;
    int densities;
    int factor;
    string path;
    if (argc == 6) {
		side = atof(argv[1]);
        delta = atof(argv[2]);
        densities = stoi(argv[3]);
        factor = stoi(argv[4]);
        path = argv[5];

        
	} 
	else {
		cout << "Error: Expected arguments are: " << endl;
        cout << "   1. Simulation cube side (um)" << endl;
        cout << "   2. Discretization  (um)" << endl;
        cout << "   3. Number of densities" << endl;
        cout << "   4. Granurality" << endl;
        cout << "   5. Path to csv timing" << endl;
		exit(1);
	}
    int size, wrank;
    /* Declaracion de variables */
    MPI_Init(&argc, &argv);
    /* Reparto de contenido */
    /* Bucle principal del programa */
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);
    
    if (wrank == 0 ) {    
        cout << "Cube side:  " << side << " um" << endl;
        cout << "Number of voxels per dimension " << side/delta << endl;
        cout << "Total voxels " << side/delta * side/delta * side/delta << endl;
        cout << "Number of substrates: " << densities << endl;
        cout << "Path to CSV " << path << " number of threads" << endl;
        cout << "Execution with " << omp_get_max_threads() << " number of threads" << endl;
    }

    
    microenvironment_omp m_omp;

    m_omp.init_diffusion_coeficients(side, delta, densities);
    m_omp.init_densities();
    microenvironment m;
    
    auto start = std::chrono::high_resolution_clock::now();
    //m.factor = factor;
    m.granurality = factor;
    m.timing_csv = path;
    m.init_diffusion_coeficients(side, delta, densities, wrank, size,1, MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    if (wrank == 0) std::cout << "Time taken by init microenvironment: " << duration.count() << " ms" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    m.init_densities();
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    if (wrank == 0) std::cout << "Time taken by init densities: " << duration.count() << " ms" << std::endl;

    bool identical = m.compare_microenvironment(m_omp);
    std::cout << "Rank " << wrank << " is identical " << identical << std::endl;

    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; ++i)
        m_omp.diffusion_decay_3D_solver();
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    if (wrank == 0) std::cout << "  BioFVM: Diffusion decay 100 executions: " << duration.count() << " ms" << std::endl;

     #ifdef EXTRAE_LINK
        Extrae_restart();
    #else
        start = std::chrono::high_resolution_clock::now();
    #endif

    for (int i = 0; i < 100; ++i)
        m.diffusion_decay_3D_solver();
    #ifdef EXTRAE_LINK
        Extrae_shutdown();
    #else
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        if (wrank == 0) std::cout << "  BioFVM-B: Diffusion decay 100 executions: " << duration.count() << " ms" << std::endl;
    #endif
    

    identical = m.compare_microenvironment(m_omp);
    std::cout << "Rank " << wrank << " is identical " << identical << std::endl;

    if (identical != 1) {
        string filename = "./m_mpi.txt";
        m.print_voxels_densities(&filename);
        filename = "./m_omp.txt";
        m_omp.print_voxels_densities(&filename);
    }

    MPI_Finalize();
    return 0;
}
