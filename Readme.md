BioFVM-B Diffusion-decay solver 3D test bench 

    Author: Jose Luis Estragués Muñoz

Repository structure

BioFVM-B DDT/tests -> 
    *Tridiagonals.cpp: Test to evaluate the impact of new datastructure
    *Loops.cpp: Test to evaluate loop reordering
    *vector_256d.cpp: Test to evaluate the enhanced vectorization with AVX-256D
    *vector_512d.cpp: Test to evaluate the enhanced vectorization with AVX-512D

BioFVM-B DDT/output -> Folder to output from the simulations

BioFVM-B DDT/timing -> CSV target folder

BioFVM-B DDT/scripts -> Useful scripts for Marenostrum 5 supercomputer @ BSC

BioFVM-B DDT/kernels -> Here your kernels
