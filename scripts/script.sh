#!/bin/bash



export OMP_NUM_THREADS=$threads
./tridiagonals 5000 10 2 1 ./outputs/${threads}_threads.csv #1> small.log 2> small.log
