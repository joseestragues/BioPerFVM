PROGRAM_NAME := y_diffusion
PROGRAM_NAME_EX := y_diffusion_extrae

CC := mpic++

ARCH := native

LIBEX=-L/apps/BSCTOOLS/extrae/latest/openmpi_4_0_1/lib/ -lomptrace

FLAGS := -march=$(ARCH) -O3 -fomit-frame-pointer -fopenmp  -std=c++11 -fpermissive -g

COMPILE_COMMAND := $(CC) $(FLAGS)

Tridiagonals_OBJECTS := utils.o microenvironment.o kernels.o


tridiagonals: ./tests/tridiagonals.cpp $(Tridiagonals_OBJECTS)
	    $(COMPILE_COMMAND) -o tridiagonals tridiagonals.cpp $(Tridiagonals_OBJECTS)

#tridiagonals_ex: tridiagonals.cpp 
#	    $(COMPILE_COMMAND) -DEXTRAE_LINK -o tridiagonals_ex tridiagonals.cpp ${LIBEX}

tridiagonals_512: ./tests/tridiagonals_512.cpp $(Tridiagonals_OBJECTS)
	    $(COMPILE_COMMAND) -o tridiagonals_512 tridiagonals_512.cpp $(Tridiagonals_OBJECTS)

#tridiagonals_512_ex: tridiagonals_512.cpp 
#	    $(COMPILE_COMMAND) -DEXTRAE_LINK -o tridiagonals_512_ex tridiagonals_512.cpp ${LIBEX}

novector: novector.cpp $(Tridiagonals_OBJECTS)
	    $(COMPILE_COMMAND) -o novector novector.cpp $(Tridiagonals_OBJECTS)

#novector_ex: novector.cpp 
#	    $(COMPILE_COMMAND) -DEXTRAE_LINK -o novector_ex novector.cpp ${LIBEX}

utils.o: utils.cpp
	$(COMPILE_COMMAND) -c utils.cpp

microenvironment.o: microenvironment.cpp
	$(COMPILE_COMMAND) -c microenvironment.cpp


#Compile kernels for Diffussion Decay 
solver.o: ./kernels/solver.cpp
	$(COMPILE_COMMAND) -c ./kernels/solver.cpp

solver_256d.o: ./kernels/solver_256d.cpp
	$(COMPILE_COMMAND) -c ./kernels/solver_256d.cpp

solver_512d.o: ./kernels/solver_512d.cpp
	$(COMPILE_COMMAND) -c ./kernels/kernels_512d.cpp



clean:
	rm -f tridiagonals
	rm -f *.o
