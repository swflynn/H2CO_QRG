#==============================================================================#
FC=gfortran
FFLAGS=-O3
#FFLAGS= -g -O0 -Wall -Wextra -pedantic -fcheck=all -fbacktrace
#==============================================================================#
all: grid
#==============================================================================#
#                       Compile all Fortran Files
#==============================================================================#
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^
#==============================================================================#
#                               DGB Cluster
#==============================================================================#
#grid: sobol.o morse_greedy.o		#removed sobol for convenient testing
grid: morse_greedy.o
	$(FC) -o $@ $^ 
#==============================================================================#
clean:
	rm -rf *.o *.mod grid
