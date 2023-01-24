# Compiler
MPIFC = mpifort   
FC = gfortran  
IC = ifort

# The CFLAGS variable sets compiler flags:
#  -g                            compile with debug information
#  -Wall                         give verbose compiler warnings
#  -O0                           do not optimize generated code
#  -O2 -ftree-vectorize          apply level 2 optimization with auto-vectorization only
#  -O3                           apply level 3 optimizations
#  -std=gnu99                    use the GNU99 standard language definition
#  -fopt-info-all=opt_report.all save optimization report to file
# -fopenmp (-qopenmp)            enable openmp for gnu (intel)
# -fcheck=all                    array bounds checking (gnu compiler)
# -CB                            array bounds checking (intel compiler)
# -ffpe-trap=invalid             stops program at first NaN (gnu compiler)
# -fpe0                          stops program at first NaN (intel compiler)

CFLAGS = -O3  -qopenmp #-fpe0 

OBJS = constants_mod.o grid_arrays_mod.o fp_solver_mod.o driver.o

fp_test : $(OBJS)
	$(IC) $(CFLAGS) -o fp_test $(OBJS)
	
constants_mod.o : ../src/constants_mod.f90 
	$(IC) $(CFLAGS) -c ../src/constants_mod.f90

grid_arrays_mod.o : ../src/grid_arrays_mod.f90 constants_mod.o 
	$(IC) $(CFLAGS) -c ../src/grid_arrays_mod.f90
    
fp_solver_mod.o : ../src/fp_solver_mod.f90  constants_mod.o
	$(IC) $(CFLAGS) -c ../src/fp_solver_mod.f90
    
driver.o : ../src/driver.f90 constants_mod.o grid_arrays_mod.o fp_solver_mod.o
	$(IC) $(CFLAGS) -c ../src/driver.f90
    
clean:
	rm  -f *.o *.mod 