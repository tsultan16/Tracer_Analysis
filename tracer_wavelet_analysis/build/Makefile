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

CFLAGS = -O3  -qopenmp -fpe0 

OBJS = constants_mod.o grid_arrays_mod.o fft_mod.o fwt_mod.o readfile_mod.o tracer_rec_mod.o analysis_main.o

mode_decomp : $(OBJS)
	$(IC) $(CFLAGS) -o mode_decomp $(OBJS)
	
constants_mod.o : ../src/constants_mod.f90 
	$(IC) $(CFLAGS) -c ../src/constants_mod.f90

grid_arrays_mod.o : ../src/grid_arrays_mod.f90 constants_mod.o
	$(IC) $(CFLAGS) -c ../src/grid_arrays_mod.f90
				
fft_mod.o : ../src/fft_mod.f90 constants_mod.o 
	$(IC) $(CFLAGS) -c ../src/fft_mod.f90	

fwt_mod.o : ../src/fwt_mod.f90 constants_mod.o grid_arrays_mod.o
	$(IC) $(CFLAGS) -c ../src/fwt_mod.f90

readfile_mod.o : ../src/readfile_mod.f90  constants_mod.o  
	$(IC) $(CFLAGS) -c ../src/readfile_mod.f90
    
tracer_rec_mod.o : ../src/tracer_rec_mod.f90  constants_mod.o  grid_arrays_mod.o fwt_mod.o
	$(IC) $(CFLAGS) -c ../src/tracer_rec_mod.f90
    

analysis_main.o : ../src/analysis_main.f90 constants_mod.o grid_arrays_mod.o fft_mod.o fwt_mod.o readfile_mod.o tracer_rec_mod.o
	$(IC) $(CFLAGS) -c ../src/analysis_main.f90
    
clean:
	rm  -f *.o *.mod 