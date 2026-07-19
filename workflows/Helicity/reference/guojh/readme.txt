1. install mkl

2. compile mkl_poisson.f90 and mkl_dfti

    ifort -c /share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/include/mkl_dfti.f90 -I/share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/include
    ifort -c /share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/include/mkl_poisson.f90 -I/share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/include

3. change parameters in helicity.f90

4. ifort -o helicity helicity.f90 mkl_poisson.o mkl_dfti.o -I/share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/include -L/share/apps/intel/compilers_and_libraries_2017.1.132/linux/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

5.  ./helicity
