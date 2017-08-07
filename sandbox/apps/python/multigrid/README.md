**Benchmarks**: 

 1. Jacobi 2D

 2. Jacobi 3D

 3. NAS-MG

**Required Parameters**

The following parameters can be changed in the Makefile for each benchmark. 

 1. Multigrid Cycle (V / W) : CYCLE='V'

 2. Multigrid levels : L=3      

 3. coarse-grid size along each dimension : SIZE=31  

 4. V / W Cycle iterations : NIT=10   

 5. pre-smoothing steps : NU1=10   

 6. post-smoothing steps : NU2=0    

 7. corase-smoothing steps : NUC=0    


**How to run**

In order to run any benchmark, navigate to its folder and run make command.

Example:

$ cd jacobi2d

$ cd make

The optimized code is available in the \*.cpp file and it also creates a shared library (\*.so file).


