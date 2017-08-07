**Benchmarks**: 

 1. Jacobi 2D

 2. Jacobi 3D

 3. NAS-MG

**Required Parameters**

The following parameters can be changed in the Makefile for each benchmark. 

 # Multigrid Cycle : V / W  &nbsp; &nbsp;          CYCLE='V'

 # Multigrid levels  &nbsp; &nbsp;                 L=3      

 # coarse-grid size along each dimension &nbsp;   SIZE=31  

 # V / W Cycle iterations   &nbsp; &nbsp;          NIT=10   

 # pre-smoothing steps      &nbsp; &nbsp;          NU1=10   

 # post-smoothing steps      &nbsp; &nbsp;         NU2=0    

 # corase-smoothing steps    &nbsp; &nbsp;         NUC=0    


**How to run**

In order to run any benchmark, navigate to its folder and run make command.

Example:

$ cd jacobi2d

$ cd make

The optimized code is available in the \*.cpp file and it also creates a shared library (\*.so file).


