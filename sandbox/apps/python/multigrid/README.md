**Benchmarks**: 
 1. Jacobi 2D
 2. Jacobi 3D
 3. NAS-MG

**Required Parameters**
The following parameters can be changed in the Makefile for each benchmark. 
 # Multigrid Cycle : V / W                 CYCLE='V'
 # Multigrid levels                        L=3      
 # coarse-grid size along each dimension   SIZE=31  
 # V / W Cycle iterations                  NIT=10   
 # pre-smoothing steps                     NU1=10   
 # post-smoothing steps                    NU2=0    
 # corase-smoothing steps                  NUC=0    

**How to run**
In order to run any benchmark, navigate to its folder and run make command.
Example:

$ cd jacobi2d
$ cd make

The optimized code is available in the \*.cpp file and it also creates a shared library (\*.so file).


