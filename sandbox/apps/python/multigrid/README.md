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

 7. coarse-smoothing steps : NUC=0    


**How to run**

Run the following commands before running any of the experiments.

$ export KMP\_PLACE\_THREADS=24c,1t

$ export OMP\_NUM\_THREADS=24

In order to run any benchmark, navigate to its directory and run 'make'. 

Example:

$ cd jacobi2d

$ make

The make command optionally takes the following parameters:

1. all - same as running make. (reads the polymage code, creates a new C file and compiles and executes the code.)

2. new - same as the previous option.

3. tune - uses the app\_tuner.py file and tunes for various configurations of group and tile size.

4. existing - compiles and runs an existing optimized C file.

5. ready - runs an already compiled version of the optimized C code.

6. clean/cleaner - deletes the shared object and the optimized C file. 

The optimized code is written to the \*.cpp file, and it is also turned into 
a shared library (\*.so file).

The polymage DSL code is written in the file named polymage\_\*.py
