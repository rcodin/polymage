**How to run Multigrid benchmarks**
-------------------------------


**Install polymage**

The instructions to install Polymage is present here - https://bitbucket.org/udayb/polymage

*Note*:
- Install Pluto to run the benchmarks in polymg-dtile-opt+ configuration
- In case you get pygraphviz not found error, uninstall pygraphviz and install with the below options
> sudo pip3 install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"


**Benchmarks**: 

- Jacobi 2D

- Jacobi 3D

- NAS-MG

**Required Parameters**

The various configurations as mentioned in the paper can be obtained by varying the following parameters in the Makefile.

 1. Multigrid Cycle (V / W) : CYCLE='V'

 2. Multigrid levels : L=3      

 3. coarse-grid size along each dimension : SIZE=31  

 4. V / W Cycle iterations : NIT=10   

 5. pre-smoothing steps : NU1=10   

 6. post-smoothing steps : NU2=0    

 7. coarse-smoothing steps : NUC=0    


**How to run**

Run the following commands before running any of the experiments.

> $ export KMP\_PLACE\_THREADS=24c,1t

> $ export OMP\_NUM\_THREADS=24

In order to run any benchmark, navigate to its directory and run 'make'. 

Example:

> $ cd sandbox/apps/python/multigrid

> $ cd jacobi2d

> $ make

The make command optionally takes the following parameters:

- all - same as running make. (reads the polymage code, creates a new C file and compiles and executes the code.)

- new - same as the previous option.

- tune - uses the app\_tuner.py file and tunes for various configurations of group and tile size.

- existing - compiles and runs an existing optimized C file.

- ready - runs an already compiled version of the optimized C code.

- clean/cleaner - deletes the shared object and the optimized C file. 

The optimized code is written to the \*.cpp file, and it is also turned into 
a shared library (\*.so file).

The polymage DSL code is written in the file named polymage\_\*.py

To get results for polymg-opt+:

- Storage optimizations are enabled by default. Running make command gives the result for polymg-opt+

To get results for polymg-opt:

- comment OPT\_ARGS (line no: 19) from the Makefile and run the make command.

To get results for polymg-dtile-opt+:

- change the branch to Tstencil.

>   $ git fetch && git checkout origin/Tstencils

>   $ make
