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

'make' runs polymg-opt+, polymg-opt, handopt and handopt-pluto configurations.

Example:

> $ cd sandbox/apps/python/multigrid

> $ cd jacobi2d

> $ make

The make command can be run with the following arguments:

- all - same as running make. (runs polymg-opt+, polymg-opt, handopt and handopt-pluto configurations)

- new - runs polymg-opt+ configurations.

- polymg-opt-plus - runs polymg-opt+ configurations.

- polymg-opt - runs polymg-opt configurations.

- handopt - runs hand optimized code (for Jacobi2D and Jacobi3D only).

- handopt-pluto - runs hand optimized code with Pluto's diamond tiling (for Jacobi2D and Jacobi3D only).

- tune - uses the app\_tuner.py file and tunes for various group and tile size for polymg code.

- clean/cleaner - deletes the shared object and the optimized Polymage generated C file. 

The optimized Polymage code is written to the \*.cpp file, and it is also turned into 
a shared library (\*.so file).

The Polymage DSL code is written in the file named polymage\_\*.py

To get results for polymg-dtile-opt+:

- change the branch to Tstencil.

>   $ git fetch && git checkout origin/Tstencils

>   $ make
