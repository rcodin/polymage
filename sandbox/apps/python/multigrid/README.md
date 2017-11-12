**How to run Multigrid benchmarks**
-------------------------------


**Install PolyMage**

Instructions on general pre-requisites to run Polymage apps are here: 
https://bitbucket.org/udayb/polymage

*Note*:

- If you encounter a 'pygraphviz not found' error, uninstall pygraphviz and 
  install with the below options

> sudo pip3 install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"


**Benchmarks**: 

- Jacobi 2D

- Jacobi 3D

- NAS-MG

**Required Parameters**

The various configurations as mentioned in the paper can be obtained by varying the following parameters in the Makefile.

 1. Multigrid Cycle (V / W): CYCLE='V'

 2. Multigrid levels: L=3      

 3. Coarse-grid size along each dimension: SIZE=31  

 4. V / W Cycle iterations: NIT=10   

 5. Pre-smoothing steps: NU1=10   

 6. Post-smoothing steps: NU2=0    

 7. Coarse-smoothing steps: NUC=0    

 8. Number of runs (for timing purposes): RUNS=1 (default)


**How to run**

Run the following commands before running any of the experiments.

> $ export KMP\_PLACE\_THREADS=24c,1t

> $ export OMP\_NUM\_THREADS=24

In order to run any benchmark, navigate to its directory and run 'make'. 

'make' executes polymg-opt+, polymg-opt, handopt and handopt-pluto 
configurations.

Example:

> $ cd sandbox/apps/python/multigrid

> $ cd jacobi2d

> $ make

The various make targets available and what they correspond to are listed 
below:

- all (executes polymg-opt+, polymg-opt, handopt and handopt-pluto 
  configurations)

- polymg-opt-plus - executes polymg-opt+ configurations.

- polymg-opt - executes polymg-opt configurations.

- handopt - executes hand optimized code (for Jacobi2D and Jacobi3D 
  only).

- handopt-pluto - executes hand optimized code with Pluto's diamond 
  tiling (for Jacobi2D and Jacobi3D only).

- tune - uses the app\_tuner.py file and tunes for various group and tile size for polymg code.

- clean/cleaner - deletes the shared object and the optimized Polymage generated C file. 

The optimized Polymage code is written to the \*.cpp file, and it is also turned into 
a shared library (\*.so file).

The Polymage DSL code is written in the file named polymage\_\*.py

To get results for polymg-dtile-opt+ (Multigrid smoothing optimized with 
diamond tiling), libpluto (included as a submodule) is needed.

- change the branch to Tstencil.

>   $ git checkout origin/Tstencils

>   $ make

**Output**

Comment out the "--timer" option in the Makefile to see error for each 
V-cycle iteration. In such a scenario (where performance timing isn't 
being performed), it's meaningful to use RUNS=1.
