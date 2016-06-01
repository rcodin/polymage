# from ctypes import cdll, Structure, c_int, c_double, c_uint
from cffi import FFI
import islpy as isl


class PlutoFFI(object):
    def __init__(self):
        self.ffi = FFI()

        self.ffi.cdef("""
            struct plutoOptions{

    /* To tile or not? */
    int tile;

    /* Intra-tile optimization */
    int intratileopt;

    /* Load-balanced tiling */
    int lbtile;

    /* Load-balanced tiling (one dimensional concurrent start)*/
    int partlbtile;

    /* Extract scop information from libpet*/
    int pet;

    /* dynamic scheduling 
     * using Synthesized Runtime Interface */
    int dynschedule;

    /* dynamic scheduling - previous technique of 
     * building the entire task graph in memory 
     * using Intel TBB Flow Graph scheduler */
    int dynschedule_graph;

    /* dynamic scheduling - previous technique of 
     * building the entire task graph in memory 
     * using a custom DAG scheduler */
    // no longer maintained
    int dynschedule_graph_old;

    /* consider transitive dependences between tasks */
    int dyn_trans_deps_tasks;

    /* parallelization */
    int parallel;

    /* prefer pure inner parallelism to pipelined parallelism */
    int innerpar;

    /* Automatic unroll/unroll-jamming of loops */
    int unroll;

    /* unroll/jam factor */
    int ufactor;

    /* Enable or disable post-transformations to make code amenable to
     * vectorization (default - enabled) */
    int prevector;

    /* consider RAR dependences */
    int rar;

    /* Decides the fusion algorithm (MAXIMAL_FUSE, NO_FUSE, or SMART_FUSE) */
    int fuse;

    /* for debugging - print default cloog-style total */
    int scancount;

    /* parameters will be assumed to be at least this much */
    /* This is appended to the context passed to cloog */
    int codegen_context;

    /* Loop depth (1-indexed) to force as parallel */
    int forceparallel;

    /* multiple (currently two) degrees of pipelined parallelism */
    int multipar;

    /* Tile for L2 too */
    /* By default, only L1 tiling is done; under parallel execution, every
     * processor executes a sequence of L1 tiles (OpenMP adds another blocking
     * on the parallel loop). With L2 tiling, each processor executes a
     * sequence of L2 tiles and barrier is done after a group of L2 tiles is
     * exectuted -- causes load imbalance due to pipe startup when problem
     * sizes are not huge */
    int l2tile;


    /* NOTE: --ft and --lt are to manually force tiling depths */
    /* First depth to tile (starting from 0) */
    int ft;
    /* Last depth to tile (indexed from 0)  */
    int lt;

    /* Output for debugging */
    int debug;

    /* More debugging output */
    int moredebug;

    /* Not implemented yet: Don't output anything unless something fails */
    int quiet;

    /* Pure polyhedral unrolling (instead of postpass) */
    int polyunroll;

    /* Identity transformation */
    int identity;

    /* Identity transformation */
    int identity_data_dist;

    /* Generate scheduling pragmas for Bee+Cl@k */
    int bee;

    /* Force this for cloog's -f */
    int cloogf;

    /* Force this for cloog's -l */
    int cloogl;

    /* Enable cloog's -sh (simple convex hull) */
    int cloogsh;

    /* Enable cloog's -backtrack */
    int cloogbacktrack;

    /* Use isl to compute dependences (default) */
    int isldep;

    /* Use candl to compute dependences */
    int candldep;

    /* Access-wise dependences with ISL */
    int isldepaccesswise;

    /* Coalesce ISL deps */
    int isldepcoalesce;

    /* Compute lastwriter for dependences */
    int lastwriter;

    /* DEV: Don't use cost function */
    int nodepbound;

    /* hard upper bound for transformation coefficients */
    int coeff_bound;

    /* Ask candl to privatize */
    int scalpriv;

    /* No output from Pluto if everything goes right */
    int silent;

    /* Read input from a .scop file */
    int readscop;

    /* Use PIP as ilp solver. */
    int pipsolve;

    /* Use isl as ilp solver. */
    int islsolve;

    int glpksolve;

    /* Index set splitting */
    int iss;

    int distmem;

    /*  adding support to generate opencl code */
    int opencl; 

    /* use multi-level distribution function */
    /* for dynamic scheduling or distributed-memory code */
    /* OFF by default */
    int multi_level_distribution;

    int commopt;

    /*Communication code generation using flow-out partitioning */
    int commopt_fop;
    /* generate code to choose between unicast pack and multicast pack 
     * for each partition at runtime */
    int fop_unicast_runtime;

    /*Communication code generation using flow-out intersection flow-in */
    int commopt_foifi;

    /*Report communication for distributed memory*/
    int timereport;

    /* if true, variables are not declared globally
     * but each variable's declaration is provided 
     * through the macro '#define __DECLARATION_OF_<variable-name> <declaration>'*/
    int variables_not_global;

    int data_dist;
    int verify_output;

    int mpiomp;
    int fusesends;
    int blockcyclic;
    int cyclesize;

    //enables mod eliminate and data ptr optimization for data tiling
    int data_tile_opt;

    //Propagates the bounding box constraints across non fused loops
    int global_opt;

    //auto compute pi
    int compute_pi;

    //max number of tiles to be used while computing pi
    int num_tiles_per_dim;

    //number of initial partitions used while computing pi
    int num_inital_partitions;

    /* Output file name supplied from -o */
    char *out_file;

    /* Polyhedral compile time stats */
    int time;

    /* Experimental optimizations to make Pluto faster/scalable */
    int fast;

    /* Eliminate Farkas multipliers using PolyLib */
    int efup;

    /* fast linear independence check */
    int flic;

    /* SCoP number when processing multiple SCoPs per file */
    int scopnum;
};
typedef struct plutoOptions PlutoOptions;


/* Fusion options for options->fuse */

/* Do not fuse across SCCs */
#define NO_FUSE 0
/* Geared towards maximal fusion, but not really maximal fusion */
#define MAXIMAL_FUSE 1
/* Something in between the above two */
#define SMART_FUSE 2

PlutoOptions *pluto_options_alloc();
void pluto_options_free(PlutoOptions *);

void pluto_schedule_str(const char *domains_str,
        const char *dependences_str,
        char** schedule_str_buffer_ptr,
        int *schedule_strlen,
        PlutoOptions *options);

""")
        self.sharedobj = self.ffi.dlopen('libpluto.so')

        print('Loaded lib {0}'.format(self.sharedobj))

    def options_alloc(self):
        return self.sharedobj.pluto_options_alloc()

    def schedule(self, ctx, domains, dependences, options):
        assert isinstance(domains, isl.UnionSet)
        assert isinstance(dependences, isl.UnionMap)
        #domains_str = self.ffi.new("char[]", domains.to_str().encode('utf-8'))
        #dependences_str = self.ffi.new("char[]", dependences.to_str().encode('utf-8'))
        
        domains_str = domains.to_str().encode('utf-8')
        dependences_str = dependences.to_str().encode('utf-8')
        schedule_strbuf = self.ffi.new("char **")
        schedule_strlen = self.ffi.new("int *")

        self.sharedobj.pluto_schedule_str(domains_str, dependences_str,
                                          schedule_strbuf, schedule_strlen,
                                          options)

        # import pdb
        # pdb.set_trace()
        assert schedule_strbuf[0] != self.ffi.NULL, \
            ("unable to get schedule from PLUTO")

        schedule_str = self.ffi.string(schedule_strbuf[0]).decode('utf-8')
        schedule = isl.UnionMap.read_from_str(ctx, schedule_str)
        return schedule


if __name__ == "__main__":
    ffi = PlutoFFI()
    ctx = isl.Context.alloc()
    opts = ffi.options_alloc()
    domains = isl.UnionSet.read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= 0 and i0 <= p_0}")
    deps = isl.UnionMap.read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : (exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and 8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 >= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 = [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 + 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 >= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + 1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 >= 0 and 8i1 >= 8192 - p_5 + p_7))}")
    sched = ffi.schedule(ctx, domains, deps, opts)
