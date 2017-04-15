# from ctypes import cdll, Structure, c_int, c_double, c_uint
from cffi import FFI
import logging
import islpy as isl

# LOG CONFIG #
libpluto_logger = logging.getLogger("libpluto.py")
libpluto_logger.setLevel(logging.INFO)
LOG = libpluto_logger.log

_pluto_header_str = \
    """
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

        /* Decides the fusion algorithm (MAXIMAL_FUSE, NO_FUSE,
        or SMART_FUSE) */
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
         * processor executes a sequence of L1 tiles (OpenMP adds another
         * blocking on the parallel loop). With L2 tiling, each processor
         * executes a sequence of L2 tiles and barrier is done after a
         * group of L2 tiles is exectuted -- causes load imbalance due to pipe
         * startup when problem sizes are not huge */
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
         * through the macro '#define __DECLARATION_OF_<variable-name>
         * <declaration>'*/
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
            PlutoOptions *options);

    void pluto_schedules_strbuf_free(char *schedules_str_buffer);
"""


class PlutoOptions(object):
    """
    Options to be passed to PLUTO during schedules creation

    TODO: add more PLUTO options as desired
    """
    def __init__(self, pluto_ffi, raw_options_ptr):
        self._pluto_ffi = pluto_ffi
        self._raw_options_ptr = raw_options_ptr
        self._raw_options_ptr.parallel = 1
        self._raw_options_ptr.partlbtile = 1
        self._raw_options_ptr.lbtile = 1
        self._raw_options_ptr.tile = 1

    @property
    def partlbtile(self):
        """
        Load-balanced tiling (one dimensional concurrent start)
        """
        return self._raw_options_ptr.partlbtile

    @partlbtile.setter
    def partlbtile(self, partlbtile):
        """
        Set load-balanced tiling (one dimensional concurrent start) status

        Parameters
        ----------
        partlbtile : Bool
        """
        self._raw_options_ptr.partlbtile = 1 if partlbtile else 0

    def __del__(self):
        self._pluto_ffi._destroy_raw_options_ptr(self._raw_options_ptr)


class LibPluto(object):
    """
    Represents the FFI to libpluto. On construction, this loads
    libpluto.so and maintains a reference to it as long as it lives
    """
    def __init__(self):
        self.ffi = FFI()

        self.ffi.cdef(_pluto_header_str)
        self.sharedobj = self.ffi.dlopen('libpluto.so')

        return

    def create_options(self):
        """
        Creates a PlutoOptions object, which allows configuring PLUTO.

        Parameters
        ----------
        None

        Returns
        -------
        PlutoOptions to configure
        """
        return PlutoOptions(self, self.sharedobj.pluto_options_alloc())

    def _destroy_raw_options_ptr(self, raw_options_ptr):
        """
        Frees a raw C PlutoOptions* owned by Python PlutoOptions

        NOTE
        ----
        This function is internal, and should *only* be called by
        PlutoOptions
        """
        # HACK: don't free, try to figure out why there's a segfault
        # self.sharedobj.pluto_options_free(raw_options_ptr)

    def schedule(self, ctx, domains, dependences, options):
        self.map_input_translation = {}
        self.map_domain_tuples_translation = {}
        self.map_output_translations = {}

        if isinstance(domains, isl.BasicSet):
            domains = isl.UnionSet.from_basic_set(domains)

        assert isinstance(domains, isl.UnionSet)

        if isinstance(dependences, isl.BasicMap):
            dependences = isl.UnionMap.from_basic_map(dependences)
        assert isinstance(dependences, isl.UnionMap)
        assert isinstance(options, PlutoOptions)

        log_string1 = "domains:\n" + str(domains)
        log_string2 = "depdendences:\n" + str(dependences)
        LOG(log_level, log_string1)
        LOG(log_level, log_string2)

        domains_str = domains.to_str().encode('utf-8')
        dependences_str = dependences.to_str().encode('utf-8')
        schedule_strbuf_ptr = self.ffi.new("char **")

        self.sharedobj.pluto_schedule_str(domains_str, dependences_str,
                                          schedule_strbuf_ptr,
                                          options._raw_options_ptr)

        assert schedule_strbuf_ptr[0] != self.ffi.NULL, \
            ("unable to get schedule from PLUTO")

        schedule_str = self.ffi.string(schedule_strbuf_ptr[0]).decode('utf-8')
        schedule = isl.UnionMap.read_from_str(ctx, schedule_str)

        self.sharedobj.pluto_schedules_strbuf_free(schedule_strbuf_ptr[0])

        return schedule


# This is somewhat of a hack, just to run a "test" if this file is
# executed separately.
# TODO: move this to a separate file
if __name__ == "__main__":
    pluto = LibPluto()

    ctx = isl.Context.alloc()
    opts = pluto.create_options()
    opts.partlbtile = 1
    domains = isl.UnionSet.read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= 0 and i0 <= p_0}")
    deps = isl.UnionMap.read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : (exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and 8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 >= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 = [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 + 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 >= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + 1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 >= 0 and 8i1 >= 8192 - p_5 + p_7))}")
    sched = pluto.schedule(ctx, domains, deps, opts)
