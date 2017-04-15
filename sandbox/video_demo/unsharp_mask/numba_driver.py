from unsharp_numba import *
#from numba import jit

#@jit("float32[::](uint8[::],none)",nogil = True, cache = True)
def unsharp_numba_driver(frame,lib_func):
    res = pipeline_numba(frame)
    return res
