import ctypes
import numpy as np
from cv2 import *
import sys

sys.path.insert(0, "../")
from structs import *

def bilateral_polymage(frame, lib_func):
    rows = frame.shape[0]
    cols = frame.shape[1]
    total_pad = 56
    res = np.empty((rows, cols), np.float32)
    lib_func(ctypes.c_int(cols+total_pad), \
             ctypes.c_int(rows+total_pad), \
             ctypes.c_void_p(frame.ctypes.data), \
             ctypes.c_void_p(res.ctypes.data))
    return res

def add_bilateral_app(app_id):
    # 1. modes
    modes = [ModeType.P_NAIVE, ModeType.P_OPT]
    # 2. modes that need shared library
    lib_modes = modes

    # 3. set python functions from frame_process.py
    app_func_map = {}
    app_func_map[ModeType.P_NAIVE] = bilateral_polymage
    app_func_map[ModeType.P_OPT] = bilateral_polymage

    # 4. create an App object
    app_dir = os.path.dirname(os.path.realpath(__file__))
    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app
