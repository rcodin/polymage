import ctypes
import numpy as np
from cv2 import *
import sys

sys.path.insert(0, "../")
from structs import *

def harris_cv(frame, lib_func):
    gray = cvtColor(frame, COLOR_BGR2GRAY)
    gray = np.float32(gray) / 4.0
    res = cornerHarris(gray, 3, 3, 0.04)
    return res

def harris_polymage(frame, lib_func):
    rows = frame.shape[0]
    cols = frame.shape[1]
    res = np.empty((rows, cols), np.float32)
    lib_func(ctypes.c_int(cols-2), \
             ctypes.c_int(rows-2), \
             ctypes.c_void_p(frame.ctypes.data), \
             ctypes.c_void_p(res.ctypes.data))
    return res

def add_harris_app(app_id):
    # 1. modes
    modes = [ModeType.CV2, ModeType.P_NAIVE, ModeType.P_OPT]
    # 2. modes that need shared library
    lib_modes = [ModeType.P_NAIVE, ModeType.P_OPT]

    # 3. set python functions from frame_process.py
    app_func_map = {}
    app_func_map[ModeType.CV2] = harris_cv
    app_func_map[ModeType.P_NAIVE] = harris_polymage
    app_func_map[ModeType.P_OPT] = harris_polymage

    # 4. create an App object
    app_dir = os.path.dirname(os.path.realpath(__file__))
    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app

