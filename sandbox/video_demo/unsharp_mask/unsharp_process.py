import ctypes
import numpy as np
from cv2 import *
import sys

sys.path.insert(0, "../")
from structs import *

# unsharp mask parameters
thresh = 0.001
weight = 3

def unsharp_cv(frame, lib_func):
    kernel = np.array([[1,4,6,4,1]], np.float32) / 16.0
    blurx = filter2D(frame, -1, kernel)
    kernel = np.array([[0,0,1,0,0], [0,0,4,0,0], [0,0,6,0,0], [0,0,4,0,0], [0,0,1,0,0]], np.float32) / 16.0
    blury = filter2D(blurx, -1, kernel)
    sharpen = addWeighted(frame, (1+weight), blury, (-weight), 0)
    th, mask = threshold(absdiff(frame, blury), thresh, 1, THRESH_BINARY)
    mask = mask.astype(bool)
    res = frame
    np.copyto(res, sharpen, 'same_kind', mask)

    return res

def unsharp_polymage(frame, lib_func):
    rows = frame.shape[0]
    cols = frame.shape[1]
    res = np.empty((rows-4, cols-4, 3), np.float32)
    lib_func(ctypes.c_int(cols - 4), \
             ctypes.c_int(rows - 4), \
             ctypes.c_float(thresh), \
             ctypes.c_float(weight), \
             ctypes.c_void_p(frame.ctypes.data), \
             ctypes.c_void_p(res.ctypes.data))
    return res

def add_unsharp_app(app_id):
    # 1. modes
    modes = [ModeType.CV2, ModeType.P_NAIVE, ModeType.P_OPT]
    # 2. modes that need shared library
    lib_modes = [ModeType.P_NAIVE, ModeType.P_OPT]

    # 3. set python functions from frame_process.py
    app_func_map = {}
    app_func_map[ModeType.CV2] = unsharp_cv
    app_func_map[ModeType.P_NAIVE] = unsharp_polymage
    app_func_map[ModeType.P_OPT] = unsharp_polymage

    # 4. create an App object
    app_dir = os.path.dirname(os.path.realpath(__file__))
    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app

