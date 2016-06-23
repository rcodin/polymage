import ctypes
import numpy as np
from cv2 import *
import sys
from structs import *
from PIL import Image, ImageFilter
from unsharp_numba_version import unsharp_numba
from numba import jit

sys.path.insert(0, "../")

# unsharp mask parameters
thresh = 0.001
weight = 3

# PIL version
@jit("uint8[::](uint8[::], int64)", cache = True, nogil = True)
def unsharp_pil(frame, lib_func):
    im = Image.fromarray(frame)
    kernelx = (0,0,0,0,0,0,0,0,0,0,1,4,6,4,1,0,0,0,0,0,0,0,0,0,0)
    kernely = (0,0,1,0,0,0,0,4,0,0,0,0,6,0,0,0,0,4,0,0,0,0,1,0,0)
    blurx = im.filter(ImageFilter.Kernel((5,5),kernelx,scale = None, offset = 0))
    blury = blurx.filter(ImageFilter.Kernel((5,5),kernely,scale = None, offset = 0))
    sharpen = Image.blend(im,blury,-weight)
    #diff = ImageChops.difference(im,blury)
    """m = im.filter(ImageFilter.UnsharpMask(radius = 2, \
            percent = 130, threshold = 1))"""
    res = np.array(sharpen)
    return res

# OpenCV version
@jit("float32[::](uint8[::], int64)", cache = True, nogil = True)
def unsharp_cv(frame, lib_func):
    frame_f = np.float32(frame) / 255.0
    res = frame_f
    kernelx = np.array([1, 4, 6, 4, 1], np.float32) / 16
    kernely = np.array([[1], [4], [6], [4], [1]], np.float32) / 16
    blury = sepFilter2D(frame_f, -1, kernelx, kernely)
    sharpen = addWeighted(frame_f, (1 + weight), blury, (-weight), 0)
    th, choose = threshold(absdiff(frame_f, blury), thresh, 1, THRESH_BINARY)
    choose = choose.astype(bool)
    np.copyto(res, sharpen, 'same_kind', choose)
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
    modes = [ModeType.CV2, ModeType.P_NAIVE, ModeType.P_OPT, \
        ModeType.PIL, ModeType.NUMBA]
    # 2. modes that need shared library
    lib_modes = [ModeType.P_NAIVE, ModeType.P_OPT]

    # 3. set python functions from frame_process.py
    app_func_map = {}
    app_func_map[ModeType.CV2] = unsharp_cv
    app_func_map[ModeType.P_NAIVE] = unsharp_polymage
    app_func_map[ModeType.P_OPT] = unsharp_polymage
    app_func_map[ModeType.PIL] = unsharp_pil
    app_func_map[ModeType.NUMBA] = unsharp_numba

    # 4. create an App object
    app_dir = os.path.dirname(os.path.realpath(__file__))
    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app
