import ctypes
import numpy as np
from cv2 import *
import sys

sys.path.insert(0, "../")
from structs import *

from ctypes import *
import math
import random

def sample(probs):
    s = sum(probs)
    probs = [a/s for a in probs]
    r = random.uniform(0, 1)
    for i in range(len(probs)):
        r = r - probs[i]
        if r <= 0:
            return i
    return len(probs)-1

def c_array(ctype, values):
    arr = (ctype*len(values))()
    arr[:] = values
    return arr

class BOX(Structure):
    _fields_ = [("x", c_float),
                ("y", c_float),
                ("w", c_float),
                ("h", c_float)]

class IMAGE(Structure):
    _fields_ = [("w", c_int),
                ("h", c_int),
                ("c", c_int),
                ("data", POINTER(c_float))]

class METADATA(Structure):
    _fields_ = [("classes", c_int),
                ("names", POINTER(c_char_p))]

    

#lib = CDLL("/home/pjreddie/documents/darknet/libdarknet.so", RTLD_GLOBAL)
lib = CDLL("darknet/files/libdarknet.so", RTLD_GLOBAL)
lib.network_width.argtypes = [c_void_p]
lib.network_width.restype = c_int
lib.network_height.argtypes = [c_void_p]
lib.network_height.restype = c_int

predict = lib.network_predict
predict.argtypes = [c_void_p, POINTER(c_float)]
predict.restype = POINTER(c_float)

set_gpu = lib.cuda_set_device
set_gpu.argtypes = [c_int]

make_image = lib.make_image
make_image.argtypes = [c_int, c_int, c_int]
make_image.restype = IMAGE

make_boxes = lib.make_boxes
make_boxes.argtypes = [c_void_p]
make_boxes.restype = POINTER(BOX)

free_ptrs = lib.free_ptrs
free_ptrs.argtypes = [POINTER(c_void_p), c_int]

num_boxes = lib.num_boxes
num_boxes.argtypes = [c_void_p]
num_boxes.restype = c_int

make_probs = lib.make_probs
make_probs.argtypes = [c_void_p]
make_probs.restype = POINTER(POINTER(c_float))

detect = lib.network_predict
detect.argtypes = [c_void_p, IMAGE, c_float, c_float, c_float, POINTER(BOX), POINTER(POINTER(c_float))]

reset_rnn = lib.reset_rnn
reset_rnn.argtypes = [c_void_p]

load_net = lib.load_network
load_net.argtypes = [c_char_p, c_char_p, c_int]
load_net.restype = c_void_p

free_image = lib.free_image
free_image.argtypes = [IMAGE]

letterbox_image = lib.letterbox_image
letterbox_image.argtypes = [IMAGE, c_int, c_int]
letterbox_image.restype = IMAGE

load_meta = lib.get_metadata
lib.get_metadata.argtypes = [c_char_p]
lib.get_metadata.restype = METADATA

load_image = lib.load_image_color
load_image.argtypes = [c_char_p, c_int, c_int]
load_image.restype = IMAGE

rgbgr_image = lib.rgbgr_image
rgbgr_image.argtypes = [IMAGE]

predict_image = lib.network_predict_image
predict_image.argtypes = [c_void_p, IMAGE]
predict_image.restype = POINTER(c_float)

network_detect = lib.network_detect
network_detect.argtypes = [c_void_p, IMAGE, c_float, c_float, c_float, POINTER(BOX), POINTER(POINTER(c_float))]

test_detector = lib.test_detector
test_detector.argtypes = [c_char_p, c_void_p, c_char_p, c_float, c_float, c_char_p, c_int]

# void test_detector(char *datacfg, char *cfgfile, char *weightfile, char *filename, float thresh, float hier_thresh, char *outfile, int fullscreen)

net = load_net("darknet/files/yolo.cfg", "darknet/files/yolo.weights", 0)
# meta = load_meta("/home/e0358-3/asst-2/darknet/cfg/coco.data")

# void test_detector(char *datacfg, void *net_ptr, char *filename, float thresh, float hier_thresh, char *outfile, int fullscreen)

def classify(net, meta, im):
    out = predict_image(net, im)
    res = []
    for i in range(meta.classes):
        res.append((meta.names[i], out[i]))
    res = sorted(res, key=lambda x: -x[1])
    return res

def detect(net, meta, image, thresh=.5, hier_thresh=.5, nms=.45):
    im = load_image(image, 0, 0)
    boxes = make_boxes(net)
    probs = make_probs(net)
    num =   num_boxes(net)
    network_detect(net, im, thresh, hier_thresh, nms, boxes, probs)
    res = []
    for j in range(num):
        for i in range(meta.classes):
            if probs[j][i] > 0:
                res.append((meta.names[i], probs[j][i], (boxes[j].x, boxes[j].y, boxes[j].w, boxes[j].h)))
    res = sorted(res, key=lambda x: -x[1])

    free_image(im)
    free_ptrs(cast(probs, POINTER(c_void_p)), num)
    return res

def darknet_polymage(frame, lib_func):
    rows = frame.shape[0]
    cols = frame.shape[1]
    total_pad = 56
    imwrite('messigray.png', frame)
    # IMAGE im;
    # im["w"] = frame[0]
    # im[""] = frame[1]
    # im->c = 3
    # im->data = ctypes.c_char_p()
    
    # res = detect(net, meta, 'messigray.png')
    test_detector("darknet/files/coco.data", net, '/home/e0358-3/asst-2/polymage/sandbox/video_demo/messigray.png', .5, .5, "try",0)
    # lib_func(ctypes.c_int(cols), \
    #          ctypes.c_int(rows), \
    #          net, \
    #          ctypes.c_char_p("/home/e0358-3/asst-2/darknet/cfg/coco.data"))
    # ret = np.empty((frame.shape[0], frame.shape[1], 3), np.uint8)
    ret = imread('try.jpg', IMREAD_COLOR)
    return ret

def add_darknet_app(app_id):

    # 1. modes
    modes = [ModeType.P_OPT]
    # 2. modes that need shared library
    lib_modes = modes

    # 3. set python functions from frame_process.py
    app_func_map = {}
    app_func_map[ModeType.P_OPT] = darknet_polymage

    # 4. create an App object
    app_dir = os.path.dirname(os.path.realpath(__file__))
    print app_dir
    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app
