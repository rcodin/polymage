import sys
import os.path
from PIL import Image
import numpy as np
from arg_parser import parse_args

from printer import print_header, print_usage, print_line

def init_images(app_data):
    print("[init.py] : initializing images...")

    app_args = app_data['app_args']

    # input image: 
    img_path = app_args.img_file
    img = np.array(Image.open(img_path))
    rows, cols, c = img.shape

    alpha = float(app_args.alpha)
    beta = float(app_args.beta)
    levels = 8
    L = 8

    if rows != 2560 or cols != 1536:
        print("Please use 1536x2560 image size")
        sys.exit(0)

    alpha /= (levels-1)

    image_ghost = np.zeros((rows+4, cols+4, 3), np.uint16)
    image_ghost[2:rows+2, 2:cols+2, 0:3] = \
        np.array(img * 256, np.uint16)

    # move colour dimension outside
    image_flip = np.rollaxis(image_ghost, 2)
    image_flip = image_flip.ravel()

    # result array
    OUT = np.zeros((3, rows, cols), np.uint16).ravel()


    # final output image
    #OUT = np.zeros((rows, cols), np.float32).ravel()

    img_data = {}
    img_data['IN'] = image_flip
    img_data['OUT'] = OUT

    app_data['img_data'] = img_data
    app_data['rows'] = rows
    app_data['cols'] = cols
    app_data['alpha'] = alpha
    app_data['beta'] = beta
    return

def get_input(app_data):
    # parse the command-line arguments
    app_args = parse_args()
    app_data['app_args'] = app_args

    app_data['mode'] = app_args.mode
    app_data['runs'] = int(app_args.runs)
    app_data['graph_gen'] = bool(app_args.graph_gen)
    app_data['timer'] = app_args.timer

    # storage optimization
    app_data['optimize_storage'] = bool(app_args.optimize_storage)
    # early freeing of allocated arrays
    app_data['early_free'] = bool(app_args.early_free)
    # pool allocate option
    app_data['pool_alloc'] = bool(app_args.pool_alloc)

    return

def init_all(app_data):
    pipe_data = {}
    app_data['pipe_data'] = pipe_data

    get_input(app_data)

    init_images(app_data)

    return

