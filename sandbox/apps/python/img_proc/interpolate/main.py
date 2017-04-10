import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_interpolate
from exec_pipe import interpolate
#from app_tuner import auto_tune

app = "multiscale_interpolate"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune+':
        for g_size in [3, 5, 7, 10, 15, 20, 30, 200]:
            create_lib(build_interpolate, app, app_data, g_size)
            for t in range (0, 5):
                print ("Running for iteration #", t)
                interpolate(app_data)
    elif app_data['mode'] == 'tune':
        pass
    else:
        create_lib(build_interpolate, app, app_data)
        for t in range (0, 10):
                print ("Running for iteration #", t)
                interpolate(app_data)

    return

main()
