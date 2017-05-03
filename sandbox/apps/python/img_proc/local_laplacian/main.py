import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_laplacian
from exec_pipe import laplacian
#from app_tuner import auto_tune

app = "laplacian"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune':
        print("Tuning...")
        #auto_tune(app_data)
    else:
        create_lib(build_laplacian, app, app_data)
        for i in range (0, 5):
            laplacian(app_data)

    return

main()
