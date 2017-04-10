import numpy as np
import time
import sys

from __init__ import *

from init import init_all
from printer import print_header, print_config, print_line
from builder import create_lib,build_harris
from exec_pipe import harrispipe
from app_tuner import auto_tune

app = "harris"

def main():
    print_header()

    app_data = {}
    app_data['app'] = app
    app_data['ROOT'] = ROOT

    init_all(app_data)
    print_config(app_data)

    if app_data['mode'] == 'tune+':
        for g_size in [3, 5, 7, 200]:
            create_lib(build_harris, app, app_data, g_size)
            for t in range (0, 5):
                print ("Running for iteration #", t)
                harrispipe(app_data)
    elif app_data['mode'] == 'tune':
        pass
    else:
        create_lib(build_harris, app, app_data)
        harrispipe(app_data)

    return

main()
