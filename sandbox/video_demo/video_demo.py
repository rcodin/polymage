import sys
from common import clock, draw_str
from structs import *
from key_bindings import *

sys.path.insert(0, "harris_corner")
sys.path.insert(0, "unsharp_mask")
sys.path.insert(0, "bilateral_grid")
sys.path.insert(0, "local_laplacian")

from harris_process import *
from unsharp_process import *
from bilateral_process import *
from laplacian_process import *

def none_app_process(frame, lib_func):
    return frame

def add_none_app(app_id):
    modes = [ModeType.NONE]
    lib_modes = []
    app_dir = None
    app_func_map = {}
    app_func_map[modes[0]] = none_app_process

    app = App(app_id, app_dir, modes, lib_modes, app_func_map)

    return app

def app_init():
    # global map for all apps
    app_map = {}
    app_modes_map = {}
    apps_list = AppType.apps_list

    ''' Add Unhsarp Mask app '''
    app_id = AppType.UNSHARP
    assert app_id in apps_list
    app_map[app_id] = add_unsharp_app(app_id)

    ''' Add Harris Corner app '''
    app_id = AppType.HARRIS
    assert app_id in apps_list
    app_map[app_id] = add_harris_app(app_id)

    ''' Add Bilateral Grid app '''
    app_id = AppType.BILATERAL
    assert app_id in apps_list
    app_map[app_id] = add_bilateral_app(app_id)

    ''' Add Local Laplacian app '''
    app_id = AppType.LAPLACIAN
    assert app_id in apps_list
    app_map[app_id] = add_laplacian_app(app_id)

    ''' Add a NONE type app '''
    app_id = AppType.NONE
    assert app_id in apps_list
    app_map[app_id] = add_none_app(app_id)

    for app_id in app_map.keys():
        app_modes_map[app_id] = app_map[app_id].modes.keys()

    return app_map, app_modes_map

def app_destroy(apps_map):
    for app_id in apps_map:
        apps_map[app_id].destroy()
    return

def main():
    assert len(sys.argv) >= 2, "Video file path not specifed"
    video_file = sys.argv[1]

    # initialize all apps
    apps_map, app_modes_map = app_init()

    # key bindings for video mode switching
    key_bindings = generate_key_bindings(app_modes_map)

    # create a video processor object
    video_proc = VideoProcessor(apps_map, video_file, key_bindings)

    # start videp processing
    video_proc.process()

    # terminate
    video_proc.finish()

    # destroy lib apps and memory pool
    app_destroy(apps_map)

main()
