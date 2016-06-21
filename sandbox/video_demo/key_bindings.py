from structs import *

def generate_key_bindings(def_modes_map):
    QUIT = ModeType.QUIT
    PREV = ModeType.PREV
    NONE_MODE = ModeType.NONE
    CURRENT = ModeType.CURRENT
    CV2 = ModeType.CV2
    P_NAIVE = ModeType.P_NAIVE
    P_OPT = ModeType.P_OPT
    NUMBA = ModeType.NUMBA
    PIL = ModeType.PIL

    NONE_APP = AppType.NONE
    UNSHARP = AppType.UNSHARP
    HARRIS = AppType.HARRIS
    BILATERAL = AppType.BILATERAL
    LAPLACIAN = AppType.LAPLACIAN

    # mode chars
    none = 255
    q = ord('q') # QUIT
    n = ord('n') # P_NAIVE
    o = ord('o') # P_OPT
    c = ord(' ') # CV2
    p = ord('p') # PIL
    j = ord('j') # NUMBA (jit)

    mode_keys = {}
    mode_keys[QUIT] = q
    mode_keys[NONE_MODE] = none
    mode_keys[P_NAIVE] = n
    mode_keys[P_OPT] = o
    mode_keys[CV2] = c
    mode_keys[PIL] = p
    mode_keys[NUMBA] = j

    # app chars
    esc = ord('\x1b')
    u = ord('u')
    h = ord('h')
    b = ord('b')
    l = ord('l')

    app_keys = {}
    # TODO : add key for NONE : ESC
    app_keys[NONE_APP] = esc
    app_keys[HARRIS] = h
    app_keys[UNSHARP] = u
    app_keys[BILATERAL] = b
    app_keys[LAPLACIAN] = l

    if NONE_APP not in def_modes_map:
        def_modes_map[NONE_APP] = []

    # Map :
    # IN : ((APP, MODE), key)
    # OUT : (APP, MODE)
    key_bind = {}

    # switch between apps
    # all to all : irrespective of the current mode
    for app_id in def_modes_map:
        for mode_id in def_modes_map[app_id]:
            for app_id_k in def_modes_map:
                app_key = app_keys[app_id_k]
                key_bind[((app_id, mode_id), app_keys[app_id_k])] = (app_id_k, CURRENT)
            # Do nothing
            key_bind[((app_id, mode_id), app_keys[app_id])] = (NONE_APP, CURRENT)

    # switch between modes in an app
    # toggle to and from PREV mode
    for app_id in def_modes_map:
        modes = def_modes_map[app_id]
        # all to all
        for k_mode_id in modes:
            for v_mode_id in modes:
                key_bind[((app_id, k_mode_id), mode_keys[v_mode_id])] = \
                    (app_id, v_mode_id)
            # switch from mode X to itself results in a toggle
            key_bind[((app_id, k_mode_id), mode_keys[k_mode_id])] = \
                (app_id, PREV)

    # Do nothing / Undefined:
    for app_id in def_modes_map:
        def_modes = def_modes_map[app_id]
        undef_modes = ModeType.modes_list.difference(set(def_modes))
        for d_mode_id in def_modes:
            for u_mode_id in undef_modes:
                key_bind[((app_id, d_mode_id), mode_keys[u_mode_id])] = (app_id, CURRENT)

    # QUIT mode and NONE mode
    for app_id in def_modes_map:
        for mode_id in def_modes_map[app_id]:
            key_bind[((app_id, mode_id), q)] = (app_id, QUIT)
            key_bind[((app_id, mode_id), none)] = (app_id, CURRENT)

    all_keys = app_keys.values()
    all_keys.extend(mode_keys.values())
    key_bind[0] = all_keys

    return key_bind

