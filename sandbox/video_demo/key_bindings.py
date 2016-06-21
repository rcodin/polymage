from structs import *

def generate_key_bindings(def_modes_map):
    QUIT = ModeType.QUIT
    PREV = ModeType.PREV
    CURRENT = ModeType.CURRENT
    CV2 = ModeType.CV2
    P_NAIVE = ModeType.P_NAIVE
    P_OPT = ModeType.P_OPT
    NUMBA = ModeType.NUMBA
    PIL = ModeType.PIL

    NONE = AppType.NONE
    UNSHARP = AppType.UNSHARP
    HARRIS = AppType.HARRIS
    BILATERAL = AppType.BILATERAL
    LAPLACIAN = AppType.LAPLACIAN

    # mode chars
    n = ord('n') # P_NAIVE / P_OPT
    q = ord('q') # QUIT
    o = ord(' ') # CV2
    p = ord('p') # PIL
    j = ord('j') # NUMBA (jit)

    mode_keys = {}
    mode_keys[QUIT] = q
    mode_keys[P_NAIVE] = n
    mode_keys[P_OPT] = n
    mode_keys[CV2] = o
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
    app_keys[NONE] = esc
    app_keys[HARRIS] = h
    app_keys[UNSHARP] = u
    app_keys[BILATERAL] = b
    app_keys[LAPLACIAN] = l

    if NONE not in def_modes_map:
        def_modes_map[NONE] = []

    # Map :
    # IN : ((APP, MODE), key)
    # OUT : (APP, MODE)
    key_bind = {}

    # toggle P_NAIVE and P_OPT mode
    for app_id in def_modes_map:
        key_bind[((app_id, P_NAIVE), n)] = (app_id, P_OPT)
        key_bind[((app_id, P_OPT), n)] = (app_id, P_NAIVE)

    # switch between apps
    print(def_modes_map)
    for app_id in def_modes_map:
        for mode_id in def_modes_map[app_id]:
            for app_id_k in def_modes_map:
                app_key = app_keys[app_id_k]
                key_bind[((app_id, mode_id), app_keys[app_id_k])] = (app_id_k, CURRENT)
            # Do nothing
            key_bind[((app_id, mode_id), app_keys[app_id])] = (NONE, CURRENT)

    toggle_modes = [CV2, NUMBA, PIL]
    # toggle to and from PREV mode
    for app_id in def_modes_map:
        t_modes = set(def_modes_map[app_id]).intersection(set(toggle_modes))
        for t_mode_id in t_modes:
            for mode_id in def_modes_map[app_id]:
                key_bind[((app_id, mode_id), mode_keys[t_mode_id])] = \
                    (app_id, t_mode_id)
            key_bind[((app_id, t_mode_id), mode_keys[t_mode_id])] = \
                (app_id, PREV)

    # Do nothing / Undefined:
    for app_id in def_modes_map:
        def_modes = def_modes_map[app_id]
        undef_modes = ModeType.modes_list.difference(set(def_modes))
        for d_mode_id in def_modes:
            for u_mode_id in undef_modes:
                key_bind[((app_id, d_mode_id), mode_keys[u_mode_id])] = \
                    (app_id, CURRENT)

    # QUIT mode
    for app_id in def_modes_map:
        for mode_id in def_modes_map[app_id]:
            key_bind[((app_id, mode_id), q)] = (app_id, QUIT)

    return key_bind

