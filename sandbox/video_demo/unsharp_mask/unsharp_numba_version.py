import numpy as np
#from numpy import jit

# unsharp mask parameters
thresh = 0.001
weight = 3

#@jit(nogil = True, cache = True, nopython = True)
def mini(x, y):
    return x if x < y else y

#@jit(nogil = True, cache = True, nopython = True)
def maxi(x, y):
    return x if x > y else y

#@jit("float32[::](uint8[::], int64)", nogil = True, cache = True)
def unsharp_numba(frame, lib_func):
    r = frame.shape[0] - 4
    c = frame.shape[1] - 4
    image_f = np.float32(frame) / 255.0
    res = image_f
    T_x = 32
    T_y = 256
    im = np.zeros((3, (r+4), (c+4)), np.float32)

    # TODO: copy numpy style
    for i in range(r+4):
        for j in range(c+4):
            for k in range(3):
                im[k,i,j] = image_f[i,j,k]

    blurx = np.zeros((3, T_x, T_y + 10), dtype=np.float32)
    blury = np.zeros((3, T_x, T_y + 10), dtype=np.float32)
    sharpen = np.zeros((3, T_x, T_y + 10), dtype=np.float32)

    for ti1 in range((int(r + 1) / T_x) + 1):
        ct0 = mini(r + 1, (T_x * ti1) + T_x - 1)
        ct1 = maxi(2, (T_x * ti1))
        ct4 = mini(r + 1, (T_x * ti1) + T_x - 1)
        ct5 = maxi(2, (T_x * ti1))
        ct8 = mini(r + 1, (T_x * ti1) + T_x - 1)
        ct9 = maxi(2, T_x * ti1)
        ct12 = mini(r + 1, (T_x * ti1) + T_x - 1)
        ct13 = maxi(2, T_x * ti1)
        for ti2 in range(-1, (int(c + 3) / T_y) + 1):
            ct2 = mini(c + 3, (T_y*ti2) + T_y + 5)
            ct3 = maxi(0, (T_y * ti2))
            ct6 = mini(c + 1, (T_y * ti2) + T_y + 4)
            ct7 = maxi(2, (T_y * ti2) + 1)
            ct10 = mini(c + 1, (T_y * ti2) + T_y + 3)
            ct11 = maxi(2, (T_y * ti2) + 2)
            ct14 = mini(c + 1, (T_y * ti2) + T_y + 2)
            ct15 = maxi(2, (T_y * ti2) + 3)

            # compute 'blurx'
            for i0 in range(3):
                for i1 in range(ct1, ct0 + 1):
                    for i2 in range(ct3, ct2 + 1):
                        blurx[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2] = \
                            im[i0, i1 - 2, i2] * 0.0625 + \
                            im[i0, i1 - 1, i2] * 0.25 + \
                            im[i0 ,i1 + 0, i2] * 0.375 + \
                            im[i0, i1 + 1, i2] * 0.25 + \
                            im[i0, i1 + 2, i2] * 0.0625

            # compute blury
            for i0 in range(3):
                for i1 in range(ct5, ct4 + 1):
                    for i2 in range(ct7, ct6 + 1):
                        blury[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2] = \
                            blurx[i0, (-T_x*ti1) + i1, -2 + (-T_y*ti2) + i2] * 0.0625 + \
                            blurx[i0, (-T_x*ti1) + i1, -1 + (-T_y*ti2) + i2] * 0.25 + \
                            blurx[i0, (-T_x*ti1) + i1,  0 + (-T_y*ti2) + i2] * 0.375 + \
                            blurx[i0, (-T_x*ti1) + i1,  1 + (-T_y*ti2) + i2] * 0.25 + \
                            blurx[i0, (-T_x*ti1) + i1,  2 + (-T_y*ti2) + i2] * 0.0625

            # compute 'sharp'
            for i0 in range(3):
                for i1 in range(ct9, ct8 + 1):
                    for i2 in range(ct11, ct10 + 1):
                        sharpen[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2] = \
                            im[i0, i1, i2] * (1 + weight) - \
                            blury[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2] * weight

            # compute final 'res'
            for i0 in range(3):
                for i1 in range(ct13, ct12 + 1):
                    for i2 in range(ct15, ct14 + 1):
                        res[i1, i2, i0] = im[i0, i1, i2] \
                            if im[i0, i1, i2] - \
                                blury[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2] < thresh \
                        else sharpen[i0, (-T_x*ti1) + i1, (-T_y*ti2) + i2]

    return res
