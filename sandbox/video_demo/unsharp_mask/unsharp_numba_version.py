import numpy as np
from numba import jit

thresh = 0.001
weight = 3

@jit("float32[::](uint8[::], none)", nogil = True, cache = True)
def unsharp_numba(frame, lib_func):
    image_f = np.float32(frame) / 255.0
    r = frame.shape[0] - 4
    c = frame.shape[1] - 4
    im = np.zeros((3 * (r+4) * (c+4)),dtype = np.float32)
    res = image_f
    #im = np.rollaxis(image_f,2).flatten()
    #Loops seem to be faster than Numpy functions when optimized with Numba
    for i in range(r+4):
        for j in range(c+4):
            for k in range(3):
                im[k*(r+4)*(c+4) + i*(c+4) + j] = image_f[i,j,k]
    T_x = 32
    T_y = 256

    for ti1 in range(((r + 1) / T_x) + 1):
        blurx = np.zeros((3, T_x + 10, T_y + 10), dtype = np.float32)
        blury = np.zeros((3, T_x + 10, T_y + 10), dtype = np.float32)
        sharpen = np.zeros((3, T_x + 10, T_y + 10), dtype = np.float32)
        ct0 = min(r + 1, (T_x * ti1) + T_x - 1)
        ct1 = max(2, (T_x * ti1))
        ct4 = min(r + 1, (T_x * ti1) + T_x - 1)
        ct5 = max(2, (T_x * ti1))
        ct8 = min(r + 1, (T_x * ti1) + T_x - 1)
        ct9 = max(2, T_x * ti1)
        ct12 = min(r + 1, (T_x * ti1) + T_x - 1)
        ct13 = max(2, T_x * ti1)
        for ti2 in range(-1, ((c + 3) / T_y) + 1):
            ct2 = min(c + 3, (T_y * ti2) + T_y + 5)
            ct3 = max(0, (T_y * ti2))
            ct6 = min(c + 1, (T_y * ti2) + T_y + 4)
            ct7 = max(2, (T_y * ti2) + 1)
            ct10 = min(c + 1, (T_y * ti2) + T_y + 3)
            ct11 = max(2, (T_y * ti2) + 2)
            ct14 = min(c + 1, (T_y * ti2) + T_y + 2)
            ct15 = max(2, (T_y * ti2) + 3)

            # compute 'blurx'
            for i0 in range(3):
                for i1 in range(ct1, ct0 + 1):
                    for i2 in range(ct3, ct2 + 1):
                        blurx[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2] = \
                            im[i0 * (r+4) * (c+4) +  (i1 - 2) * (c+4) + i2] * 0.0625 + \
                            im[i0 * (r+4) * (c+4) +  (i1 - 1) * (c+4) + i2] * 0.25 + \
                            im[i0 * (r+4) * (c+4) +  (i1 - 0) * (c+4) + i2] * 0.375 + \
                            im[i0 * (r+4) * (c+4) +  (i1 + 1) * (c+4) + i2] * 0.25 + \
                            im[i0 * (r+4) * (c+4) +  (i1 + 2) * (c+4) + i2] * 0.0625

            # compute 'blury'
            for i0 in range(3):
                for i1 in range(ct5, ct4 + 1):
                    for i2 in range(ct7, ct6 + 1):
                        blury[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2] = \
                            blurx[i0, (-T_x * ti1) + i1, -2 + (-T_y * ti2) + i2] * 0.0625 + \
                            blurx[i0, (-T_x * ti1) + i1, -1 + (-T_y * ti2) + i2] * 0.25 + \
                            blurx[i0, (-T_x * ti1) + i1,  0 + (-T_y * ti2) + i2] * 0.375 + \
                            blurx[i0, (-T_x * ti1) + i1,  1 + (-T_y * ti2) + i2] * 0.25 + \
                            blurx[i0, (-T_x * ti1) + i1,  2 + (-T_y * ti2) + i2] * 0.0625

            # compute 'sharp'
            for i0 in range(3):
                for i1 in range(ct9, ct8 + 1):
                    for i2 in range(ct11, ct10 + 1):
                        sharpen[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2] = \
                        im[i0 * (r+4) * (c+4) +  i1 * (c+4) + i2] * (1 + weight) \
                        - blury[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2] * weight

            # compute final 'res'
            for i0 in range(3):
                for i1 in range(ct13, ct12 + 1):
                    for i2 in range(ct15, ct14 + 1):
                        if abs(im[i0 * (r+4) * (c+4) +  i1 * (c+4) + i2] - blury[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2]) >= thresh:
                            res[i1,i2,i0] = sharpen[i0, (-T_x * ti1) + i1, (-T_y * ti2) + i2]
                        



    return res
