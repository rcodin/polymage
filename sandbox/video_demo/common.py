#!/usr/bin/env python

import numpy as np
import cv2
import os

def draw_str(dst, pixel, s):
    x = pixel[0]
    y = pixel[1]
    cv2.putText(dst, s, (x+1, y+1), cv2.FONT_HERSHEY_PLAIN, 2.5, (0, 0, 0), thickness = 3)
    cv2.putText(dst, s, (x, y), cv2.FONT_HERSHEY_PLAIN, 2.5, (0, 0, 255), thickness = 3)

def clock():
    return cv2.getTickCount() / cv2.getTickFrequency()

def image_clamp(image_in, image_out, \
          R, C, K, \
          dtype, dfactor, \
          left, total):
    if K > 1:
        # mid of top, bottom, left and right resp.
        image_out[0:left, left:C+left, 0:K] = \
            np.array(image_in[0, 0:C, 0:K] * dfactor, dtype)
        image_out[R+left:R+total, left:C+left, 0:K] = \
            np.array(image_in[R-1, 0:C, 0:K] * dfactor, dtype)
        image_out[left:left+R, 0:left, 0:K] = \
            np.array(image_in[0:R, 0, 0:K].reshape(R, 1, 3) * dfactor, dtype)
        image_out[left:left+R, left+C:C+total, 0:K] = \
            np.array(image_in[0:R, C-1, 0:K].reshape(R, 1, 3) * dfactor, dtype)
        # corners :
        image_out[0:left, 0:left, 0:K] = \
            image_out[left, 0:left, 0:K]
        image_out[0:left, left+C:C+total, 0:K] = \
            image_out[left, left+C:C+total, 0:K]
        image_out[left+R:R+total, 0:left, 0:K] = \
            image_out[left+R-1, 0:left, 0:K]
        image_out[left+R:R+total, left+C:C+total, 0:K] = \
            image_out[left+R-1, left+C:C+total, 0:K]
    else:
        # mid of top, bottom, left and right resp.
        image_out[0:left, left:C+left] = \
            np.array(image_in[0, 0:C] * dfactor, dtype)
        image_out[R+left:R+total, left:C+left] = \
            np.array(image_in[R-1, 0:C] * dfactor, dtype)
        image_out[left:left+R, 0:left] = \
            np.array(image_in[0:R, 0].reshape(R, 1, 3) * dfactor, dtype)
        image_out[left:left+R, left+C:C+total] = \
            np.array(image_in[0:R, C-1].reshape(R, 1, 3) * dfactor, dtype)
        # corners :
        image_out[0:left, 0:left] = \
            image_out[left, 0:left]
        image_out[0:left, left+C:C+total] = \
            image_out[left, left+C:C+total]
        image_out[left+R:R+total, 0:left] = \
            image_out[left+R-1, 0:left]
        image_out[left+R:R+total, left+C:C+total] = \
            image_out[left+R-1, left+C:C+total]
