from __init__ import *

import sys
import subprocess
import numpy as np
from fractions import Fraction

sys.path.insert(0, ROOT)

from compiler import *
from constructs import *

# PolyMage Specification
# ======================

def local_laplacian(pipe_data):
    L = 8 # number of pyramid levels
    levels = 8 # number of intensity levels

    # Input Parameters
    R = Parameter(Int, "R") # image rows
    C = Parameter(Int, "C") # image cols
    alpha = Parameter(Float, "alpha") # alpha
    beta = Parameter(Float, "beta") # beta
    # Images
    img_colour = Image(UShort, "img_colour", [3, R+4, C+4])

    pipe_data['R'] = R
    pipe_data['C'] = C

    # Vars
    x = Variable(Int, "x") # for rows
    y = Variable(Int, "y") # for cols
    c = Variable(Int, "c") # for colours
    r = Variable(Int, "r") # for remap LUT
    k = Variable(Int, "k") # for LUT levels

    # Intervals
    rgb = Interval(Int, 0, 2) # for colours
    rows = Interval(Int, 0, R+3)
    cols = Interval(Int, 0, C+3)
    rows2 = Interval(Int, 0, R-1)
    cols2 = Interval(Int, 0, C-1)
    lutLev = Interval(Int, 0, levels-1) # for LUT levels
    rLut = Interval(Int, -256*(levels-1), 256*(levels-1)-1) # for LUT


    #####################################################################################

    # DOWNSAMPLE
    def pyrDown(f, l, kl, name):
        # decrement factor
        decFactor = 1 << l
        # original factor
        orgFactor = 1 << (l-1)

        # domain (downsampled)
        decRowr = Interval(Int, 0, (R/decFactor)+3)
        decColr = Interval(Int, 0, (C/decFactor)+3)
        # domain (original)
        colr = Interval(Int, 0, (C/orgFactor)+3)

        # downsample in 'x' dimesnion (using a [1 3 3 1]' filter)

        # body case
        condx = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/decFactor)+1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/orgFactor)+1)
        # boundary cases
        condxLeft   = Condition(y, '<=', 1)
        condxRight  = Condition(y, '>=', (C/orgFactor)+2)
        condxTop    = Condition(x, '<=', 1)
        condxBottom = Condition(x, '>=', (R/decFactor)+2)
        
        # if the function does not have dimension 'k'
        if kl:
            downx = Function(([x, y], [decRowr, colr]), Float, "Dx_" + name + "_L" + str(l))
            downx.defn = [Case(condx, (f(2*(x-1)-1, y) +
                                                   f(2*(x-1)  , y) * 3.0 +
                                                   f(2*(x-1)+1, y) * 3.0 +
                                                   f(2*(x-1)+2, y)
                                                  ) * 0.125),
                          Case(condxLeft, 0),
                          Case(condxRight, 0),
                          Case(condxBottom, 0),
                          Case(condxTop, 0)]
        # if the function has dimension 'k'
        else:
            downx = Function(([k, x, y], [lutLev, decRowr, colr]), Float, "Dx_" + name + "_L" + str(l))
            downx.defn = [Case(condx, (f(k, 2*(x-1)-1, y) +
                                                   f(k, 2*(x-1)  , y) * 3.0 +
                                                   f(k, 2*(x-1)+1, y) * 3.0 +
                                                   f(k, 2*(x-1)+2, y)
                                                  ) * 0.125),
                          Case(condxLeft, 0),
                          Case(condxRight, 0),
                          Case(condxBottom, 0),
                          Case(condxTop, 0)]

        #fi

        # set the boundaries to zero


        # downsample in 'y' dimension (using a [1 3 3 1] filter)

        condy = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/decFactor)+1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/decFactor)+1)
        condyLeft   = Condition(y, '<=', 1)
        condyRight  = Condition(y, '>=', (C/decFactor)+2)
        condyTop    = Condition(x, '<=', 1)
        condyBottom = Condition(x, '>=', (R/decFactor)+2)


        if kl:
            downy = Function(([x, y], [decRowr, decColr]),Float, "D_" + name + "_L" + str(l))
            downy.defn = [Case(condy, (downx(x, 2*(y-1)-1) +
                                                   downx(x, 2*(y-1)  ) * 3.0 + 
                                                   downx(x, 2*(y-1)+1) * 3.0 +
                                                   downx(x, 2*(y-1)+2)
                                                  ) * 0.125),
                          Case(condyLeft, 0),
                          Case(condyRight, 0),
                          Case(condyBottom, 0),
                          Case(condyTop, 0)]
        else:
            downy = Function(([k, x, y], [lutLev, decRowr, decColr]),Float, "D_" + name + "_L" + str(l))
            downy.defn = [Case(condy, (downx(k, x, 2*(y-1)-1) +
                                                   downx(k, x, 2*(y-1)  ) * 3.0 +
                                                   downx(k, x, 2*(y-1)+1) * 3.0 +
                                                   downx(k, x, 2*(y-1)+2)
                                                  ) * 0.125),
                          Case(condyLeft, 0),
                          Case(condyRight, 0),
                          Case(condyBottom, 0),
                          Case(condyTop, 0)]
        #fi


        return downy


    # UPSAMPLE
    def pyrUp(f, l, kl, name):
        decFactor = 1 << l+1
        orgFactor = 1 << l

        # domain (original)
        decColr = Interval(Int, 0, (C/decFactor)+3)

        # domain (upsampled)
        colr = Interval(Int, 0, (C/orgFactor)+3)
        rowr = Interval(Int, 0, (R/orgFactor)+3)

        # upsample in 'x' dimension

        condx = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/orgFactor)+1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/decFactor)+1)
        condxLeft   = Condition(y, '<=', 1)
        condxRight  = Condition(y, '>=', (C/decFactor)+2)
        condxTop    = Condition(x, '<=', 1)
        condxBottom = Condition(x, '>=', (R/orgFactor)+2)


        if kl:
            upx = Function(([x, y], [rowr, decColr]),Float, "Ux_" + name + "_L" + str(l))
            evenXexpr = 0.25 * f(x/2  , y) + 0.75 * f(x/2+1, y)
            oddXexpr =  0.25 * f(x/2+2, y) + 0.75 * f(x/2+1, y)
        else:
            upx = Function(([k, x, y], [lutLev, rowr, decColr]),Float, "Ux_" + name + "_L" + str(l))
            evenXexpr = 0.25 * f(k, x/2  , y) + 0.75 * f(k, x/2+1, y)
            oddXexpr =  0.25 * f(k, x/2+2, y) + 0.75 * f(k, x/2+1, y)
        #fi
        upXexpr = Select(Condition(x%2, "==", 0),
                                evenXexpr,
                                oddXexpr)

        upx.defn = [Case(condx, upXexpr),
                    Case(condxLeft, 0),
                    Case(condxRight, 0),
                    Case(condxBottom, 0),
                    Case(condxTop, 0)]


        # upsample in 'y' dimension

        condy = Condition(x, '>=', 2) & \
                Condition(x, '<=', (R/orgFactor)+1) & \
                Condition(y, '>=', 2) & \
                Condition(y, '<=', (C/orgFactor)+1)
        condyLeft   = Condition(y, '<=', 1)
        condyRight  = Condition(y, '>=', (C/orgFactor)+2)
        condyTop    = Condition(x, '<=', 1)
        condyBottom = Condition(x, '>=', (R/orgFactor)+2)


        if kl:
            upy = Function(([x, y], [rowr, colr]),Float, "U_" + name + "_L" + str(l))
            evenYexpr = 0.25 * upx(x, y/2    ) + \
                        0.75 * upx(x, y/2 + 1)
            oddYexpr =  0.25 * upx(x, y/2 + 2) + \
                        0.75 * upx(x, y/2 + 1)
        else:
            upy = Function(([k, x, y], [lutLev, rowr, colr]),Float, "U_" + name + "_L" + str(l))
            evenYexpr = 0.25 * upx(k, x, y/2    ) + \
                        0.75 * upx(k, x, y/2 + 1)
            oddYexpr =  0.25 * upx(k, x, y/2 + 2) + \
                        0.75 * upx(k, x, y/2 + 1)
        #fi
        upYexpr = Select(Condition(y%2, "==", 0),
                                evenYexpr,
                                oddYexpr)

        upy.defn = [Case(condy, upYexpr),
                    Case(condyLeft, 0),
                    Case(condyRight, 0),
                    Case(condyBottom, 0),
                    Case(condyTop, 0)]

        return upy


    #####################################################################################
    # LOCAL LAPLACIAN :

    # 0. Convert to Gray
    img = Function(([x, y], [rows, cols]),Float, "img")
    img.defn = [(
                     0.299 * img_colour(2, x, y) +\
                     0.587 * img_colour(1, x, y) +\
                     0.114 * img_colour(0, x, y)
                     ) / 65535.0] #255.0]


    # 1. Make Gaussian Pyramid of the input.
    inGPyramid = {}
    inGPyramid[0] = img
    for l in range(1, L):
        inGPyramid[l] = pyrDown(inGPyramid[l-1], l, True, "inGPyramid")


    # 2. Remapping function
    '''
    # not as a LUT
    def remap(thing):
        fx = Cast(Float, thing) / 256.0
        val = alpha * fx * Exp( -fx*fx / 2.0 )
        return val
    '''

    # as a LUT
    fx = Cast(Float, r) / 256.0
    val = alpha * fx * Exp( -fx*fx / 2.0 ) # pow(e, _)

    remap = Function(([r], [rLut]),Float, "remapLUT")
    remap.defn = [val]


    # 3. Make the processed Gaussian Pyramid
    gPyramid = {}
    # Do a lookup into a lut with 256 entires per intensity level
    iLevel = k * (1.0 / (levels-1))
    idx = img(x, y) * Cast(Float, (levels-1)*256.0)
    # clamp idx
    idx = Min(
            Max(
              Cast(Int, idx),
              0),
            (levels-1)*256)

    gPyramid[0] = Function(([k, x, y], [lutLev, rows, cols]), Float, "gPyramid_L0")
    gPyramid[0].defn = [beta * ( img(x, y) - iLevel ) + 
                             iLevel + 
                             remap(idx - 256*k)]
    for l in range(1, L):
        gPyramid[l] = pyrDown(gPyramid[l-1], l, False, "gPyramid")


    # 4. Get the laplacian of processed Gaussian Pyramid
    lPyramid = {}
    # level: L-1
    lPyramid[L-1] = gPyramid[L-1]
    # level: everything else
    for l in range(L-2, -1, -1):
        rowr = Interval(Int, 0, (R/(1<<l))+3)
        colr = Interval(Int, 0, (C/(1<<l))+3)
        lPyramid[l] = Function(([k, x, y], [lutLev, rowr, colr]),Float, "lPyramid" + "_L"  + str(l))
        lPyramid[l].defn = [gPyramid[l](k, x, y) - 
                                 pyrUp(gPyramid[l+1], l, False, "gPyramid")(k, x, y)]


    # 5. Make the Laplacian Pyramid for the output
    outLPyramid = {}
    for l in range(0, L):
        # Split input pyramid value into integer(LUT index) and floating(weight) parts
        lev = inGPyramid[l](x, y) * Cast(Float, (levels-1))
        # LUT index - clamped:[0, levels-2]
        li = Min( Max( Cast(Int, lev), 0), levels-2 )
        # weight
        lf = lev - li

        rowr = Interval(Int, 0, (R/(1<<l))+3)
        colr = Interval(Int, 0, (C/(1<<l))+3)

        # Linearly interpolate between the nearest processed pyramid levels
        outLPyramid[l] = Function(([x, y], [rowr, colr]),Float, "outLPyramid_L" + str(l))
        outLPyramid[l].defn = [lPyramid[l](li  , x, y) * (1.0 - lf) + 
                                    lPyramid[l](li+1, x, y) * lf]


    # 6. Make the Gaussian pyramid of the output
    outGPyramid = {}
    outGPyramid[L-1] = outLPyramid[L-1]
    for l in range(L-2, -1, -1):
        if l == 0:
            outGPyramid[l] = Function(([x, y], [rows, cols]),Float, "result_ref_gray")
        else:
            rowr = Interval(Int, 0, (R/(1<<l))+3)
            colr = Interval(Int, 0, (C/(1<<l))+3)
            outGPyramid[l] = Function(([x, y], [rowr, colr]),Float, "outGPyramid_L" + str(l))
        #fi
        outGPyramid[l].defn = [pyrUp(outGPyramid[l+1], l, True, "outGPyramid")(x, y) + 
                                    outLPyramid[l](x, y)]

    result_ref_gray = outGPyramid[0]


    # 7. 
    # Halide :
    # "Reintroduce color (Connelly: use eps to avoid scaling up noise
    # w/ apollo3.png input)"
    eps = 0.01
    colourMap = outGPyramid[0](x+2, y+2) * \
                ((img_colour(c, x+2, y+2)/65535.0) + eps) / (img(x+2, y+2) + eps)
    colourMap = Min( Max( colourMap, 0.0 ), 1.0 ) # clamp
    colourMap = Cast( UShort, colourMap * 65535.0)

    colour = Function(([c, x, y], [rgb, rows2, cols2]),UShort, "laplacian")
    colour.defn = [colourMap]

    result_ref_colour = colour
    return result_ref_colour

    #####################################################################################
# END
