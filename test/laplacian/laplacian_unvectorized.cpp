#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include "simple_pool_allocator.h"
#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))
extern "C" void pipeline_laplacian(int C, int R, float alpha, float beta, void * img_colour_void_arg, void * laplacian_void_arg)
{
  short unsigned int * img_colour;
  img_colour = (short unsigned int *) (img_colour_void_arg);
  short unsigned int * laplacian;
  laplacian = (short unsigned int *) (laplacian_void_arg);
  /* users : ['remapLUT'] */
  float * _arr_4_0;
  _arr_4_0 = (float *) (pool_allocate((sizeof(float) * 3584)));
  for (int _i0 = -1792; (_i0 <= 1791); _i0 = (_i0 + 1))
  {
    _arr_4_0[(_i0 - -1792)] = ((alpha * ((float) (_i0) / 256.0)) * std::exp(((-(((float) (_i0) / 256.0)) * ((float) (_i0) / 256.0)) / 2.0)));
  }
  #pragma omp parallel for schedule(static) collapse(2)
  for (int _T_i1 = 0; (_T_i1 <= 2563); _T_i1 = (_T_i1 + 1))
  {
    for (int _T_i2 = 0; (_T_i2 <= 6); _T_i2 = (_T_i2 + 1))
    {
      /* users : ['img'] */
      float _buf_1_1[(1 * 256)];
      /* users : ['kgPyramid_L0'] */
      float _buf_2_2[(8 * 256)];
      /* users : ['outLPyramid_L0'] */
      float _buf_1_3[(1 * 256)];
      int _ct0 = ((1539 < ((256 * _T_i2) + 255))? 1539: ((256 * _T_i2) + 255));
      #pragma ivdep
      for (int _i2 = (256 * _T_i2); (_i2 <= _ct0); _i2 = (_i2 + 1))
      {
        _buf_1_1[(0 + ((-256 * _T_i2) + _i2))] = ((((0.299 * img_colour[(((0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)]) + (0.587 * img_colour[(((1 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)])) + (0.114 * img_colour[(((2 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)])) / 65535.0);
      }
      for (int _i0 = 0; (_i0 <= 7); _i0 = (_i0 + 1))
      {
        int _ct1 = ((1539 < ((256 * _T_i2) + 255))? 1539: ((256 * _T_i2) + 255));
        #pragma ivdep
        for (int _i2 = (256 * _T_i2); (_i2 <= _ct1); _i2 = (_i2 + 1))
        {
          int _ct2 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (1792.0)));
          int _ct3 = 0;
          int _ct4 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (1792.0))) > 0)? _ct2: _ct3);
          int _ct5 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (1792.0)));
          int _ct6 = 0;
          int _ct7 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (1792.0))) > 0)? _ct5: _ct6);
          int _ct8 = _ct7;
          int _ct9 = 1792;
          int _ct10 = ((_ct4 < 1792)? _ct8: _ct9);
          _buf_2_2[((_i0 * (1 * 256)) + ((-256 * _T_i2) + _i2))] = (((beta * (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] - (_i0 * 0.14285714285714285))) + (_i0 * 0.14285714285714285)) + _arr_4_0[((_ct10 - (256 * _i0)) - -1792)]);
        }
      }
      int _ct11 = ((1539 < ((256 * _T_i2) + 255))? 1539: ((256 * _T_i2) + 255));
      #pragma ivdep
      for (int _i2 = (256 * _T_i2); (_i2 <= _ct11); _i2 = (_i2 + 1))
      {
        int _ct12 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct13 = 0;
        int _ct14 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct12: _ct13);
        int _ct15 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct16 = 0;
        int _ct17 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct15: _ct16);
        int _ct18 = _ct17;
        int _ct19 = 6;
        int _ct20 = ((_ct14 < 6)? _ct18: _ct19);
        int _ct21 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct22 = 0;
        int _ct23 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct21: _ct22);
        int _ct24 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct25 = 0;
        int _ct26 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct24: _ct25);
        int _ct27 = _ct26;
        int _ct28 = 6;
        int _ct29 = ((_ct23 < 6)? _ct27: _ct28);
        int _ct30 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct31 = 0;
        int _ct32 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct30: _ct31);
        int _ct33 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct34 = 0;
        int _ct35 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct33: _ct34);
        int _ct36 = _ct35;
        int _ct37 = 6;
        int _ct38 = ((_ct32 < 6)? _ct36: _ct37);
        int _ct39 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct40 = 0;
        int _ct41 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct39: _ct40);
        int _ct42 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)));
        int _ct43 = 0;
        int _ct44 = (((int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7))) > 0)? _ct42: _ct43);
        int _ct45 = _ct44;
        int _ct46 = 6;
        int _ct47 = ((_ct41 < 6)? _ct45: _ct46);
        _buf_1_3[(0 + ((-256 * _T_i2) + _i2))] = ((_buf_2_2[((_ct20 * (1 * 256)) + ((-256 * _T_i2) + _i2))] * (1.0 - ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)) - _ct29))) + (_buf_2_2[(((_ct38 + 1) * (1 * 256)) + ((-256 * _T_i2) + _i2))] * ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)) - _ct47)));
      }
      if (((_T_i1 <= 2559) && (_T_i2 <= 5)))
      {
        for (int _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
        {
          #pragma ivdep
          for (int _i2 = (256 * _T_i2); (_i2 <= ((256 * _T_i2) + 255)); _i2 = (_i2 + 1))
          {
            float _ct48 = ((_buf_1_3[(0 + ((-256 * _T_i2) + _i2))] * ((img_colour[(((_i0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)] / 65535.0) + 0.01)) / (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] + 0.01));
            float _ct49 = 0.0;
            float _ct50 = ((((_buf_1_3[(0 + ((-256 * _T_i2) + _i2))] * ((img_colour[(((_i0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)] / 65535.0) + 0.01)) / (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] + 0.01)) > 0.0)? _ct48: _ct49);
            float _ct51 = ((_buf_1_3[(0 + ((-256 * _T_i2) + _i2))] * ((img_colour[(((_i0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)] / 65535.0) + 0.01)) / (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] + 0.01));
            float _ct52 = 0.0;
            float _ct53 = ((((_buf_1_3[(0 + ((-256 * _T_i2) + _i2))] * ((img_colour[(((_i0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)] / 65535.0) + 0.01)) / (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] + 0.01)) > 0.0)? _ct51: _ct52);
            float _ct54 = _ct53;
            float _ct55 = 1.0;
            float _ct56 = ((_ct50 < 1.0)? _ct54: _ct55);
            laplacian[(((_i0 * (R * C)) + (_T_i1 * C)) + _i2)] = (unsigned short) ((_ct56 * 65535.0));
          }
        }
      }
    }
  }
  pool_deallocate(_arr_4_0);
}
