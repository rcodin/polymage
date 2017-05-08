#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <cmath>
#include <string.h>
#include <immintrin.h>
#include <mmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <tmmintrin.h>
#include <smmintrin.h>
#include <nmmintrin.h>
//#include <ammintrin.h>
#include <x86intrin.h>

#include "simple_pool_allocator.h"

#define isl_min(x,y) ((x) < (y) ? (x) : (y))
#define isl_max(x,y) ((x) > (y) ? (x) : (y))
#define isl_floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))

#define vSIZE 8
typedef unsigned short int v8si __attribute__ ((vector_size(sizeof(short int)*8)));
__m256 cv65535 = {65535.0f, 65535.0f, 65535.0f, 65535.0f, 65535.0f, 65535.0f, 65535.0f, 65535.0f};
__m256 const1792 = {1792.0, 1792.0, 1792.0, 1792.0, 1792.0, 1792.0, 1792.0, 1792.0};
__m256 const7f = {7.0f, 7.0f, 7.0f, 7.0f, 7.0f, 7.0f, 7.0f, 7.0f};
__m256i _cv0i = {0};

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

      int _i1 = _T_i1;
      if (_ct0 == 256 * _T_i2 + 255)
        {
            __m256 const299 = {0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299, 0.299};
            __m256 const587 = {0.587, 0.587, 0.587, 0.587, 0.587, 0.587, 0.587, 0.587};
            __m256 const114 = {0.114, 0.114, 0.114, 0.114, 0.114, 0.114, 0.114, 0.114};
            
            for (int _i2 = (256 * _T_i2); (_i2 <= _ct0); _i2 = (_i2 + 8))
            {
                __m128i img_vec0 = _mm_loadu_si128 ((__m128i *)&img_colour[(((0 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)]);
                __m128i xlo = _mm_unpacklo_epi16(img_vec0, _mm_set1_epi16(0));
                __m128i xhi = _mm_unpackhi_epi16(img_vec0, _mm_set1_epi16(0));
                __m128 ylo = _mm_cvtepi32_ps(xlo);
                __m128 yhi = _mm_cvtepi32_ps(xhi);
                __m256 img_vec_f0 = _mm256_set_m128 (yhi, ylo);
                
                __m128i img_vec1 = _mm_loadu_si128 ((__m128i *)&img_colour[(((1 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)]);
                xlo = _mm_unpacklo_epi16(img_vec0, _mm_set1_epi16(0));
                xhi = _mm_unpackhi_epi16(img_vec0, _mm_set1_epi16(0));
                ylo = _mm_cvtepi32_ps(xlo);
                yhi = _mm_cvtepi32_ps(xhi);
                __m256 img_vec_f1 = _mm256_set_m128 (yhi, ylo);
                
                __m128i img_vec2 = _mm_loadu_si128 ((__m128i *)&img_colour[(((2 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)]);
                xlo = _mm_unpacklo_epi16(img_vec0, _mm_set1_epi16(0));
                xhi = _mm_unpackhi_epi16(img_vec0, _mm_set1_epi16(0));
                ylo = _mm_cvtepi32_ps(xlo);
                yhi = _mm_cvtepi32_ps(xhi);
                __m256 img_vec_f2 = _mm256_set_m128 (yhi, ylo);
                
                __m256 t1 = _mm256_add_ps (_mm256_mul_ps (img_vec_f0, const299), 
                                           _mm256_add_ps (_mm256_mul_ps (img_vec_f1, const587), 
                                                          _mm256_mul_ps (img_vec_f2, const114)));
                
                __m256 t2 = _mm256_div_ps (t1, cv65535);
                _mm256_storeu_ps (&_buf_1_1[(0 + ((-256 * _T_i2) + _i2))], t2);
            }
        }
        else
        {
            #pragma ivdep
            for (int _i2 = (256 * _T_i2); (_i2 <= _ct0); _i2 = (_i2 + 1))
            {
                _buf_1_1[(0 + ((-256 * _T_i2) + _i2))] = 
                ((((0.299 * img_colour[(((0 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)]) + 
                (0.587 * img_colour[(((1 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)])) + 
                (0.114 * img_colour[(((2 * ((4 + R) * (4 + C))) + (_i1 * (4 + C))) + _i2)])) / 65535.0);
            }
        }
        
      for (int _i0 = 0; (_i0 <= 7); _i0 = (_i0 + 1))
      {
        int _ct1 = ((1539 < ((256 * _T_i2) + 255))? 1539: ((256 * _T_i2) + 255));
        int _ct3 = _ct1;
        
        if (_ct3 == (256*_T_i2 + 255))
          {
              __m256i i0_vec = _mm256_set1_epi32 (256*_i0);
              __m256i const1792i = _mm256_set1_epi32 (1792.0);
              
              for (int _i2 = (256 * _T_i2); (_i2 <= _ct3); _i2 = (_i2 + 8))
              {
                int _ct4 = (int) ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (1792.0)));
                int _ct6 = ((_ct4 > 0)? _ct4: 0);
                int _ct12 = ((_ct6 < 1792)? _ct6: 1792);
                __m256 _buf_1_1_vec = _mm256_loadu_ps(&_buf_1_1[(0 + ((-256 * _T_i2) + _i2))]);
                __m256i _ct4_vec = _mm256_cvttps_epi32 (_mm256_mul_ps(_buf_1_1_vec, const1792));
                __m256i _ct6_vec =  _mm256_max_epu32 (_ct4_vec, _cv0i);
                __m256i _ct12_vec = _mm256_min_epu32 (_ct6_vec, const1792i);
                _ct12_vec = _mm256_add_epi32 (_ct12_vec, const1792i);
                _ct12_vec = _mm256_sub_epi32 (_ct12_vec, i0_vec);
                
                __m256 _arr_4_0_vec;
                float _arr_4_0_idx [8]; 
                for (int j = 0; j < 8; j++)
                {
                    int _ct12 = _mm256_extract_epi32 (_ct12_vec, j);
                    _arr_4_0_idx[j] = _arr_4_0[_ct12];
                }
                
                _arr_4_0_vec = _mm256_set_ps (_arr_4_0_idx[7], _arr_4_0_idx[6], _arr_4_0_idx[5], _arr_4_0_idx[4],
                                              _arr_4_0_idx[3], _arr_4_0_idx[2], _arr_4_0_idx[1], _arr_4_0_idx[0]);
                
                _mm256_storeu_ps (&_buf_2_2[((_i0 * (1 * 256)) + ((-256 * _T_i2) + _i2))],
                                  _mm256_add_ps (_buf_1_1_vec, _arr_4_0_vec));
              }
          }
          else
          {
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
                _buf_2_2[((_i0 * (1 * 256)) + ((-256 * _T_i2) + _i2))] = 
                (((beta * (_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] - (_i0 * 0.14285714285714285))) + 
                (_i0 * 0.14285714285714285)) + _arr_4_0[((_ct10 - (256 * _i0)) - -1792)]);
            }
          }
      }
      int _ct11 = ((1539 < ((256 * _T_i2) + 255))? 1539: ((256 * _T_i2) + 255));
      
      if (_ct11 == 256*_T_i2 + 255)
        {
            __m256 const1f = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            __m256i const6i = _mm256_set1_epi32 (6);
            __m256i const1x256 = _mm256_set1_epi32 (256);
            __m256i c50index = _mm256_set1_epi32 (-256 * _T_i2);
            __m256i _cv1i = _mm256_set1_epi32 (1);
            
            for (int _i2 = (256 * _T_i2); (_i2 < _ct11+1); _i2 = (_i2 + 8))
            {
                __m256 _buf_1_1v = _mm256_loadu_ps (&_buf_1_1[(0 + ((-256 * _T_i2) + _i2))]);
                __m256 _ct15vecf = _mm256_mul_ps (_buf_1_1v, const7f);
                
                __m256i _ct15 = _mm256_cvttps_epi32 (_ct15vecf);
                __m256i _ct17 = _mm256_max_epu32 (_ct15, _cv0i);
                __m256i _ct50 = _mm256_min_epu32 (_ct17, const6i);
                
                __m256 _ct15sub50 = _mm256_sub_ps (_ct15vecf, _mm256_cvtepi32_ps (_ct50));
                __m256 _1subct15sub50 = _mm256_sub_ps (const1f, _ct15sub50);
                
                __m256i _ct50_index = _mm256_mullo_epi32 (const1x256, _ct50);
                __m256i _i2vec = _mm256_set_epi32 (_i2 + 7, _i2 + 6, _i2 + 5, _i2 + 4, _i2 + 3, _i2 + 2, _i2 + 1, _i2);
                _i2vec = _mm256_add_epi32 (c50index, _i2vec);
                _ct50_index = _mm256_add_epi32 (_i2vec, _ct50_index);
                
                __m256i _ct50_1_index = _mm256_mullo_epi32 (const1x256, 
                                                            _mm256_add_epi32 (_ct50, _cv1i));
                _ct50_1_index = _mm256_add_epi32 (_i2vec, _ct50_1_index);
                
                __m256 _buf_2_2_index, _buf_2_2_index_1;
                float _buf_2_2_arr [8], _buf_2_2_arr_1[8];
                
                for (int j = 0; j < 7; j++)
                {
                    int _ct50 = _mm256_extract_epi32 (_ct50_index, j);
                    _buf_2_2_arr[j] = _buf_2_2[_ct50];
                }
                
                for (int j = 0; j < 7; j++)
                {
                    int _ct50 = _mm256_extract_epi32 (_ct50_1_index, j);
                    _buf_2_2_arr_1[j] = _buf_2_2[_ct50];
                }
                
                _buf_2_2_index = _mm256_set_ps (_buf_2_2_arr[7], _buf_2_2_arr[6], _buf_2_2_arr[5], _buf_2_2_arr[4],
                                                _buf_2_2_arr[3], _buf_2_2_arr[2], _buf_2_2_arr[1], _buf_2_2_arr[0]);
                
                _buf_2_2_index_1 = _mm256_set_ps (_buf_2_2_arr_1[7], _buf_2_2_arr_1[6], _buf_2_2_arr_1[5], _buf_2_2_arr_1[4],
                                                  _buf_2_2_arr_1[3], _buf_2_2_arr_1[2], _buf_2_2_arr_1[1], _buf_2_2_arr_1[0]);
                
                __m256 res = _mm256_add_ps(_mm256_mul_ps (_buf_2_2_index, _1subct15sub50), 
                                           _mm256_mul_ps (_buf_2_2_index_1, _ct15sub50));
                _mm256_storeu_ps (&_buf_1_3[(0 + ((-256 * _T_i2) + _i2))], res);
            }
        }
        else
        {
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
            _buf_1_3[(0 + ((-256 * _T_i2) + _i2))] = 
            ((_buf_2_2[((_ct20 * (1 * 256)) + ((-256 * _T_i2) + _i2))] * 
            (1.0 - ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)) - _ct29))) + 
            (_buf_2_2[(((_ct38 + 1) * (1 * 256)) + ((-256 * _T_i2) + _i2))] * ((_buf_1_1[(0 + ((-256 * _T_i2) + _i2))] * (float) (7)) - _ct47)));
            }
      }
      if (((_T_i1 <= 2559) && (_T_i2 <= 5)))
      {
        for (int _i0 = 0; (_i0 <= 2); _i0 = (_i0 + 1))
        {
            int __i2 = (256 * _T_i2);
            __m256 cv01f = {0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f};
            __m256 cv0f = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
            __m256 cv1f = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
            
            for (int _i2 = __i2; _i2 < __i2 + 256; _i2 = (_i2 + vSIZE))
            {
              __m256 _buf_1_3_vec = _mm256_loadu_ps(&_buf_1_3[(0 + ((-256 * _T_i2) + _i2))]);
              __m128i img_vec = _mm_loadu_si128 ((__m128i *)&img_colour[(((_i0 * ((4 + R) * (4 + C))) + (_T_i1 * (4 + C))) + _i2)]);
              __m128i xlo = _mm_unpacklo_epi16(img_vec, _mm_set1_epi16(0));
              __m128i xhi = _mm_unpackhi_epi16(img_vec, _mm_set1_epi16(0));
              __m128 ylo = _mm_cvtepi32_ps(xlo);
              __m128 yhi = _mm_cvtepi32_ps(xhi);
              __m256 img_vec_f = _mm256_set_m128 (yhi, ylo);
              __m256 _buf_1_1_vec = _mm256_loadu_ps (&_buf_1_1[(0 + ((-256 * _T_i2) + _i2))]);
              __m256 _ct51 =  _buf_1_3_vec * 
                       (( (img_vec_f) / cv65535) + cv01f) / 
                        (_buf_1_1_vec + cv01f);
              
              __m256 _ct53 = _mm256_max_ps (_ct51, cv0f);
              __m256 _ct59 = _mm256_min_ps (_ct53, cv1f);
              _ct59 = _ct59*cv65535;
              __m256i _res32i = _mm256_cvtps_epi32 (_ct59);
              __m128i _r = _mm256_castsi256_si128 (_res32i);
              _mm_storeu_si128 ((__m128i*)&laplacian[(((_i0 * (R * C)) + (_i1 * C)) + _i2)], _r);
            }
        }
      }
    }
  }
  pool_deallocate(_arr_4_0);
}
