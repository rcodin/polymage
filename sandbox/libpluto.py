from ctypes import cdll, Structure, c_int, c_double, c_uint


pluto = cdll.LoadLibrary('/usr/local/lib/libpluto.so')
print('Loaded lib {0}'.format(pluto))
