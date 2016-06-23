Compiling the C++ file into a shared library
============================================

To compile all demo apps, run

$ make

To play the video demo, run the python script as

$ python video_demo.py path/to/video/file

The implementations work for a generic resolution.

For a sample video, try:
https://peach.blender.org/download/
or
http://www.divx.com/en/devices/profiles/video


Switch Options (key strokes)

Toggle between apps:
===================
'h'   - Harris Corner Detection
'u'   - Unsharp Mask
'b'   - Bilateral Grid
'l'   - Local Laplacian Filters
'ESC' - Original

Toggle between modes:
====================
' ' - OpenCV (*1)
'n' - PolyMage Naive
'o' - PolyMage Opt
'p' - PIL (*2)
'j' - Numba [Numpy jit] (*2)

'q' - QUIT


--------------------------------
* currently available for :
--------------------------------
*1 : Unsharp Mask, Harris Corner
*2 : Unsharp Mask
--------------------------------


-- To be able to toggle between two modes in an app, visit at least
   two modes of the app.
