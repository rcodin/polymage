## Requirements

* python-2.7
* ctypes
    * ** $ sudo apt-get install python-ctypeslib**
* numpy
    * ** $ sudo apt-get install python-numpy**
* numba
    * Get 32 or 64 bit version of [Anaconda](https://continuum.io/downloads.html)
        * **$ bash Anaconda2-4.0.0-Linux-x86_64.sh**  [... _follow instructions_]
        * **$ conda install numba**
* OpenCV
    * Install Library
        * [[instructions for installation on ubuntu]](https://help.ubuntu.com/community/OpenCV)
    * Install cv2
        * **$ sudo apt-get install python-opencv**
* PIL (Python Imaging Library)
    * **$ pip install pillow**
* tabulate
    * **$ pip install tabulate**

## How To Use The Video Demo

##### To compile all the demo apps, run :
> **$ make**

##### To play the video demo, run the python script as :
> **$ python video_demo.py /path/to/video/file**

The implementations work for a generic resolution.  
For a sample video, try:
[Big Buck Bunny](https://peach.blender.org/download/)
or
[DivX Video Samples](http://www.divx.com/en/devices/profiles/video)


## How To Switch Between Video Modes

#### Toggle between apps:
* **'h'**   : **H**arris Corner Detection
* **'u'**   : **U**nsharp Mask
* **'b'**   : **B**ilateral Grid
* **'l'**   : **L**ocal Laplacian Filters
    * _lowercase **L**_
* **'ESC'** : Original

#### Toggle between app modes:
* **' '** : OpenCV [1]
    * (_space_)
* **'n'** : PolyMage **N**aive
* **'o'** : PolyMage **O**pt
* **'p'** : **P**IL [2]
* **'j'** : Numba [2]
    *  Numpy **J**it
* **'q'** - **Q**UIT

--
#####Currently available for :
[1] Unsharp Mask, Harris Corner  
[2] Unsharp Mask


#####To be able to toggle between modes in an app, at least two modes of the app must be visited.