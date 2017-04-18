## Requirements

* python-2.7
* ctypes
    * ** $ sudo apt-get install python-ctypeslib**
* numpy
    * ** $ sudo apt-get install python-numpy**
* numba
    Option 1
    * Get 32 or 64 bit version of [Anaconda](https://continuum.io/downloads.html)
        * **$ bash Anaconda2-4.0.0-Linux-x86_64.sh**  [... _follow instructions_]
        * **$ conda install numba**

    Option 2
       Alternatively, one can install numba through pip.


* OpenCV
    * Install Library
        * [[instructions for installation on ubuntu]](https://help.ubuntu.com/community/OpenCV)
    * Install cv2
        * On an Ubuntu: **$ sudo apt-get install python-opencv**
        * On a Fedora: **$ sudo yum install opencv-python
* PIL (Python Imaging Library)
    * **$ pip install pillow**
* tabulate
    * **$ pip install tabulate**

## How To Use The Video Demo

##### To compile all the demo apps, run :
> **$ make**

##### To play the video demo, run the python script as :
> **$ python video_demo.py /path/to/video/file**

**Sample (1) from Demo**  
![Sample 1][pic1]  

**Sample (2) from Demo**  
![Sample 2][pic2]

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

[pic1]: https://bytebucket.org/udayb/polymage/raw/3cc5b191d0e4b8c4f85bbae2253d56b697c723f8/sandbox/video_demo/screenshots/video_demo_screenshot1.png
[pic2]: https://bytebucket.org/udayb/polymage/raw/3cc5b191d0e4b8c4f85bbae2253d56b697c723f8/sandbox/video_demo/screenshots/video_demo_screenshot2.png
