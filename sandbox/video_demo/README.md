### Requirements

* python-2.7
* ctypes
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> ** $ sudo apt-get install python-ctypeslib**</span>
* numpy
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> ** $ sudo apt-get install python-numpy**</span>
* numba
 * Get 32 or 64 bit version of [Anaconda](https://continuum.io/downloads.html)
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ bash Anaconda2-4.0.0-Linux-x86_64.sh**</span>  [... _follow instructions_]
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ conda install numba**</span>
* OpenCV
 * Install Library [<sup>instructions for ubuntu users</sup>](https://help.ubuntu.com/community/OpenCV)
 * Install cv2 <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ sudo apt-get install python-opencv**</span>
* PIL (Python Imaging Library)
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ pip install pillow**</span>
* tabulate
 * <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ pip install tabulate**</span>

### How To Use The Video Demo

##### To compile all the demo apps, run :
> <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ make** </span>

##### To play the video demo, run the python script as :
> <span style="color:lime; background-color:black; font-family:Courier; font-size:1em"> **$ python video_demo.py /path/to/video/file** </span>

The implementations work for a generic resolution.  
For a sample video, try:
[Big Buck Bunny](https://peach.blender.org/download/)
or
[DivX Video Samples](http://www.divx.com/en/devices/profiles/video)


### How To Switch Between Video Modes

#### Toggle between apps:
* **'h'**   : **H**arris Corner Detection
* **'u'**   : **U**nsharp Mask
* **'b'**   : **B**ilateral Grid
* **'l'**   : **L**ocal Laplacian Filters
 * _lowercase **L**_
* **'ESC'** : Original

#### Toggle between app modes:
* **' '** : OpenCV <sup>1</sup>
 * (_space_)
* **'n'** : PolyMage **N**aive
* **'o'** : PolyMage **O**pt
* **'p'** : **P**IL <sup>2</sup>
* **'j'** : Numba <sup>2</sup>
 *  Numpy **J**it
* **'q'** - **Q**UIT

--
#####Currently available for :
<sup>1</sup> Unsharp Mask, Harris Corner  
<sup>2</sup> Unsharp Mask


#####To be able to toggle between modes in an app, at least two modes of the app must be visited.