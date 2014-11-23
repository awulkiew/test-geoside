test-geoside
============

This program draws a grid of geographical points [-50,-50]x[50,50]. Colors depends on the side of a point WRT geographical segment (-51,-51)->(51,51).

The colors are:
* dark green - left
* dark red - right
* cyan - left on spheroid, right on sphere
* magenta - right on spheroid, left on sphere

Sphere or spheroid with flattening = 0

![f=0](f0.png)

Spheroid with flattening = 1/300 (Earth)

![f=1/300(earth)](fearth.png)

Spheroid with flattening = 0.25

![f=0.25](f0.25.png)

Spheroid with flattening = 0.5

![f=0.5](f0.5.png)
