test-geoside
============

This repository contains 2 test programs for the purpose of testing spherical and geographical algorithms and strategies.

#### spheroid.cpp

This program is a 3D visualizer of a spheroid, two points with their longitudes and latitudes and curves on the surface representing a segment between those points. Curves are calculated using various methods.

Axes:
* X - red
* Y - green
* Z - blue

Points:
* p1 - orange
* p2 - yellow

Curves:
* experimental - orange -> yellow
* geodetic mapping - red
* geocentric mapping - green
* reduced mapping - blue
* great ellipse - magenta
* vincenty - gray->white

(NOTE: great ellipse and mappings overlaps)

UI
* Mouse - navigation
* WSAD  - move p1
* IKJL  - move p2
* h     - display help
* ,     - decrease b, increase flattening
* .     - increase b, decrease flattening
* m     - experimental method switch
* 1     - experimental curve on/off
* 2     - geodetic curve on/off
* 3     - geocentric curve on/off
* 4     - reduced curve on/off
* 5     - great ellipse curve on/off
* 6     - vincenty curve on/off

Spheroid with flattening = 0.25

![sph](sph.png)

#### geoside.cpp

This program draws a grid of geographical points. Colors depends on the side of a point WRT geographical segment. The side is calculated using 2 methods:

1. spherical side formula transforming spherical coordinates into 3D cartesian vectors and checking the side using vector algebra.
2. comparison of azimuths calculated using Vincenty's inverse formula, azimuth1 of segment points SP1 and SP2, azimuth2 of SP1 and P.

The colors below indicates:
* gray - left
* yellow - right
* red - right on spheroid, left on sphere
* green - left on spheroid, right on sphere

Spheroid with flattening = 0.25

![f=0.25](f0.25.png)
