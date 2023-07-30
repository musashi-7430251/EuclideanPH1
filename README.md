# Euclidean PH1
This software will compute the one-dimensional Vietoris-Rips Persistent Homology of a point cloud in low dimensional Euclidean Space. 

# How to use this software
This software has only been confirmed to be working on Ubuntu Linux, you may run into issues if you use other operating systems or distros. There are two files, one for 2D and one for 3D. 

# Requirements for both 2D and 3D 
You will need the following files from nanoflann  https://github.com/jlblancoc/nanoflann
"nanoflann.hpp"
"KDTreeVectorOfVectorsAdaptor.h"

You will also need the Boost Library. https://www.boost.org. You will want to have the folder called "boost" in the same folder as this repository. However, if you want to use the 3D version you can ignore this step as CGAL will install this library for you automatically. 

# Requirement only for 2D
You will need 
