# Euclidean PH1
This software will compute the one-dimensional Vietoris-Rips Persistent Homology of a point cloud in low dimensional Euclidean Space. 

# How to use this software and warnings
This software has only been confirmed to be working on Ubuntu Linux, you may run into issues if you use other operating systems or distros. There are two files, one for 2D and one for 3D. One is called final_algorithm_2D.cpp, the other is called final_algorithm_3D.cpp

# Requirements for both 2D and 3D 
You will need the following files from nanoflann  https://github.com/jlblancoc/nanoflann
"nanoflann.hpp"
"KDTreeVectorOfVectorsAdaptor.h"

You will also need the Boost Library. https://www.boost.org. You will want to have the folder called "boost" in the same folder as this repository. However, if you want to use the 3D version you can ignore this step as CGAL will install this library for you automatically. 

# Requirement only for 2D
You will need Delaunator in order to compute the Delaunay triangulation for a 2D point cloud which can be found here https://github.com/abellgithub/delaunator-cpp. In particular, you will need the file "delaunator.hpp". 

# Requirement only for 3D 
You will need to install CGAL, which can be found here https://www.cgal.org/download.html. 

# If you want to use both the 2D and 3D version. Your folder from which you will be executing the code from should look like the following. 

Final_algorithm_2D.cpp

README.md

final_algorithm_3D.cpp

delaunator.hpp

nanoflann.hpp

KDTreeVectorOfVectorsAdaptor.h

boost (not necessary if you installed CGAL). 

# How to use the 2D file. 

