# Euclidean PH1
This software will compute the degree-1 Vietoris-Rips Persistent Homology of a point cloud in low dimensional Euclidean Space. So far we have code for computing the degree-1 Vietoris-Rips persistent homology of point clouds in 2D and 3D. 

# How to use this software and warnings
This software has only been confirmed to be working on Ubuntu Linux, you may run into issues if you use other operating systems or distros. There are two files, one for 2D and one for 3D. One is called final_algorithm_2D.cpp, the other is called final_algorithm_3D.cpp

The algorithm this software is based on does rely on the assumption of pairwise unique distances. However, the software is quite resistant to rounding errors. 

# Limitations
How many points you are able to compute PH1 for is dependant on:
- The shape of your point cloud
- The dimension your point cloud resides in

# Requirements 
You will need the following files from nanoflann  https://github.com/jlblancoc/nanoflann
"nanoflann.hpp"
"KDTreeVectorOfVectorsAdaptor.h"

You will need to install CGAL, which can be found here https://www.cgal.org/download.html.
You will also need the Boost Library. https://www.boost.org. But since we are using CGAL this should be installed automatically along with CGAL.

# Your folder from which you will be executing the code from should look like the following. 

Final_algorithm_2D.cpp

README.md

final_algorithm_3D.cpp

nanoflann.hpp

KDTreeVectorOfVectorsAdaptor.h

point_cloud_to_analyse.csv

# Input file
This software takes input in the form of a csv file. Your csv file will either consist of 2 columns of numbers corresponding to 
the x and y co-ordinates of your points or will consist of 3 columns of numbers corresponding to the x,y and z co-ordinates of 
your points. 

# How to use the 2D file. 
You will need to make changes to lines 365 and 366 of Final_algorithm_2D.cpp. 
On line 365 you will need to change the number of neighbors.
On line 366 you will need to replace the file name with the name of the CSV file you want to use. 

To compile, you will need to type

g++ -std=c++14 -O3 Final_algorithm_2D.cpp -lmpfr -lgmp

and then run the corresponding executable. 

# How to use the 3D file. 
You will need to make changes to lines 365 and 366 of Final_algorithm_2D.cpp. 
On line 365 you will need to change the number of neighbors.
On line 366 you will need to replace the file name with the name of the CSV file you want to use. 

To compile, you will need to type

g++ -std=c++11 -O3 Final_algorithm_3D.cpp -lmpfr -lgmp

and then run the corresponding executable. 



