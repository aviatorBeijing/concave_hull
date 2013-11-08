concave_hull
============

C/Python code to generate 2D alpha-shape (concave hull). I am only concerned about the 2D points set. No higher dimension is considered. 
In the following references, all C routines and Python script are provided by other authors. But it just doesn't work right away. 
I need to modify the code here and there. So to save the effort in the future, I modified them, and uploaded here.

Reference:
http://www.netlib.org/voronoi/hull.html
http://stackoverflow.com/questions/6833243/how-can-i-find-the-alpha-shape-concave-hull-of-a-2d-point-cloud

Build:
On Linux platform:
1) unzip concave_hull-master.zip
2) cd concave_hull-master
3) make
4) make
5) python hull.py
6) Insert the output in the "google.html"
7) View "google.html" in your favorite HTML viewer. Vola!
