# UnstEuler2D
Completely Validated Unstructured Two-Dimensional Euler Solver

# Quick Start

1- Make sure you have python and matplotlib installed on your computer.

2- You need to have a C++ compiler. I recommend Gcc compilers but you can use any compiler by selecting your compiler in the makefile.

3- At terminal type make. This bulids the binary file "unst2d"

4- Run the code for the following arguments

$> ./unst2d ./grids/naca0012.mesh .8 1.25

The code runs the NACA0012 airfoil at Machnumber 0.8 with angle of attack 1.25 degrees. It converges and finally plots the following solution on screen. Then you are ready to go ahead and change the code and use it for your own work!

![alt tag](https://raw.github.com/arrgasm/UnstEuler2D/master/grids/density.png)
