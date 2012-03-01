import sys
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
import math

# ============================================================================
# A simple 2D trimesh plotter and contourer utility that can be embedded in the C code
# syntax :
#        ariplot tri_mesh_file.dat data_file.dat plt_type xmin xmax ymin ymax title xlabel ylabel OUTPUT
# where :
#        tri_mesh_file.dat is a file containing connectivity integers (winding direction is arbitrary)
#             example:
#                     0 4 16
#                     2 1 6
#                     : : :
#                     6 18 0
#             where the first row is a triangle containing nodes 0,4 and 16. The second connects 2 1 6 
#             and so on.
#        data_file.dat is a file containing coordinates and solution for that point.
#             example:                     
#                     1.3 0.4 .7
#                     0.2 .4  0.19 
#                      :   :   :
#                     0.4 0.11 2.5
#             where for point (1.3,0.4) the value of solution is 0.7 and so on.
#         plt_type is the type of plot
#              for plt_type = Contour it plots contours efficiently
#              for plt_type = ColorTri it plots colored triangles with averages of solution at three vertices.
#         xmin, ymin, ymin, ymax, title, xlabel, ylabel are clear. No need for description.
#         OUTPUT shows the level of interation with user. For OUTPUT = display it shows the interactive toolbar
#         to the user. For OUTPUT = png,eps,... it writes the image to a file with name = title and quits.
#         
#
#
#  This is :   
   #             ######                      
  # #   #####  # #     # #       ####  ##### 
 #   #  #    # # #     # #      #    #   #   
#     # #    # # ######  #      #    #   #   
####### #####  # #       #      #    #   #   
#     # #   #  # #       #      #    #   #   
#     # #    # # #       ######  ####    #   

#                ,        ,         
#               /(        )`        
#               \ \___   / |        
#               /- _  `-/  '        
#              (/\/ \ \   /\        
#              / /   | `    \       
#              O O   ) /    |       
#              `-^--'`<     '       
#             (_.)  _  )   /        
#              `.___/`    /         
#                `-----' /          
#   <----.     __ / __   \          
#   <----|====O)))==) \) /====      
#   <----'    `--' `.__,' \         
#                |        |         
#                 \       /       /\
#            ______( (_  / \______/ 
#          ,'  ,-----'   |          
#          `--{__________)         
# Brought to you by Arash Ghasemi, ghasemi.arash@gmail.com
# License : BSD (see attached file)
#                2012
# ============================================================================
#                            main function
# ============================================================================ 

def main():
    gridf = sys.argv[1] #grid file
    dataf = sys.argv[2] #data file
    plt_type = sys.argv[3] #plot type
    xmin = float(sys.argv[4])
    xmax = float(sys.argv[5])
    ymin = float(sys.argv[6])
    ymax = float(sys.argv[7])
    title = sys.argv[8]
    xlabel = sys.argv[9]
    ylabel = sys.argv[10]
    OUTPUT = sys.argv[11]

#extracting data from data file ...
    infid = open(dataf,'r')
    lines = infid.readlines()
    x = []
    y = []
    z = []
    for i in lines:
        row = i.split()
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))
    infid.close()
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

#extracting triangle connectivity from grid file ...
    infid = open(gridf,'r')
    lines = infid.readlines()
    triangles = []
    for i in lines:
        row = i.split()
        triangles.append( [int(row[0]), int(row[1]), int(row[2])] )
    infid.close()
    triangles = np.asarray(triangles)
#preparing plot
    plt.figure()
    plt.gca().set_aspect('equal')
    if plt_type == 'Contour':
        plt.tricontourf(x, y, triangles, z)
    elif plt_type == 'ColorTri':
        plt.tripcolor(x, y, triangles, z, shading='faceted')
    else:
        print 'Unknown plot type!'
    plt.colorbar()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.ylim([ymin,ymax])
    plt.xlim([xmin,xmax])
    if OUTPUT == 'display':
        plt.show()
    else:
        plt.savefig(title +'.'+ OUTPUT,format = OUTPUT)

if __name__ == "__main__":
    main()
