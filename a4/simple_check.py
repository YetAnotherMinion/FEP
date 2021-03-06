# Simple 2 element quad mesh with quadratic shape functions
# fixed in the x about the left edge, and fixed in the y
# on the bottom edge, with a traction applied to the right
# edge in the positive x direction
#
# Each element is a unit square in the global coordinates
# 
#
# (22,23)-----(14,15)-----(12,13)-----( 2, 3)-----( 0, 1)
#     |                       |                       |
#     |                       |                       |
#     |                       |                       |
# (24,25)     (16,17)     (18,19)     ( 4, 5)     ( 6, 7)
#     |                       |                       |
#     |                       |                       |
#     |                       |                       |    
# (28,29)-----(26,27)-----(20,21)-----(10,11)-----( 8, 9)
#
#
# Element 1       Element 2   
# 1 { 0  -> 28    1 { 0  -> 20
#     1  -> 29        1  -> 21
# 2 { 2  -> 20    2 { 2  ->  8
#     3  -> 21        3  ->  9
# 3 { 4  -> 12    3 { 4  ->  0
#     5  -> 13        5  ->  1
# 4 { 6  -> 22    4 { 6  -> 12
#     7  -> 23        7  -> 13
# 5 { 8  -> 26    5 { 8  -> 10
#     9  -> 27        9  -> 11
# 6 { 10 -> 18    6 { 10 ->  6
#     11 -> 19        11 ->  7
# 7 { 12 -> 14    7 { 12 ->  2
#     13 -> 15        13 ->  3
# 8 { 14 -> 24    8 { 14 -> 18
#     15 -> 25        15 -> 19
# 9 { 16 -> 16    9 { 16 ->  4
#     17 -> 17        17 ->  5

import scipy
import numpy
import matplotlib.pyplot

IEN=(
        (   (28, 29),
            (20, 21),
            (12, 13),
            (22, 23),
            (26, 27),
            (18, 19),
            (14, 15),
            (24, 25),
            (16, 17)
        ), (
            (20, 21),
            (8, 9),
            (0, 1),
            (12, 13),
            (10, 11),
            (6, 7),
            (2, 3),
            (18, 19),
            (4, 5)
        )
    )
#the location of each element is given nodewise
x_offset = 3.0
y_offset = 4.0
cell_size = 1.0

coordinates = (x_offset+cell_size*2,    y_offset+cell_size,
                x_offset+1.5*cell_size, y_offset+cell_size,
                x_offset+1.5*cell_size, y_offset+cell_size/2,
                x_offset+2*cell_size,   y_offset+cell_size/2,
                x_offset+2*cell_size,   y_offset,
                x_offset+1.5*cell_size, y_offset,
                x_offset+cell_size,     y_offset+cell_size,
                x_offset+cell_size/2,   y_offset+cell_size,
                x_offset+cell_size/2,   y_offset+cell_size/2,
                x_offset+cell_size,     y_offset+cell_size/2,
                x_offset+cell_size,     y_offset,
                x_offset,               y_offset+cell_size,
                x_offset,               y_offset+cell_size/2,
                x_offset+cell_size/2,   y_offset,
                x_offset,               y_offset
    )
#we know the global number of dofs for this toy problem
nGlobalDOFs = 30
K = numpy.zeros((nGlobalDOFs,nGlobalDOFs))
F = numpy.zeros((nGlobalDOFs,1))
#plot the mesh coordinates for visual check of numbering order
elm_x = numpy.zeros(nGlobalDOFs/2)
elm_y = numpy.zeros(nGlobalDOFs/2)
#prepare the shape functions and prepare the plot at the same time
shapeVals = (
    lambda e,n : (e*n)*(e-1)*(n-1)/4.0,
    lambda e,n : (e*n)*(e+1)*(n-1)/4.0,
    lambda e,n : (e*n)*(e+1)*(n+1)/4.0,
    lambda e,n : (e*n)*(e-1)*(n+1)/4.0,
    lambda e,n : (1-(e*e))*(n-1)*n/2.0,
    lambda e,n : (e+1)*(1-(n*n))*e/2.0,
    lambda e,n : (1-(e*e))*(n+1)*n/2.0,
    lambda e,n : (e-1)*(1-(n*n))*e/2.0,
    lambda e,n : (1-(e*e))*(1-(n*n))
    )

shapeGrads = (
    lambda e,n : ( (e - 0.5)*(n*(n-1))/2.0, (n - 0.5)*(e*(e-1))/2.0, 0.0),
    lambda e,n : ( (e + 0.5)*(n*(n-1))/2.0, (n - 0.5)*(e*(e+1))/2.0, 0.0),
    lambda e,n : ( (e + 0.5)*(n*(n+1))/2.0, (n + 0.5)*(e*(e+1))/2.0, 0.0),
    lambda e,n : ( (e - 0.5)*(n*(n+1))/2.0, (n + 0.5)*(e*(e-1))/2.0, 0.0),
    lambda e,n : ( -e*(n*(n-1)), (n - 0.5)*(1-(e*e)), 0.0),
    lambda e,n : ( (e + 0.5)*(1-(n*n)), -n*(e*(e+1)), 0.0),
    lambda e,n : ( -e*(n*(n+1)), (n + 0.5)*(1-(e*e)), 0.0),
    lambda e,n : ( (e - 0.5)*(1-(n*n)), -n*(e*(e-1)), 0.0),
    lambda e,n : ( -2.0*e*(1-(n*n)), -2.0*n*(1-(e*e)), 0.0)
    )

K = numpy.zeros(shape = (nGlobalDOFs,nGlobalDOFs))

for elm_indx, elm in enumerate(IEN):
    elm_sz = len(elm)
    ke = numpy.zeros(shape =(elm_sz, elm_sz))
    for node_indx, node in enumerate(elm):
        for dim_indx, dim in enumerate(node):
            #convention assumes that first dimension is x
            if 0 == dim_indx:
                elm_x[dim/2] = coordinates[dim]
            elif 1 == dim_indx:
                elm_y[(dim-1)/2] = coordinates[dim]
            else:
                #should never reach her
                raise(RuntimeError)

# UNCOMMENT below to visually check mesh node numberings
# matplotlib.pyplot.scatter(elm_x, elm_y)
# labels = ['({0},{1})'.format(ii,ii+1) for ii in range(0, nGlobalDOFs,2)]
# for label, x, y in zip(labels, elm_x, elm_y):
#     matplotlib.pyplot.annotate(label, xy = (x,y), xytext = (15,15),
#         textcoords = 'offset points', ha = 'left', va = 'bottom')
# matplotlib.pyplot.show()

