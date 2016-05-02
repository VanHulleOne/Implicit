# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:00:00 2016

@author: Luke
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import geometry as g

# Circle
#func = 'x**2+y**2-10'
# Cardiod
#func = '((x-3)**2+y**2+5*(x-3))**2-5**2*((x-3)**2+y**2)'
#Cassini oval https://en.wikipedia.org/wiki/Implicit_curve
func = '((x**2+y**2)**2-2*5**2*(x**2-y**2)-(5**4-5**4))'
# wavey surface
#func = 'np.sin(x+y)-np.cos(x*y)+1'
# Example
#func = '5*x**3 -17.3 * y**2 + np.sin(x*y)'

delta = 0.5
size = 10

x = np.arange(-size, size, delta)
y = np.arange(-size, size, delta)

X, Y = np.meshgrid(x, y)

FN = lambda x,y: eval(func)

F = FN(X,Y)

cs = plt.contour(X, Y, F, [0])
paths = cs.collections[0].get_paths()

plt.show()

print 'Number of connected components is: %d' %len(paths)

def pairwise(l1):
    l1Iter = iter(l1)
    first = curr = next(l1Iter)
    for p in l1Iter:
        yield [curr, p]
        curr = p
    yield[curr, first]

shape = g.Polygon(*paths[0].vertices)  




#GRID_LENGTH = 100
#import matplotlib.pyplot
#from numpy import arange
#from numpy import meshgrid
#
#def detectShape(fn, h):
##Build Grid
#    x_range = arange(-50.0, 50.0, h)
#    y_range = arange(-50.0, 50.0, h)
#    X, Y = meshgrid(x_range,y_range)
#    
#    FN = lambda x,y: eval(fn)
#    
#    F = FN(X,Y)
#    
#    cs = matplotlib.pyplot.contour(X, Y, F, [0])
#    path = cs.collections[0].get_paths()
    

#Run tests

#Print # of connected components
#Determine if S is simply connected
#Determine if S has non-empty interior


#Max distance from approximation to actual boundary
#print "Max Distance from approximation to actual boundary:"
#print h*2**.5 #The furthest distance our test could be off is where
 #the intersection touches the corner of a cell, but the actual
 #curve extends to the far corner
#print


#Draw shape
#    matplotlib.pyplot.show()