# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:00:00 2016

@author: Luke
"""

import numpy as np
import matplotlib.pyplot as plt

# Circle
#func = 'x**2+y**2-10'
# Cardiod
# func = '((x-3)**2+y**2+5*(x-3))**2-5**2*((x-3)**2+y**2)'
#Cassini oval https://en.wikipedia.org/wiki/Implicit_curve
#func = '((x**2+y**2)**2-2*5**2*(x**2-y**2)-(5**4-5**4))'
# wavey sruface
#func = 'np.sin(x+y)-np.cos(x*y)+1'
# Example
func = '5*x**3 -17.3 * y**2 + np.sin(x*y)'

delta = 0.5

x = np.arange(-5, 5, delta)
y = np.arange(-5, 5, delta)

X, Y = np.meshgrid(x, y)

FN = lambda x,y: eval(func)

F = FN(X,Y)

cs = plt.contour(X, Y, F, [0])
paths = cs.collections[0].get_paths()

plt.show()

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