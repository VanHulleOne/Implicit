# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:00:00 2016

@author: Luke
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import geometry as g
from collections import Counter
from Geom import Polygon, Point, Line

# The X and Y offsets in a point
X, Y = 0, 1
# The return results for a left right test
LEFT, ON_EDGE, RIGHT = -1,0,1
# For the line sweep event queue these are flags in the Point class
# that tell if the point is the start, end, or intersection of a line
START, END, INTER = 0,1,2
# The distance tolerance given for our points
EPSILON = 10e-6
# The arctan of a small number is just itself so use the same value
# for angular tests
ANGLE_EPS = EPSILON

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

prev = None
pointList = []
for coord in paths[0].vertices:
    point = Point(coord)
    if point != prev:
        pointList.append(point)
        prev = point
    

def pairwise(l1):
    l1Iter = iter(l1)
    first = prev = next(l1Iter)
    for curr in l1Iter:
        yield (prev, curr)
        prev = curr
    yield (prev, first)
    
poly = Polygon(pairwise(pointList))

#shape = g.Polygon(*paths[0].vertices)  





    
