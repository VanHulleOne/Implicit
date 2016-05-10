# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 16:00:00 2016

@author: Luke
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import geometry as g
from Geom import Polygon, Point, Line
from collections import deque

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
#func = '((x**2+y**2)**2-2*5**2*(x**2-y**2)-(5**4-5**4))'
# wavey surface
#func = 'np.sin(x+y)-np.cos(x*y)+1'
# Example
#func = '5*x**3 -17.3 * y**2 + np.sin(x*y)'
# Distance
func = '(-4+(x**2+y**2)**0.5)'

delta = 0.5
maxError = np.sqrt(2*delta**2)
size = 10

x = np.arange(-size, size, delta)
y = np.arange(-size, size, delta)

xGrid, yGrid = np.meshgrid(x, y)

FN = lambda x,y: eval(func)

F = FN(xGrid, yGrid)

cs = plt.contour(xGrid, yGrid, F, [0], colors = 'r')
paths = cs.collections[0].get_paths()

plt.show()
print '\n' + func

print '\nNumber of connected components is: %d' %len(paths)

prev = None
pointList = []
for coord in paths[len(paths)/2].vertices:
    point = Point(coord)
    if point != prev:
        pointList.append(point)
        prev = point
    

def pairwise(l1):
    l1Iter = iter(l1)
    prev = next(l1Iter)
    for curr in l1Iter:
        yield (prev, curr)
        prev = curr
    
poly = Polygon(pairwise(pointList))

print poly.simplyConnected()
print ''

shape = g.Polygon(*paths[0].vertices)  
csHigher = plt.contour(xGrid, yGrid, F,[delta], colors = 'g')
pathsHigher = csHigher.collections[0].get_paths()
shapeHigher = g.Polygon(*pathsHigher[0].vertices)

isEmpty = abs(shapeHigher.area) > abs(shape.area)

if isEmpty:
    print 'The figure has an empty interior'
else:
    print 'The figure has a non-empty interior'

print ''
print 'The max distance from the approximation is {:.3f}'.format(maxError)
print ''

def quick2ndDeriv_coro(num, tolerance):
    first = rollingAvg_coro(num)
    next(first)
    isSharp = False
    for i in xrange(num):
        temp = yield isSharp
        first.send(temp)
    second = rollingAvg_coro(num)
    next(second)
    q = deque()
    for i in xrange(num):
        temp = yield isSharp
        second.send(temp)
        q.append(temp)
    while 1:
        temp = yield isSharp
        q.append(temp)
        firstAvg = first.send(q.popleft())
        secondAvg = second.send(temp)
        if(abs(firstAvg-secondAvg) > tolerance):
            isSharp = True
            
def rollingAvg_coro(num):
    q = deque()
    total = 0
    average = None
    for i in xrange(num):
        temp = yield average 
        total += temp
        q.append(temp)
        average = 1.0*total/len(q)        
        
    while 1:
        temp = yield 1.0*total/num
        total += temp - q.popleft()
        q.append(temp)

def sharpCorner(lines, num, tolerance):
    sharpX = quick2ndDeriv_coro(num, tolerance)
    next(sharpX)
    sharpY = quick2ndDeriv_coro(num, tolerance)
    next(sharpY)
    for line in lines + lines[:6]:
        normalizedSlope = line.normalizedSlope()
        if(sharpX.send(normalizedSlope[X]) or sharpY.send(normalizedSlope[Y])):
            return True
    return False
    
def isDistance(line1, line2):
    midPoint = line1.end
    ''' To find the bisector I take the normalized slopes (length of line = 1)
    and add that to middle point.'''
    p1Hat = Point(midPoint.pointVector-line1.normalizedSlope())
    p3Hat = Point(midPoint.pointVector+line2.normalizedSlope())
    ''' half way between pHats is the bisector '''
    testPoint = Point((p1Hat.pointVector+p3Hat.pointVector)/2.0)
    ''' put that point into our function '''
    result = FN(testPoint.x, testPoint.y)
    if isEmpty and result < 0:
        ''' Here i test to see if the point is inside or outside and then
        switch it to outside if necessary to make sure i'm not getting closer
        to a different part of the line.'''
        bisectLine = Line([testPoint, midPoint])
        testPoint = Point(midPoint.pointVector + bisectLine.normalizedSlope())
    ''' if the abs distance is less than my grid's max error return True '''
    if abs(midPoint-testPoint - FN(testPoint.x, testPoint.y)) < maxError:
        return True
    return False
    
lineArray = np.array(poly.lines[:len(poly.lines) if len(poly.lines)%2 == 0
                        else len(poly.lines)-1])
lineArray = lineArray.reshape((len(lineArray)/2,2))
      
if sharpCorner(poly.lines, 2, 0.4):
    print 'The shape has at least one sharp corner.'
else:
    print 'The shape does not have sharp corners.'
    
print ''
if all(isDistance(*linePair) for linePair in lineArray):
    print 'The function is the distance function.'
else:
    print 'The function is not the distance function.'
      
        
        
        
        
        
        
        
        
        
        