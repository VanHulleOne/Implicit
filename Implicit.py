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
func = '((x**2+y**2)**2-2*5**2*(x**2-y**2)-(5**4-5**4))'
# wavey surface
#func = 'np.sin(x+y)-np.cos(x*y)+1'
# Example
#func = '5*x**3 -17.3 * y**2 + np.sin(x*y)'
# Distance
#func = '(-4+(x**2+y**2)**0.5)'
h = 0.5
delta = np.sqrt(h**2/2)
maxError = h#np.sqrt(2*delta**2)
size = 10

x = np.arange(-size, size, delta)
y = np.arange(-size, size, delta)

xGrid, yGrid = np.meshgrid(x, y)

FN = lambda x,y: eval(func)

F = FN(xGrid, yGrid)

cs = plt.contour(xGrid, yGrid, F, [0], colors = 'r')
csf = plt.contourf(xGrid, yGrid, F, [0, np.inf], colors = 'c')
paths = cs.collections[0].get_paths()

plt.show()
print '\n' + func

print '\nNumber of connected components is: %d' %len(paths)

prev = None
pointList = [] # the list wich will store all of the points in the shape
for coord in paths[len(paths)/2].vertices:
    """
    Due to the grid layout some points are repeated so we remove the
    duplicates when creating the point list. For functions which produce multiple
    contours we are only testing the middle contour.
    """
    point = Point(coord)
    if point != prev:
        pointList.append(point)
        prev = point

def pairwise(l1):
    """
    A generator to yield points in line pairs.
    
    The points are the boundry of the shape, to turn those into lines
    we need to get them into pairs where the shared point between two
    lines is repeated. This method does that.
    
    Parameter
    ---------
    l1 - the list to pair
    
    Yield
    -----
    A tuple of the elements paired
    """
    l1Iter = iter(l1)
    prev = next(l1Iter)
    for curr in l1Iter:
        yield (prev, curr)
        prev = curr

""" Create the shape from my Geom module which is a pared down version of
my submissions for Program 2 """    
poly = Polygon(pairwise(pointList))

""" The shape is simply connected if all of the points appear exactly twice."""
print poly.simplyConnected()
print ''

""" Use SymPy geometry module to create the shape so we can get its area. """
shape = g.Polygon(*paths[0].vertices)  

"""
The contour is drawn with the positive Z direction being to the left
so if shape has a positive area then the interior is non-empty. If shape
has a negative area then the interior is empty.
"""
if shape.area < 0:
    print 'The figure has an empty interior'
else:
    print 'The figure has a non-empty interior'

print ''
print 'The max distance from the approximation is {:.3f}'.format(maxError)
print ''

def quick2ndDeriv_coro(num, tolerance):
    """
    yields True if the 2nd derivative is larger than tolerance.
    
    This coroutine takes yields in values and then stores two sets of
    rolling averages. The two averages are each num long with the falling
    off value of second being sent into first. If the difference between the
    two averages become larger than tolerance then True is yielded else False
    is yieled.
    
    Parameters
    ----------
    num - an int for the length of the rolling average
    tolerance - above this and True is yielded
    
    Yields
    ------
    in - the values to be tracked
    out - True if the two averages have a difference greater than tolerance
    else False.
    """
    
    """ Setting up the first rolling average. """
    first = rollingAvg_coro(num)
    next(first)
    isSharp = False
    for i in xrange(num):
        temp = yield isSharp
        first.send(temp)
    """ Setting up the second rolling average."""
    second = rollingAvg_coro(num)
    next(second)
    q = deque()
    for i in xrange(num):
        temp = yield isSharp
        second.send(temp)
        q.append(temp)
    while 1:
        """ Continuing the rolling averages after both have been initialized. """
        temp = yield isSharp
        q.append(temp)
        firstAvg = first.send(q.popleft())
        secondAvg = second.send(temp)
        if(abs(firstAvg-secondAvg) > tolerance):
            isSharp = True
            
def rollingAvg_coro(num):
    """
    A coroutine which yields a rolling average num long.
    
    The coroutine yields in values and yields out the current rolling average.
    The roll is num long.
    
    Parameter
    ---------
    num - an int for the length of the rolling average
    
    Yields
    ------
    in - the values to be averaged
    out - the current rolling average
    """
    
    q = deque()
    total = 0
    average = None
    for i in xrange(num):
        """ This sets up the initial num values. """
        temp = yield average 
        total += temp
        q.append(temp)
        average = 1.0*total/len(q)        
        
    while 1:
        temp = yield 1.0*total/num
        total += temp - q.popleft()
        q.append(temp)

def sharpCorner(lines, num, tolerance):
    """
    Returns True if lines contains a sharp corner else False.
    
    Since angles wrap around that did not seem like a good way to test
    if there were sharp angles. Instead decided to normalize the length of each
    segment to one and then find its deltaX and deltaY, what I have called
    the normalized slope. Then all we have to do is compare the normalized
    slope of one line to the next and if either delta is > tol there is a sharp
    corner. Due to some irregularities in the contour I found that taking the
    rolling average of a few lines in a row helps smooth out the data and
    provide better results. From testing a rolling overage of 2 or 3 seems
    sufficient. A tolerance of 0.4 is around 23 degrees and also seems to work.
    
    Parameters
    ----------
    lines - The list of Lines to test
    num - The length of the rolling average
    tolerance - The max different in delta averages allowed (normalized to one)
    
    Return
    ------
    True if tolerance is exceeded else False    
    """
    
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
    """
    Tests if the input function is the distance function.
    
    An angle bisector is used to pick a test point outside of the shape
    nearest an exhisting point. If the distance from the shape point to the
    test point is the same value (within a tolerance) as func(testPoint) then
    the function is the distance function.
    
    Parameters
    ----------
    line1, line2 - The two adjacent lines to use for the test
    
    Return
    ------
    True if the distance and Func return the same value (within a tolerance)
    """
    
    """ The shared point of the two lines. """
    midPoint = line1.end 
    """
    To find the bisector I take the normalized slopes (length of line = 1)
    and add that to middle point.
    """
    p1Hat = Point(midPoint.pointVector-line1.normalizedSlope())
    p3Hat = Point(midPoint.pointVector+line2.normalizedSlope())
    """ half way between pHats is the bisector the other point of the bisector """
    testPoint = Point((p1Hat.pointVector+p3Hat.pointVector)/2.0)

    result = FN(testPoint.x, testPoint.y)
    if not (isEmpty ^ result < 0):
        """
        If the shape is empty and the result is negative or if the shape is not
        empty and the result is positive then our test point is inside the shape
        so pick a new one outside the shape.
        """
        bisectLine = Line([testPoint, midPoint])
        testPoint = Point(midPoint.pointVector + bisectLine.normalizedSlope())
    ''' if the abs distance is less than my grid's max error return True '''
    if abs(midPoint-testPoint - FN(testPoint.x, testPoint.y)) < maxError:
        return True
    return False


      
if sharpCorner(poly.lines, 3, 0.4):
    print 'The shape has at least one sharp corner.'
else:
    print 'The shape does not have sharp corners.'

""" Create a NumPy array of an even length. """   
lineArray = np.array(poly.lines[:len(poly.lines) if len(poly.lines)%2 == 0
                        else len(poly.lines)-1])
""" Reshape the array to put the lines in pairs. """
lineArray = lineArray.reshape((len(lineArray)/2,2))
   
print ''
"""
If we only test one point if it matches the distance function we may have
just found an anomoly so pair the lines and test the pairs. This creates
a test point outside every other point in the shape. If any do not match
the distance then we mark it as not a distance function.
"""
if all(isDistance(*linePair) for linePair in lineArray):
    print 'The function is the distance function.'
else:
    print 'The function is not the distance function.'
      
        
        
        
        
        
        
        
        
        
        