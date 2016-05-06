# -*- coding: utf-8 -*-
"""
Created on Fri May 06 14:43:46 2016

@author: lvanhulle
"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import geometry as g
from collections import Counter

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

class Polygon:
    
    def __init__(self, pointGen):
        """
        Does the hard work and calls all of the required functions.
        """
        self.lines = set() # Store the lines in a set for unorered fast retrieval
        self.pointCount = Counter() # Counts the number of times a point appears
        self.points = [] # a List of the points
        self.addPoints(pointGen)
        self.runRest()
        
    def addPoints(self, pointGen):
        for p1, p2 in pointGen:
            """
            Each line from the file is split into two lists, one for
            each point it represents. This points are sorted so that
            the start of the line is the minimum X value.
            """
            line = Line(sorted(p1.pointVector, p2.pointVector))

            self.lines.add(line)
            for point in line:
                """
                For each point in the line due two steps. The first is
                "Add" it to the pointCount. Each coordinates for a point
                are rounded to EPSILON in the point class. From there we
                assume that each point has the exact same value. That
                XY value is turned into a string and used as a key in
                pointCount. If the key does not exist (the point has not yet
                been added) it is added and then its count is set
                to one. If the point is already there the count is incremented.
                
                The second setep is add it to the points list.
                """
                self.pointCount[str(point)] += 1
                self.points.append(point)
                
    def runRest(self):      
        if any(value != 2 for value in self.pointCount.itervalues()):
            """
            If any value in the pointCount does not equal two we know
            we do not have a valid polygon. Either three or more points
            are coincident or the shape is not closed.
            """
            
            """
            Create a list of tuples containing the point and count of
            any point whose count does not equal 2.
            """
            violators = [(point, self.pointCount[point])
                        for point in self.pointCount
                        if self.pointCount[point] != 2]
            
            # Create the message to be displayed            
            message = ['\nThe following point(s) do not occure twice:\n']
            message.append('Point' + '\t'*4 + 'Count\n')
            for point, count in violators:
                message.append(str(point) + '\t' + str(count) + '\n')
            raise Exception(''.join(message))
        
        # Sort all of the points according to Point class' __eq__ __lt__ methods
        eventQ = sorted(self.points, reverse=True)
        
        # The last point in eventQ is the lowest point so save that for later
        self.referencePoint = eventQ[-1]
        # Keep the line which belongs to the reference point
        self.referenceLine = self.referencePoint.parentLine
        # Run lineSweep algorithm to see if there are any extra intersections
        # Line sweep throws the exception if there are any.
        self.lineSweep(eventQ)
        # Sorts the lines in ccw order starting from the reference point
        self.sortedLines = self.orderLines()
        
        # Print the polygon
        self.printShape(self.sortedLines, 'polygon')
        
        # Calculated the convex hull
        self.convexHull = self.createConvexHull()
        
        # Print the convex hull
        self.printShape(self.convexHull, 'convex hull')
        
    def createConvexHull(self):
        """
        Reads the sorted lines of the polygon and creates its convex hull.
        
        Tests each vertex to see if it left or right of the previous line. If
        left add the vertex, if right pop the previous vertex off the list
        and check again.
        
        Parameters
        ----------
        None
        
        Return
        ------
        [lines] - a list of lines which is the convex hull
        
        """
        # Start out with the points from the reference line
        cvPoints = [point for point in self.referenceLine]
        for line in self.sortedLines[1:]:
            # run through the points starting at the line after the ref line
            while len(cvPoints) > 1:
                # The first point is the lowest in X so we know it is on the
                # convex hull so don't pop that one off
                
                # The test line is the last two points in the list
                testLine = Line([cvPoints[-2], cvPoints[-1]])
                # Check the side
                side = self.sideOfLine(testLine, line.end)
                # if it is right or colinear pop it off
                if side != LEFT:
                    cvPoints.pop()
                else: # It was left so exit the loop
                    break
            # append the point to the list
            cvPoints.append(line.end)
        # Turn the points into lines and return them
        return [Line([cvPoints[i], cvPoints[i+1]])
                for i in xrange(0,len(cvPoints)-1)]


    def sideOfLine(self, line, point):
        dist = self.pointToLineDist(line, point)
        if abs(dist) < EPSILON:
            return 0
        return  LEFT if dist < 0 else RIGHT
    
    def pointToLineDist(self, line, point):
        perpVect = np.array([line.vector[Y], -line.vector[X]])
        difPoint = point.pointVector-line.start.pointVector
        return np.dot(perpVect, difPoint)/np.linalg.norm(perpVect)
        
    
    def orderLines(self):
        if(self.referencePoint.x == self.referencePoint.parentLine.end.x):
            self.referencePoint = self.referencePoint.parentLine.end
            self.referenceLine.swapEnds()
        self.lines.remove(self.referenceLine)
        sortedLines = [self.referenceLine]
        while len(self.lines):
            endPoint = sortedLines[-1].end
            for line in self.lines:
                if line.start == endPoint:
                    sortedLines.append(line)
                    self.lines.remove(line)
                    break
                if line.end == endPoint:
                    line.swapEnds()
                    sortedLines.append(line)
                    self.lines.remove(line)
                    break
        return sortedLines

      
    def lineSweep(self, eventQ):
        lineQ = []
        event = eventQ.pop()
        assert(event.pointType == START)
#        lineQ.append(event.parentLine)
        while len(eventQ):            
            if event.pointType == START:
                lineQ.append(event.parentLine)
                if len(lineQ)>1:
                    self.segIntersectTest(lineQ[-1], lineQ[-2])
            else: #There should be no intersections on the eventQ
                index = lineQ.index(event.parentLine)
                lineQ.pop(index)
                if(0 < index < len(lineQ)):
                    self.segIntersectTest(lineQ[index-1], lineQ[index])
            event = eventQ.pop()

    
    def areParallel(self, line1, line2):
        """
        returns True if the two lines are parallel
        
        This method tests if two lines are parallel by finding the angle
        between the perpendicular vector of the first line and the second line.
        If the dot product between perpVect and the vect of line2 is zero then
        line1 and line2 are parallel. Farin and Hansford recommend checking within
        a physically meaningful tolerance so equation 3.14 from pg 50 of
        Farin-Hansford Geometry Toolbox is used to compute the cosine of the angle
        and compare that to our ANGLE_EPS. If cosTheda < ANGLE_EPS then the lines
        are parallel.
        
        Parameters
        ----------
        line1 - the first line
        line2 - the second line
        
        Return
        ------
        True if lines are parallel within ANGLE_EPS else False
        """
        # A vector perpendicular to line1
        perpVect = np.array([-line1.vector[Y], line1.vector[X]])
        # Farin-Hansford eq 3.14
        cosTheda = (np.dot(perpVect, line2.vector)/
                    (np.linalg.norm(perpVect)*np.linalg.norm(line2.vector)))
        # if cosTheda is < ANGLE_EPS then the lines are parallel and we return True
        return abs(cosTheda) < ANGLE_EPS

    def lineIntersect(self, cLine, line):
        if(self.areParallel(cLine, line)):
            return None
            
        return (np.cross((line.start.pointVector - cLine.start.pointVector),
                         line.vector)/(np.cross(cLine.vector, (line.vector))))
                         
    def collinearOverlapping(self, line1, line2):
        distStart = (line1.start - line2.start) + (line1.end - line2.start)
        distEnd = (line1.start - line2.end) + (line1.end - line2.end)
        
        return (abs(distStart - line1.length) < EPSILON or
                abs(distEnd - line1.length) < EPSILON)
        
    def segIntersectTest(self, line1, line2):
        """
        returns the constant aling cLine where the two lines intersect
        
        Given an input of two lines return the constant that when applied to the
        equation of the first line would equal the point of intersetion.
        cLine = p1^ + c*v1^ where the return is the value for c such that the equation
        equals the point of intersection between cLine and otherLine.
        The lines are first tested to see if they are parallel. If they are then None
        is return. From this colinear lines also return None.
        
        Parameters
        ----------
        line1 - the line along which the constant c will be calculated
        otherLine - the line that potentially crosses cLine
        
        Return
        ------
        None - if the two lines are parallel
        float - constant that when applied to the equation of cLine would produce
                the point of intersection.
    
        """
        if(self.areParallel(line1, line2)):
            if self.collinearOverlapping(line1, line2):
                raise Exception('Line intersection error.\nLines:\n' + 
                            str(line1) + '\n' + str(line2) + '\noverlap')
            return None
        
        
        t = self.lineIntersect(line1, line2)
        u = self.lineIntersect(line2, line1)
        if 0 < t < 1 and 0 < u < 1:
            intPoint = Point(line1.start.pointVector + t*line1.vector, INTER, line1)
            raise Exception('Line intersection error.\nLines:\n' + 
                            str(line1) + '\n' + str(line2) + '\nIntersect at '+
                            str(intPoint))
        return None
                

class Line(object):
    
    def __init__(self, line):
        self.start = Point(line[START], START, self)
        self.end = Point(line[END], END, self)
        self.length = self.start - self.end
        self.vector = np.array([self.end.x-self.start.x,
                                self.end.y-self.start.y])
        if self.length < EPSILON:
            raise Exception('Line: ' + str(self) + ' has no length.')
        
    def swapEnds(self):
        self.start, self.end = (Point(self.end, START, self),
                                Point(self.start, END, self))
        self.vector = np.array([self.end.x-self.start.x,
                                self.end.y-self.start.y])
    
    def includedAngle(self, point):
        reverseLine = Line([self[END], self[START]])
        testPoint = point.pointVector-reverseLine.start.pointVector
        return np.dot(reverseLine.vector, testPoint)
    
    def projectOnto(self, other):
        v = self.end.pointVector - other.start.pointVector
        return np.dot(v, other.vector)/np.linalg.norm(other.vector)
        
    def projectPoint(self, point):
        v = point.pointVector - self.start.pointVector
        return np.dot(v, self.vector)/self.length
    
    def __iter__(self):
        yield self.start
        yield self.end
        
    def __getitem__(self, index):
        if index > 1 or index < 0:
            raise IndexError()
        if index == START:
            return self.start
        return self.end
        
    def __repr__(self):
        return self.start.__repr__() +' '+ self.end.__repr__()
                
class Point(object):
    
    def __init__(self, coords, pointType=None, line=None):
        self.pointVector = np.array([round(i, 6) for i in coords])
        self.pointType = pointType
        self.parentLine = line
        self.cosAngle = None
    
    @property
    def x(self):
        return self.pointVector[X]
        
    @property
    def y(self):
        return self.pointVector[Y]
        
    def __iter__(self):
        yield self.x
        yield self.y
    
    def setCosAngle(self, refLine):
        transPointVect = self.pointVector-refLine.start.pointVector
        self.cosAngle = (np.dot(transPointVect,refLine.vector)/
                        (np.linalg.norm(transPointVect)*
                        np.linalg.norm(refLine.vector)))
        
    def __key(self):
        return (self.pointVector[X], self.pointVector[Y])
        
    def __sub__(self, other):
        return np.linalg.norm(self.pointVector - other.pointVector)
    
    def __lt__(self, other):
        return self.__key() < other.__key()
    
    def __eq__(self, other):
        return self.__key() == other.__key()
    
    def __hash__(self):
        return hash(self.__key())
        
    def __repr__(self):
        return 'X%f Y%f'%(self.pointVector[X], self.pointVector[Y])
