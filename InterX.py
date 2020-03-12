# -*- coding: utf-8 -*-
"""
Created on Wed Mar 01 10:00:09 2017

@author: Kelsey N. Lucas

Translation of NS's InterX for Matlab into Python.
Original Available at: https://www.mathworks.com/matlabcentral/fileexchange/22441-curve-intersections/content/InterX.m


Python translation (c) 2017, Kelsey N. Lucas
Original Matlab code Copyright (c) 2009, NS
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.




InterX Intersection of curves

    Input:
    -L1x - Pandas dataframe of x-values for curve 1
    -L1y - Pandas dataframe of y-values for curve 1
    Optional input:
    -L2x - Pandas dataframe of x-values for curve 2
    -L2y - Pandas dataframe of y-values for curve 2
    
    Returns:
    -P - a Pandas dataframe containing two columns - xs, the x-values of 
    intersection points, and ys, the y-values of intersection points.  If no 
    intersection is found, xs and ys are both None.


   P = InterX(L1x,L1y,L2x,L2y) returns the intersection points of two curves L1 
   and L2. The curves L1,L2 can be either closed or open.

   P = InterX(L1x,L1y) returns the self-intersection points of L1. To keep
   the code simple, the points at which the curve is tangent to itself are
   not included.
   

   Author : NS
   Translated from Matlab code Version: 3.0, 21 Sept. 2010
   Translator: Kelsey N. Lucas
   Translation version: 1.0, Mar 2017

   Two words about the algorithm: Most of the code is self-explanatory.
   The only trick lies in the calculation of C1 and C2. To be brief, this
   is essentially the two-dimensional analog of the condition that needs
   to be satisfied by a function F(x) that has a zero in the interval
   [a,b], namely
           F(a)*F(b) <= 0
   C1 and C2 exactly do this for each segment of curves 1 and 2
   respectively. If this condition is satisfied simultaneously for two
   segments then we know that they will cross at some point. 
   Each factor of the 'C' arrays is essentially a matrix containing 
   the numerators of the signed distances between points of one curve
   and line segments of the other.

"""



def InterX(L1x,L1y,L2x=None,L2y=None):
    
    import pandas as pd
    import numpy as np
    
    #Check to see if second curve exists.  If not, set up as identical to first curve.
    if L2x is None:
        L2x = L1x
        L2y = L1y
        
        #hF will tell us how to choose intersection points later - 'lt' means 'less than,'
        #meaning that only intersection and not shared points will be returned
        hF = 'lt' 
        #set up curve 2 as row-vectors
        L2x = np.mat(L2x.as_matrix()).transpose()
        L2y = np.mat(L2y.as_matrix()).transpose()
    else:
        #'le' means less than or equal to - all intersection points are returned
        hF = 'le'
        #set up curve 2 as row-vectors
        L2x = np.mat(L2x.as_matrix()).transpose()
        L2y = np.mat(L2y.as_matrix()).transpose()
    
    #set up curve 1 as column-vectors
    L1x = np.mat(L1x.as_matrix())
    L1y = np.mat(L1y.as_matrix())
    
    #combine x and y for each curve to one master matrix.  Curve 1 is column-vectors, curve 2 is row-vectors
    L1 = np.hstack([L1x,L1y])
    L2 = np.vstack([L2x,L2y])
    
    #break back to x and y values (to match naming conventions from Matlab version)
    #curve 1 is column-vectors, curve 2 is row-vectors
    x1 = L1.transpose()[0].transpose()
    y1 = L1.transpose()[1].transpose()
    x2 = L2[0]
    y2 = L2[1]

    #find the distance between adjacent x's and y's
    dx1 = np.diff(x1,axis=0)
    dx2 = np.diff(x2)
    dy1 = np.diff(y1,axis=0)
    dy2 = np.diff(y2)
    
    #Find 'signed differences'
    S1 = np.multiply(dx1,y1[0:len(y1)-1]) - np.multiply(dy1,x1[0:len(x1)-1])
    S2 = np.multiply(dx2,y2.transpose()[0:len(y2.transpose())-1].transpose()) - np.multiply(dy2,x2.transpose()[0:len(x2.transpose())-1].transpose())
    
    fact1 = dx1*y2 - dy1*x2
    fact2 = y1*dx2 - x1*dy2
    fact2T = fact2.transpose()
    S2T = S2.transpose()

    #Collected distances between points in one curve and line segments in other
    C1 = np.multiply(fact1[0:int(fact1.shape[0]),0:int(fact1.shape[1])-1]-S1,fact1[0:int(fact1.shape[0]),1:int(fact1.shape[1])]-S1)
    C2 = np.multiply(fact2T[0:int(fact2T.shape[0]),0:int(fact2T.shape[1])-1]-S2T,fact2T[0:int(fact2T.shape[0]),1:int(fact2T.shape[1])]-S2T)
    
    #if looking for self-intersections, find only points that aren't tangents between curve 1 and 'curve 2'
    if hF == 'lt':
        TF1 = C1<0
        TF2 = C2<0
        TF2 = TF2.transpose()
    #if looking for intersections between two different curves, take tangent points too
    else:
        TF1 = C1<=0
        TF2 = C2<=0
        TF2 = TF2.transpose()
    #keep indicates row and column indices of line segments where intersections between the two curves are expected
    keep = TF1 & TF2
    
    #collect row and column index values from the keep matrix
    i = []
    j = []
    for row in range(keep.shape[0]):
        for column in range(keep.shape[1]):
            if keep[row,column]:
                i.append(row)
                j.append(column)
    
    #if no intersection is found, return 'None'
    if not i:
        P = [None]
        P = pd.DataFrame(P,columns=['xs'])
        P['ys']=[None]
        return P
    
    #transpose to make data handling easier in a few steps down
    dy2 = dy2.transpose()
    dx2 = dx2.transpose()
    S2 = S2.transpose()

    #calculate some values we need to get output
    L = np.multiply(dy2[j], dx1[i])-np.multiply(dy1[i],dx2[j])

    #L will end up in the denominator a later calculation, so check to make sure none of the values are 0.
    Lnew = []
    inew = []
    jnew = []

    for num in range(L.shape[0]):
        if L[num] != 0:
            Lnew.append(L[num,0])
            inew.append(i[num])
            jnew.append(j[num])
    Lnew = np.mat(Lnew).transpose()

    #Set up numerator and denominator to solve system of equations (two line segments) for intersection points
    numerator = np.hstack([np.multiply(dx2[jnew],S1[inew]) - np.multiply(dx1[inew],S2[jnew]),np.multiply(dy2[jnew],S1[inew]) - np.multiply(dy1[inew],S2[jnew])])
    denominator = np.hstack([Lnew,Lnew])

    #Solve
    result = np.divide(numerator,denominator)

    #organize intersection points into dataframe
    points = pd.DataFrame(result,columns=['xs','ys'])
    points.drop_duplicates(inplace=True)
      
    return points