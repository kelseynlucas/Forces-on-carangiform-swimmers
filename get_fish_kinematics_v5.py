# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 08:26:08 2017

@author: Kelsey N. Lucas



Make sure to load InterX!
"""

import numpy as np
import pandas as pd


#Load some built-in file-reading functions
from os import listdir
from os.path import isfile, join

#Load local minimum /maximum finder function
from scipy.signal import argrelextrema

#Load spline function
from scipy.interpolate import splprep
from scipy.interpolate import splev




def getSequences(directory):
    """
    Reads a file path (directory) - for a test type - and finds the PIV data 
    folders within the directory.   Returns a list of the data folders.
    
    Note: assumes no other file in the directory has a name that starts with 'B'
    
    Input
    
    -directory - file path housing all the sequence folders for an individual fish.
        Note that all '\'s in the default file path 
        must be changed to '/'. No ending slash
        
    """
    
    #If have a subfolder (not a file), get name and add to list of subfolders
    seqs = [f for f in listdir(directory) if not isfile(join(directory,f))]

    #return the dictionary of PIV folders
    return seqs





def getSpeed(seq):
    """
    Look up the swimming speed of the fish in the current sequence
    
    Input:
    -seq - the name of a swimming sequence
    """
    
    #split sequence names on underscore.  Get the 2nd item.
    speed = seq.split('_')[1]
    
    #split the speed string on "B" of BLs - and get first section (speed in characters),
    #and drop units.  Accounts for speed string having different lengths
    speed = speed.split('B')[0]
    
    speed = float(speed)
    
    return speed
    
    
    
    
    

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))
    





def rms(signal):
    """
    calculates the RMS of a signal
    """
    sqs = []    
    for i in range(0,len(signal)):
        sqs.append(signal[i]**2)
     
    return (sum(sqs)/float(len(sqs)))**0.5
    
    
    
    
    
    
def caudalAngle(midX, midY, period, startFr=2):
    """
    Uses midline data to calculate the angle the caudal fin makes with the horizontal
    during period starting at startFr.
    
    Input:
    -midX = x-coordinates of midline
    -midY = y-coordinates of midline
    -period = duration of tailbeat cycle in seconds
    -startFr = frame where movement starts.  Defaults to zero.  Equals frame 
        number of first pectoral fin movement, if applicable.
    """
    #get midline data on last 10 points (the caudal fin) for period of interest
    midX = midX.iloc[90:100,startFr:int(round(period*100))+startFr]
    midY = midY.iloc[90:100,startFr:int(round(period*100))+startFr]

    #initialize storage for angles
    angles = []
    
    #for each video frame,
    stop = int(round(period*100))+startFr
    if stop > max(midY.columns)+1:
        stop = max(midY.columns)+1
        
    for i in range(startFr,stop):
        #get the slope and intercept of best fit line that describes the tail
        m, b = np.polyfit(midX[i],midY[i],1)
        
        #make Pandas dataframes of x and y values on this line
        xvals = [-1000.,1000.]
        yvals = np.multiply(xvals,m)+b

        L1x = pd.DataFrame(xvals,columns=['xvals'])
        L1y = pd.DataFrame(yvals,columns=['yvals'])

        #make Pandas dataframes of a horizontal
        L2x = pd.DataFrame(xvals,columns=['xvals'])
        yvals = [0.,0.]
        L2y = pd.DataFrame(yvals,columns=['yvals'])
        
        try:
            #Find intersection between lines
            vertex = InterX(L1x,L1y,L2x,L2y)
            
            
            #define vectors & Normalize everything by the vertex point coordinates (shift vertex to 0,0)
            #vector A - vertex to (5-norm, L1y corresponding)
            vA = [1000-vertex['xs'][0], L1y['yvals'][1]]
            #vector B - vertex to (5-norm, 0)        
            vB = [1000-vertex['xs'][0], 0.]
            
            #get numerator (dot product of vectors)
            numerator = dotproduct(vA,vB)
                
            #get denominator (magnitude of vA * magnitude of vB)
            #note magnitude = sqrt(dot product of vector with itself)
            denominator = (dotproduct(vA,vA)**0.5)*(dotproduct(vB,vB)**0.5)
            
            #get angle and add to list
            if L1y['yvals'][1] > 0:
                angles.append(np.arccos(numerator/denominator))
            else:
                angles.append(-1.*np.arccos(numerator/denominator))
        
        except:
            angles.append(0.)

    #create a list of time stamps in decimal% of period
    #get number of frames in beat
    if period*100 <= len(angles):
        numFrames = int(round(period*100))
        #make time stamp list
        t = np.multiply(np.divide(range(numFrames),float(numFrames)),100)
    else:
        numFrames = len(angles)
        #make time stamp list
        t = np.linspace(0,float(numFrames)/int(round(period*100)),numFrames)
    
    #assemble into one data frame for save out
    anglesVsTime = pd.DataFrame(t,columns=['time'])
    anglesVsTime['caudalAngle']=np.rad2deg(angles)
    
    return anglesVsTime
    
    

def binnedAngles(midX, midY, period, startFr=2, bins=[8,42,21,15,14]):
    """
    Uses midline data to calculate the angle the segments of the body (bins)
    make with the horizontal during period starting at startFr.
    
    Input:
    -midX = x-coordinates of midline
    -midY = y-coordinates of midline
    -period = duration of tailbeat cycle in seconds
    -startFr = frame where movement starts.  Defaults to zero.  Equals frame 
        number of first pectoral fin movement, if applicable.
    -bins = list of widths of each segment of the body in %BL
    """
    #get indices at the end of each bin
    postEdges = [bins[0]]
    for val in range(1,5):
        postEdges.append(bins[val]+postEdges[val-1]) 
        
    #initialize storage for angles
    angles = np.zeros((int(round(period*100)),5))
    angles = pd.DataFrame(angles,columns=[0,1,2,3,4])
    
    #for each bin,
    for segm in range(len(bins)):
        #get midline data in the bin for period of interest
        if segm == 0:
            binmidX = midX.iloc[0:postEdges[segm],startFr:int(round(period*100))+startFr]
            binmidY = midY.iloc[0:postEdges[segm],startFr:int(round(period*100))+startFr]
        else:
            binmidX = midX.iloc[postEdges[segm-1]:postEdges[segm],startFr:int(round(period*100))+startFr]
            binmidY = midY.iloc[postEdges[segm-1]:postEdges[segm],startFr:int(round(period*100))+startFr]
    
        #for each video frame,
        stop = int(round(period*100))+startFr
        if stop > max(binmidY.columns)+1:
            stop = max(binmidY.columns)+1
        
        for i in range(startFr,stop):
            #get the slope and intercept of best fit line that describes the tail
            m, b = np.polyfit(binmidX[i],binmidY[i],1)
        
            #make Pandas dataframes of x and y values on this line
            xvals = [-1000.,1000.]
            yvals = np.multiply(xvals,m)+b

            L1x = pd.DataFrame(xvals,columns=['xvals'])
            L1y = pd.DataFrame(yvals,columns=['yvals'])

            #make Pandas dataframes of a horizontal
            L2x = pd.DataFrame(xvals,columns=['xvals'])
            yvals = [0.,0.]
            L2y = pd.DataFrame(yvals,columns=['yvals'])
        
            try:
                #Find intersection between lines
                vertex = InterX(L1x,L1y,L2x,L2y)
            
            
                #define vectors & Normalize everything by the vertex point coordinates (shift vertex to 0,0)
                #vector A - vertex to (5-norm, L1y corresponding)
                vA = [1000-vertex['xs'][0], L1y['yvals'][1]]
                #vector B - vertex to (5-norm, 0)        
                vB = [1000-vertex['xs'][0], 0.]
            
                #get numerator (dot product of vectors)
                numerator = dotproduct(vA,vB)
                
                #get denominator (magnitude of vA * magnitude of vB)
                #note magnitude = sqrt(dot product of vector with itself)
                denominator = (dotproduct(vA,vA)**0.5)*(dotproduct(vB,vB)**0.5)
            
                #get angle and add to list
                if L1y['yvals'][1] > 0:
                    angles[segm][i-startFr]=np.rad2deg(np.arccos(numerator/denominator))
                else:
                    angles[segm][i-startFr]=np.rad2deg(-1.*np.arccos(numerator/denominator))
        
            except:
                angles[segm][i-startFr]=0.

    #create a list of time stamps in decimal% of period
    #get number of frames in beat
    if period*100 <= len(angles):
        numFrames = int(round(period*100))
        #make time stamp list
        t = np.multiply(np.divide(range(numFrames),float(numFrames)),100)
    else:
        numFrames = len(angles)
        #make time stamp list
        t = np.linspace(0,float(numFrames)/int(round(period*100)),numFrames)
    
    #assemble into one data frame for save out
    angles['time'] = t
    
    return angles





def kinematics_thru_curvature(midX, midY, startFr=2):
    """
    Calculates the kinematics variables by tracking travel of curvature 
    0-crossing down body
    
    Input:
    -midX - X-coordinates of the midline (pandas dataframe)
    -midY - Y-coordinates of the midline (pandas dataframe)
    -startFr - frame of the video where the sequence begins
    
    """
    #Read in midlines for the sequence
    midX = midX.iloc[:,startFr:len(midX.columns)]
    midY = midY.iloc[:,startFr:len(midY.columns)]
    
    #smooth midline with a quintic least-squares spline

    #first, get arc length value for each x and y coordinate.
    #points on a midline are equidistant, so need distance value for each step
    #Do distance formula
    x0s = midX[midX.index==0].values[0]
    x1s = midX[midX.index==1].values[0]
    xdiff = np.subtract(x1s,x0s)
    xdiffSqed = np.multiply(xdiff,xdiff)
    
    y0s = midY[midY.index==0].values[0]
    y1s = midY[midY.index==1].values[0]
    ydiff = np.subtract(y1s,y0s)
    ydiffSqed = np.multiply(ydiff,ydiff)
    
    distSqed = np.add(xdiffSqed,ydiffSqed)
    dist = np.sqrt(distSqed)
    
    #set up storage for fish length
    fishLen = []
    
    #for each frame (column)
    for col in midY.columns:
        #set up weights where N = 100, num of points on midline, and weight
        #value is 1/N^2 = 0.1
        w = np.full([100,], 0.1)
        
        #define smoothing parameter (when outline was digitized, had 0.5 pix
        #resolution error.  This is ~0.13 mm in all cases.  From smoothing parameter
        #selection, found res error^2/10 gave a good fit)
        s = (2*0.00013**2)/10 
        
        #get arc length value list
        step = dist[col-startFr]
        arcLen = np.multiply(step,range(100))
        #store measure of fish length (for normalizing later)
        fishLen.append(arcLen[99])
        
        #fit spline spline
        xvals = midX.iloc[:,col-startFr].values
        yvals = midY.iloc[:,col-startFr].values
        spl,_ = splprep(u=arcLen, x = [xvals, yvals], w=w, s=s, k=5)
        #evaluate spline
        splXvals, splYvals = splev(arcLen, spl)
        
        #store new x and y vals for curvature calculation
        if col == startFr:
            midXNew = pd.DataFrame(splXvals,columns=[col])
            midYNew = pd.DataFrame(splYvals,columns=[col])
        else:
            midXNew[col] = splXvals
            midYNew[col] = splYvals
    

    #Find x and y spacing between point n and point n+1 on the midline
    dx = pd.DataFrame(np.diff(midXNew,axis=0),columns=midXNew.columns)
    dy = pd.DataFrame(np.diff(midYNew,axis=0),columns=midYNew.columns)
    
    distSqedX = np.multiply(dx,dx)
    distSqedY = np.multiply(dy,dy)
    
    distSqed = np.add(distSqedX,distSqedY)
    
    #Find arc length spacing between point n and point n+1
    ds = pd.DataFrame(np.sqrt(distSqed), columns=midYNew.columns)
    
    #Calculate tangential angles using spacing.  Ang = arctan((dy/ds)/(dx/ds)).
    #Can ignore ds if points on a midline are equally spaced
    tangentialAngles = pd.DataFrame(np.arctan2(np.divide(dy,ds), np.divide(dx,ds)), columns=midYNew.columns)
    
    #get differences between angles (in prep for curvature calculation)
    dtAng = pd.DataFrame(np.diff(tangentialAngles,axis=0),columns=midYNew.columns)
    
    #set up storage for curvatures
    curv = np.zeros_like(midYNew)
    
    #calculate curvature. Because we've done one forward difference already, and 
    #do another here, this is equivalent to a central difference
    curv[1:-1,:] = dtAng /(2*ds[:][:-1])
    curv = pd.DataFrame(curv,columns=midYNew.columns)

    #Get curvatures only for posterior half of body (not a lot of curvature in 
    #anterior half, so this will be the more accurate part)
    curvToSearch = curv.iloc[50:99,:]
    
    
    
    #Will track 0-value of curvature, because more accurate than peak tracking
    #Find 0-crossings
      
    #make storage for list of 0-points
    zeroPts = []
    
    #Look at each column.  Use InterX to find the point of 0 crossing.
    #Need to set up x and y values for each curve for use with InterX.
    #L1 = horizontal line. L2 = curvToSearch values.
    
    L1x = pd.DataFrame([0.,150.],columns=['xvals'])
    L1y = pd.DataFrame([0.,0.],columns=['yvals'])
    
    for col in curvToSearch.columns:
        L2x = np.linspace(50.,98.,len(range(50,99)))
        L2x = pd.DataFrame(L2x,columns=['xvals'])
            
        L2y = pd.DataFrame(curvToSearch[col].values,columns=['yvals'])
            
        zPt = InterX(L1x,L1y,L2x,L2y)
            
        zeroPts.append(zPt['xs'][0])
    
    
    #will need storage for slopes and intercepts of lines fit to position vs time
    #of curvature 0.  Set them up.
    ms = []
    bs = []
    
    
    #find breaks between different wave crests passing. This could happen where
    #either no 0 crossings are found, or where the position of the 0 drops from 
    #near tail to somewhere near mid-body.
    
    #Look for where no 0 crossings were found.
    breaks = [i for i in range(len(zeroPts)) if zeroPts[i] is None]    

    segmStarts = [val+1 for val in breaks]
    segmStarts.insert(0,0)

    if breaks: #Truth value of a list is True if not empty

        #check to make sure there wasn't a sign change in curvature without
        #a no-ID frame in between
        for i in range(len(breaks)+1):
            if i == 0:
                #find spacing between point n and n+1
                spacing = np.roll(zeroPts[0:breaks[i]],-1) - zeroPts[0:breaks[i]]
                spacing = spacing[0:-1]
                
            elif i < len(breaks):
                spacing = np.roll(zeroPts[breaks[i-1]+1:breaks[i]],-1) - zeroPts[breaks[i-1]+1:breaks[i]]
                spacing = spacing[0:-1]
                
            else:
                spacing = np.roll(zeroPts[breaks[i-1]+1:len(zeroPts)],-1) - zeroPts[breaks[i-1]+1:len(zeroPts)]
                spacing = spacing[0:-1]
                
                
                
            #if spacing < 0, means we dropped back on the midline to a new 
            #crossing.  Need to note where happens.
            extraBreaks = [j+1 for j in range(len(spacing)) if spacing[j] < 0]
            
            #Convert indices in extraBreaks from refering to spacing array to
            #the zeroPts array
            
            if i == 0:
                #Then the index already is aligned with index for zeroPts.  Append.
                for num in extraBreaks:
                    segmStarts.append(num)
                    
            else:
                #need to add in number of rows before this segment
                extraBreaks = np.add(breaks[i-1]+1,extraBreaks) 
                
                #add the extra crossing point
                for num in extraBreaks:
                    segmStarts.append(num)

    else:
        #If we didn't have any frames with no 0-crossing IDed, need to search
        #for where the spacing <0 without worrying about breaks
        spacing = np.roll(zeroPts,-1) - zeroPts
        spacing = spacing[0:-1]
        extraBreaks = [j+1 for j in range(len(spacing)) if spacing[j] < -10]
        for num in extraBreaks:
            segmStarts.append(num)
    
    
    #Sort segmStarts
    segmStarts = np.sort(segmStarts)

    #Splice by where segments of curvature tracking start 
    #split up sequence into segments, and fit a line if a group has 5 or more points.
    
    for i in range(1,len(segmStarts)+1):
            
        #reset m and b
        m=0
        b=0
            
        #if we have more than 5 IDs of 0-crossing in a segment, fit line
        if i < len(segmStarts):
            #get the segment (series of BLs where 0-crossing occurs)
            segm = zeroPts[segmStarts[i-1]:segmStarts[i]]
            #set up time stamps for the segment
            t = range(segmStarts[i-1]+startFr,len(segm)+segmStarts[i-1]+startFr)
                
        else:
            segm = zeroPts[segmStarts[i-1]:len(zeroPts)]
            t = range(segmStarts[i-1]+startFr,len(segm)+segmStarts[i-1]+startFr)
            
        
        #drop any values in t corresponding to a None value in segm
        t = [t[j] for j in range(len(segm)) if segm[j] is not None]
        #drop None values from segm
        segm = [val for val in segm if val is not None]
            
        #if there are still 5 or more IDs of 0-crossing, fit a line
        if len(segm) > 4:
            m, b = np.polyfit(t,segm,1)
            
        #If we fit a line, save the slope and intercept
        if not m == 0:
            ms.append(m)
            bs.append(b)
                
       
    #waveSpeed is the average slope
    waveSpeed = np.mean(ms)
    
    #period will be 2x horizontal distance between lines.  We'll calculate 
    #distance at 85% BL
    
    period = np.divide(np.add(85.,np.multiply(bs,-1)),ms)
    period = np.roll(period,-1) - period
    period = period[0:-1]
    period = np.divide(period,100.)    
    period = np.mean(period)*2
    
    #calculate frequency
    freq = 1./period    
    
    #wavelength will be wavespeed x period
    waveLength = waveSpeed * period
    
    
    #get the average fish length
    avgCalcFishLen = np.mean(fishLen[0:int(round(period*100))])
    
    #Now that have period, use smoothed midlines to find tailbeat amplitude and
    #expected index for min and max amplitudes
    tailY = midYNew.T[99][0:int(round(period*100))]
    tailY = tailY - tailY.mean()
    
    #peak-to-peak amplitude will be 2*sqrt(2)*RMS
    tbamp = 2.*np.sqrt(2.)*rms(tailY.values)
    #normalize to lengths
    tbamp = tbamp/avgCalcFishLen
    

    #find the points where the tail crosses the central axis
    
    #use line L1 from earlier
    #make L2 from tail
    L2x = np.linspace(0,int(round(period*100-1)),len(range(int(round(period*100)))))
    L2x = pd.DataFrame(L2x,columns=['xvals'])
    
    L2y = pd.DataFrame(tailY.values,columns=['yvals'])
    
    #use InterX to get crossings
    crossings = InterX(L1x,L1y,L2x,L2y)
    
    #find out if on (-) or (+) side of curve between crossings 0 + 1, and use 
    #to set max & min indexes (based on shape of a sine wave)
    if tailY[startFr] > 0:
        minIdx = crossings['xs'][0] + period*0.25*100
        maxIdx = crossings['xs'][0] + period*0.75*100
    else:
        maxIdx = crossings['xs'][0] + period*0.25*100
        minIdx = crossings['xs'][0] + period*0.75*100
    
     
    
    #make sure max and min values IDed are within the period of interest
    if maxIdx >= period*100:
        maxIdx = maxIdx - period*100
    if minIdx >= period*100:
        minIdx = minIdx - period*100

    #make index values integers
    maxIdx = round(maxIdx)
    minIdx = round(minIdx)   
    
    #convert tangential angles into degrees
    tangentialAnglesDeg = pd.DataFrame(np.rad2deg(tangentialAngles),columns=midY.columns)
    
    return tangentialAnglesDeg, round(waveSpeed,3), round(period,2), round(freq,2), round(waveLength,3), round(tbamp,3), minIdx, maxIdx, round(avgCalcFishLen*100,2)

   


#use this section to analyze fish data
velmsBase = {'klbg2': 0.1, 'klbg3': 0.115, 'klbg4':0.097, 'klbg5':0.096, 'klbg7':0.093, 'trout4': 0.100, 'trout7': 0.110, 'trout8': 0.103}

#pathBase = 'E:/fish pressure_v2/bluegill'
pathBase = 'E:/fish pressure_v2/trout'

#individuals = ['klbg2','klbg3','klbg4','klbg5','klbg7']
individuals = ['trout4','trout7','trout8']



list_waveSpeed_curv = []
list_period_curv = []
list_freq_curv = []
list_waveLength_curv = []
list_tbamp_curv = []
list_tbMaxIdx_curv = []
list_tbMinIdx_curv = []

list_maxCaudAngle = []


list_lengths = []
list_trial = []
list_individs = []
list_speed = []
list_St = []
list_Re = []

list_totNumFrames = []




for individ in individuals:
    individPath = pathBase + '/' + individ 
    seqs = getSequences(individPath)
    
    for seq in seqs:
        
        #seq = 'trout4_2.5BLs_4_ventral_frm683-907'
        if '2.5BLs' not in seq:
            continue
        

        print(seq)
        
        midlineFile = individPath + '/' + seq + '/midlines.xlsx'
        
        midX = pd.read_excel(midlineFile, header=None, sheetname='xvals')
        midY = pd.read_excel(midlineFile, header=None, sheetname='yvals')
        
        
        numFrames = len(midX.columns)

        speed = getSpeed(seq)
        vel = speed*velmsBase[individ]
             
        
        tangentialAnglesDeg, waveSpeed_curv, period_curv, freq_curv, waveLength_curv, tbamp_curv, minIdx_curv, maxIdx_curv, avgCalcFishLen = kinematics_thru_curvature(midX, midY, startFr=2)

        St = (freq_curv*(tbamp_curv*(avgCalcFishLen/100)))/vel
        St = round(St,3)
        
        Re = (1000*vel*(avgCalcFishLen/100))/(0.001002)
        Re_rounded = round(Re/100)*100
        
        
        

        caudalAngleVsTime = caudalAngle(midX, midY, period_curv)
        caudalAnglesPath = individPath + '/' + seq + '/caudalAnglesVsTime.xlsx'
        #caudalAngleVsTime.to_excel(caudalAnglesPath)


        tanAnglesPath = individPath + '/' + seq + '/tanAnglesVsTime.xlsx'
        #tangentialAnglesDeg.to_excel(tanAnglesPath)

            
        maxCaudAngle = max(abs(caudalAngleVsTime['caudalAngle']))
            
        binAnglesVsTime = binnedAngles(midX,midY,period_curv)
        binAnglesPath = individPath + '/' + seq + '/binAnglesVsTime.xlsx'
        #binAnglesVsTime.to_excel(binAnglesPath)
            
        #break
    #break
          

        list_waveSpeed_curv.append(waveSpeed_curv)
        list_period_curv.append(period_curv)
        list_freq_curv.append(freq_curv)
        list_waveLength_curv.append(waveLength_curv)
        list_tbamp_curv.append(tbamp_curv)
        list_tbMaxIdx_curv.append(maxIdx_curv)
        list_tbMinIdx_curv.append(minIdx_curv)

        list_maxCaudAngle.append(maxCaudAngle)


        list_lengths.append(avgCalcFishLen)
        list_trial.append(seq)
        list_individs.append(individ)
        list_speed.append(speed)
        list_St.append(St)
        list_Re.append(Re_rounded)
        
        list_totNumFrames.append(numFrames)

        #break
    
    #break



dfOnePer = pd.DataFrame(list_waveSpeed_curv, columns = ['waveSpeed'])
dfOnePer['period'] = list_period_curv
dfOnePer['freq'] = list_freq_curv
dfOnePer['waveLength'] = list_waveLength_curv
dfOnePer['tbamp'] = list_tbamp_curv
dfOnePer['tbMaxIdx'] = list_tbMaxIdx_curv
dfOnePer['tbMinIdx'] = list_tbMinIdx_curv
dfOnePer['maxCaudAngle'] = list_maxCaudAngle

dfOnePer['fishLen'] = list_lengths
dfOnePer['sequence'] = list_trial
dfOnePer['individual'] = list_individs
dfOnePer['speed'] = list_speed
dfOnePer['St'] = list_St
dfOnePer['Re'] = list_Re

dfOnePer['sourceNumFrames'] = list_totNumFrames

savePath_OnePer = pathBase + '/all_kinematics_data_summary_2.5BLs.xlsx'

dfOnePer.to_excel(savePath_OnePer)









