# -*- coding: utf-8 -*-
"""
Created on Fri Aug 03 11:27:00 2018

@author: Kelsey N. Lucas

Helper functions for assembling fish kinematics, pressure, and forces into a dataframe for stats

Make sure to load InterX!


Note that in original videos, fish swims towards the left side of the screen.
So, forward is -x, backward is +x, left is down, right is up.


Code is written to look up data from a particular folder/file structure.  Edits
may be required based on user's file paths.
Assumed structure:
    master folder for results > folders for results by species > folders for results by individual ID > folders for results by video sequence
Kinematics summary files collect data for a species, and so are housed within
    the results by species folders (on same level as individual folders)

"""


import numpy as np
import pandas as pd


#Load some built-in file-reading functions
from os import listdir
from os.path import isfile, join

from scipy.signal import savgol_filter


#define a helper function to interpolate time-series data to 15 points
def interp_to_15_time_points(period,data):
    tnew = np.divide(range(15),15.)
    told = np.divide(range(int(round(period*100))),round(period*100))
    return np.interp(tnew,told,data)


def getSpeed(seq):
    """
    Look up the swimming speed of the fish in the current sequence
    
    Input:
    -seq - the name of a swimming sequence
    
    Returns:
    -speed - swimming speed in L/s
    """
    
    #split sequence names on underscore.  Get the 2nd item.
    speed = seq.split('_')[1]
    
    #split the speed string on "B" of BLs - and get first section (speed in characters),
    #and drop units.  Accounts for speed string having different lengths
    speed = speed.split('B')[0]
    
    speed = float(speed)
    
    return speed




def get_binned_amp_timeseries_for_seq(midY,length,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets mean amplitude (lateral excursion) of body segments ("bins") vs time
        for a video sequence
    
    Input:
    -midY - dataframe of y-coordinates of midline.  Columns are successive 
        video frames, rows are the list of points from head to tail.
    midY must be cut to period of interest (first column is first frame of 
        interest), and centered at zero
    
    -length - body length of fish in m
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -amps - dataframe where columns are body segments and rows are amplitude 
        over time
    
    '''
    
    #get the end point of each segment in %L
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        postEdges.append(bins[i]+postEdges[i-1])
    
    #print(postEdges)
    
    #convert midline data units from m to lengths
    midYnorm = pd.DataFrame(np.divide(midY,length),columns=midY.columns)
    
    #set up storage for data output
    amps = np.zeros([len(bins),len(midYnorm.columns)])
    
    #for each segment,
    for segm in range(len(postEdges)):
        #if this is the first segment,
        if segm == 0:
            #get the midline coordinates from head to the end of the segment
            midYsegm = midYnorm.iloc[0:postEdges[segm],:]
        #if not,
        else:
            #get the midline coordinates from the end of previous segment to 
            #end of this segment
            midYsegm = midYnorm.iloc[postEdges[segm-1]:postEdges[segm],:]

        #find the mean amplitude in each segment and store
        amps[segm,:] = midYsegm.mean(axis=0).values
    
    #format into a dataframe with columns named for video frame number
    amps = pd.DataFrame(amps,columns=midYnorm.columns)
    
    #transpose dataframe (rows are now video frame, columns are now segment)
    amps = amps.T
    #make sure everything is listed in order
    amps.sort_index(axis=0,inplace=True)
    
    #add a column with the time stamp
    amps['% Period'] = np.linspace(0.,1-1./len(midY.columns),len(midY.columns))
    
    return amps




def get_amp_vs_time_binned(pathBase,seq,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    """
    Calls get_binned_amp_timeseries_for_seq to get amplitude vs time in body
        segments ("bins") for a video sequence, then downsamples to 15 time
        points and applies a savgol filter
        
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -trialbinnedAmps - dataframe where columns are body segments and rows are 
        amplitude over time, downsampled and filtered.
    
    """
   
    #read in the kinematics data summary
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')

    #look up the fish's individual ID
    fishID = seq.split('_')[0]
    #build file path for individual's data
    individPath = pathBase + '/' + fishID
                
    #get row in kinematics dataframe where this sequence's data are
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    #look up fish length and convert from cm to mm
    length = trialMetaData['fishLen'][row]*10

    #look up tailbeat period
    period = trialMetaData['period'][row]

    #there's one column of data per 0.01s, so use period to figure number of 
    #columns of data we expect
    ncol = int(round(period*100))

    #construct file path for midline data
    midlineFile = individPath + '/' + seq + '/midlines.xlsx'

    #read y-coordinates of midline data (lateral movement)
    midY = pd.read_excel(midlineFile, header=None, sheetname='yvals')
    #convert from m to mm
    midY = midY*1000.
    #trim data to period of interest
    midY = midY.iloc[:,startFr:startFr+ncol]
    #look up when the left-most (minimum) tail tip excursion is
    minIdx = trialMetaData['tbMinIdx'][row]

    #will need to reorder data so sequence starts at left-most tail tip excursion
    #set up a dictionary to renumber data columns so left-most tail tip excursion is first
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*minIdx)
    newIdx = range(startFr,startFr+ncol)

    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]

    #get the mean y-coordinate for each point on the midline (should be the trajectory)
    meansToSubtract = midY.mean(axis=1)

    #reorder and center columns so trajectory is at y=0
    for column in midY.columns:
        midY[column] = midY[column]-meansToSubtract

    midY.rename(columns=newcols, inplace=True)
    #now that data are prepped, get the amplitude of body segments vs time
    binnedAmps = get_binned_amp_timeseries_for_seq(midY,length,bins)
    
    #go through each body segment's data
    for col in range(len(bins)):
        #for the first, need to create the dataframe storing downsampled & filtered
        #data in.  For other body segments, append to this dataframe.
        if col == 0:
            #downsample to 15 points
            trialbinnedAmps = pd.DataFrame(interp_to_15_time_points(period,binnedAmps[col]),columns=[col])
        else:
            trialbinnedAmps[col] = interp_to_15_time_points(period,binnedAmps[col])
        #apply savgol filter
        trialbinnedAmps[col] = savgol_filter(trialbinnedAmps[col],11,6)
    #add a time stamp
    trialbinnedAmps['% Period'] = np.divide(range(15),15.)
    

    return trialbinnedAmps
            
        
#next two functions are similar to amplitude functions above, but now for tangential
#angle of midline points
        
def get_binned_ang_timeseries_for_seq(tanAngles,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets mean tangential angle (angle midline makes with the trajectory at a 
        on the midline) of body segments ("bins") vs time for a video sequence
    
    Input:
    -tanAngles - dataframe of tangential angles of midline.  Columns are successive 
        video frames, rows are the list angles at midline points from head to tail.
    tanAngles must be cut to period of interest

    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -angles - dataframe where columns are body segments and rows are tangential angle
        over time
    
    '''
    
    #get the end point of each segment in %L
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        postEdges.append(bins[i]+postEdges[i-1])
    
    #set up storage for data output
    angles = np.zeros([len(bins),len(tanAngles.columns)])
    
    #for each segment,
    for segm in range(len(postEdges)):
        #if this is the first segment,
        if segm == 0:
            #get the tangential angles from head to the end of the segment
            midYsegm = tanAngles.iloc[0:postEdges[segm],:]
        #if this is the last segment,
        elif segm == len(postEdges)-1:
            #get the angles from the end of the previous segment to the end of the dataframe
            midYsegm = tanAngles.iloc[postEdges[segm-1]:99,:]
        #otherwise,
        else:
            #get the angles from the end of the previous segment to the end of this one
            midYsegm = tanAngles.iloc[postEdges[segm-1]:postEdges[segm],:]
        
        #get the mean angle over each segment and store
        angles[segm,:] = midYsegm.mean(axis=0).values
    
    #format into a dataframe with columns named for video frame number
    angles = pd.DataFrame(angles,columns=tanAngles.columns)
    
    #transpose dataframe (rows are now video frame, columns are now segment)
    angles = angles.T
    #make sure everything is listed in order
    angles.sort_index(axis=0,inplace=True)
    
    #add a column with the time stamp
    angles['% Period'] = np.linspace(0.,1-1./len(tanAngles.columns),len(tanAngles.columns))

    
    return angles    
    

def get_angles_vs_time_binned(pathBase,seq,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    """
    Calls get_binned_ang_timeseries_for_seq to get angle vs time in body
        segments ("bins") for a video sequence, then downsamples to 15 time
        points and applies a savgol filter
    
    Functions identically to get_amp_vs_time_binned, except for angles
        
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -trialbinnedAngs - dataframe where columns are body segments and rows are 
        tangential angle over time, downsampled and filtered.
    
    """
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')
    
            
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))

    tanAngleFile = individPath + '/' + seq + '/tanAnglesVsTime.xlsx'
    tanAngles = pd.read_excel(tanAngleFile, header=0)

    tanAngles = tanAngles.iloc[:,0:ncol]

    minIdx = trialMetaData['tbMinIdx'][row]


    oldIdx = np.roll(range(startFr,startFr+ncol),-1*minIdx)
    newIdx = range(startFr,startFr+ncol)

    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    tanAngles.rename(columns=newcols, inplace=True)
    binnedAngles = get_binned_ang_timeseries_for_seq(tanAngles,bins)
    
    for col in range(len(bins)):
        if col == 0:
            trialbinnedAngles = pd.DataFrame(interp_to_15_time_points(period,binnedAngles[col]),columns=[col])
        else:
            trialbinnedAngles[col] = interp_to_15_time_points(period,binnedAngles[col])
        trialbinnedAngles[col] = savgol_filter(trialbinnedAngles[col],11,6)
    trialbinnedAngles['% Period'] = np.divide(range(15),15.)
            
    return trialbinnedAngles
    

    
def calculate_angle_of_attack(amps, length, angles, speed):
    '''
    Calculates the geometric angle of attack vs time for a video sequence.
        Uses amplitude and tangential angle data (from previous functions)
        Based on equation in Read et al. (2003) J Fluids Struct
    
    Input:
    -amps - amplitude vs time in body segments 
        Uses trialbinnedAmps output of get_amp_vs_time_binned
    -length - fish body length in mm
    -angles - tangential angle vs time in body segments
        Uses trialbinnedAngs output of get_ang_vs_time_binned
    -speed - fish swimming speed in m/s
    
    Returns:
    -AoA - dataframe with angles of attack.  Columns are body segments, rows are time.
    '''
    
    #get some reference numbers - the number of time-points and the number of body segments ("bins")
    numPts = len(amps['% Period'])
    numBins = len(amps.columns)-1
    
    #Look up the length of time between data points
    deltaT = amps['% Period'][3]
    #convert fish length from mm to m
    length = length/1000.
    
    #convert amplitudes from L to m
    ampsInM = pd.DataFrame(np.multiply(amps[0],length),columns=[0])
    for col in range(1,numBins):
        ampsInM[col] = np.multiply(amps[col],length)
    
    #make reference lists - the previous (n-1) and next (n+1) amplitude
    nMinus1 = ampsInM.shift(1)
    nPlus1 = ampsInM.shift(-1)
    
    #make storage for derivatives
    derivs = np.zeros([numPts,numBins])
    
    #get central difference derivative where can
    numeratorCen = nPlus1-nMinus1
    derivs[1:numPts-1,:] = np.divide(numeratorCen.iloc[1:numPts-1,:],2*deltaT)
    
    #get forward or backward difference derivative at ends
    numeratorFwd = nPlus1 - ampsInM
    derivs[0,:] = np.divide(numeratorFwd.iloc[0,:],deltaT)
    numeratorBck = ampsInM - nMinus1
    derivs[numPts-1,:] = np.divide(numeratorBck.iloc[numPts-1,:],deltaT)
    
    #divide by speed
    derivs = np.divide(derivs,speed)
    
    #do arctan
    heaveDegrees = np.rad2deg(np.arctan(derivs))
    heaveDegrees = pd.DataFrame(heaveDegrees,columns=range(0,numBins))
    
    #get pitch-heaveDegrees
    AoA = angles - heaveDegrees
    
    return AoA
    
    
def get_AoA_vs_time_binned(pathBase, seq, bins=[10,10,20,15,15,15,15]):
    '''
    Calculates angle of attack for a video sequence from downsampled and filtered
        amplitude and tangential angle dataframes generated using the functions
        above
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    
    Returns:
    -trialbinnedAoA - dataframe with angles of attack.  Columns are body segments, rows are time.
    
    '''
    
    #read in the kinematics data
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')
    #read in the fish's swimming velocity (determined as mean flow velocity upstream of the fish)
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    #find the row where the video sequence's kinematics data are
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    #read the fish's body length, convert from cm to mm
    length = trialMetaData['fishLen'][row]*10
    #look up the swimming speed in m/s for the video
    vel = upstream_vel['upstream_vel'][row]
    
    #get downsampled, filtered amplitudes
    trialbinnedAmps = get_amp_vs_time_binned(pathBase,seq)     

    #get downsampled, filtered tangential angles
    trialbinnedAngs = get_angles_vs_time_binned(pathBase,seq)
                
    #calculate angle of attack vs time
    trialbinnedAoA = calculate_angle_of_attack(trialbinnedAmps, length, trialbinnedAngs, vel)
    
    
    return trialbinnedAoA
    
    
    
    
def get_binned_MEANCp_timeseries_for_seq_TOPSIDE(pressures,distVals,boundx,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets mean coefficient of pressure on body segments ("bins") vs time 
        on the right ("top") side of the body for a video sequence
        
    Functions similarly to get_binned_amp_timeseries_for_seq.
    
    Input:
    -pressures - dataframe of pressure coefficients at query points around the body.  
        Columns are successive video frames, rows are values at different query
        points
    -distVals - dataframe of distance along midline, measured in lengths from the
        snout.  Same structure as pressures.
    -boundx - dataframe of x-coordinates of query points; same structure as pressures
    
    pressures, distVals, & boundx must be cut to period of interest
    
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -Cps - dataframe where columns are body segments and rows are mean pressure
        coefficient over time
        
    
    '''
    
    #Get anterior and posterior boundaries of the body segments in L
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
    #set up storage for coefficients of pressure 
    Cps = np.zeros([len(bins),len(pressures.columns)])
    
    #find the tail tip (max x-coordinate) in each video frame
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    #for each video frame,
    for col in pressures.columns:
        
        #get the pressure and distance values on the right (top side) only
        halfPress = pressures[col][0:maxes[col-startFr]]
        halfDist = distVals[col][0:maxes[col-startFr]]
    
        #collect all the pressure values on each segment
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCps = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCps = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
            
            #find the mean of the pressures on the segment at each point in time
            Cps[segm,col-startFr] = segmCps.mean()
 
    #format into a dataframe with columns named for video frame number
    Cps = pd.DataFrame(Cps, columns=range(startFr,startFr+len(pressures.columns)))
    
    #transpose dataframe (rows are now video frame, columns are now segment)
    Cps = Cps.T
    #add a column with the time stamp
    Cps['% Period'] = np.linspace(0.,1-1./len(pressures.columns),len(pressures.columns))


    return Cps
    
    
def get_binned_MEANCp_timeseries_for_seq_BOTTOMSIDE(pressures,distVals,boundx,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets mean coefficient of pressure on body segments ("bins") vs time 
        on the left ("bottom") side of the body for a video sequence
        
    Functions identically to get_binned_MEANCp_timeseries_for_seq_TOPSIDE, but
        collects pressure data on the left instead of right side.
        
    '''
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    Cps = np.zeros([len(bins),len(pressures.columns)])
    
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    for col in pressures.columns:
        
        halfPress = pressures[col][maxes[col-startFr]:]
        halfDist = distVals[col][maxes[col-startFr]:]
        
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCps = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCps = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
            
            Cps[segm,col-startFr] = segmCps.mean()


    Cps = pd.DataFrame(Cps, columns=range(startFr,startFr+len(pressures.columns)))
        
    Cps = Cps.T
    Cps['% Period'] = np.linspace(0.,1-1./len(pressures.columns),len(pressures.columns))

    
    return Cps
    
    
def get_MEANpressure_vs_time_binned_all_TOP(pathBase, seq, startFr=2, bins=[10,10,20,15,15,15,15]):
    '''
    Calls get_binned_MEANCp_timeseries_for_seq_TOPSIDE to get pressure coefficient
        vs time in body segments ("bins") for a video sequence, then downsamples
        to 15 time points and applies a savgol filter
        
    Functions similarly to get_amp_vs_time_binned
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCps - dataframe with pressure coefficents.  Columns are body segments, rows are time.
    
    '''


    #read in the kinematics data
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    #read in the fish's swimming velocity (determined as mean flow velocity upstream of the fish)
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    #look up the fish's individual ID
    fishID = seq.split('_')[0]
    #build file path for individual's data
    individPath = pathBase + '/' + fishID
    
    #find the row where the video sequence's kinematics data are
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    #look up the swimming speed in m/s for the video
    vel = upstream_vel['upstream_vel'][row]

    #look up the tailbeat period
    period = trialMetaData['period'][row]

    #there's one column of data per 0.01s, so use period to figure number of 
    #columns of data we expect
    ncol = int(round(period*100))

    
    #build file path for sequence's data
    seqPath = individPath + '/' + seq

    #load pressures, normalize to coefficients, and trim dataframe so it only includes
    #data in the tailbeat period of interest
    pressures = pd.read_csv(seqPath + '/pressures.dat',delimiter='\t',header=None)/(0.5*1000*vel*vel)
    pressures = pressures.iloc[:,startFr:startFr+ncol]

    #repeat for distances along midline
    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    #repeat for x-coordinates of query points (& convert units from m to mm)
    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
                
    #will need to reorder data so sequence starts at left-most tail tip excursion
    #set up a dictionary to renumber data columns so left-most tail tip excursion is first
    mIdx = trialMetaData['tbMinIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]

    #do the reordering
    pressuresShift = pressures.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)


    #get pressure coefficient vs time in each segment
    Cps = get_binned_MEANCp_timeseries_for_seq_TOPSIDE(pressuresShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
    
    #downsample to 15 time points
    for col in range(len(bins)):
        if col == 0:
            trialCps = pd.DataFrame(interp_to_15_time_points(period,Cps[col]),columns=[col])
        else:
            trialCps[col] = interp_to_15_time_points(period,Cps[col])
    trialCps['% Period'] = np.divide(range(15),15.)
    
    #apply savgol filter
    for num in range(len(bins)):
        trialCps[num] = savgol_filter(trialCps[num].values,11,6)
    
    return trialCps



def get_MEANpressure_vs_time_binned_all_BOTTOM(pathBase, seq, startFr=2, bins=[10,10,20,15,15,15,15]):
    '''
    Calls get_binned_MEANCp_timeseries_for_seq_BOTTOMSIDE to get pressure coefficient
        vs time in body segments ("bins") for a video sequence, then downsamples
        to 15 time points and applies a savgol filter
        
    Functions identically to get_MEANpressure_vs_time_binned_all_TOP, but for the
        left side ("bottom") of body
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCps - dataframe with pressure coefficents.  Columns are body segments, rows are time.
    
    '''
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))


    seqPath = individPath + '/' + seq

    pressures = pd.read_csv(seqPath + '/pressures.dat',delimiter='\t',header=None)/(0.5*1000*vel*vel)
    pressures = pressures.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
                
                
    mIdx = trialMetaData['tbMinIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    pressuresShift = pressures.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)


                    
    Cps = get_binned_MEANCp_timeseries_for_seq_BOTTOMSIDE(pressuresShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
    
    for col in range(len(bins)):
        if col == 0:
            trialCps = pd.DataFrame(interp_to_15_time_points(period,Cps[col]),columns=[col])
        else:
            trialCps[col] = interp_to_15_time_points(period,Cps[col])
    trialCps['% Period'] = np.divide(range(15),15.)
    
    for num in range(len(bins)):
        trialCps[num] = savgol_filter(trialCps[num].values,11,6)
    
    return trialCps
    

    
def get_MEANpressure_vs_time_binned_all_BOTH(pathBase, seq, startFr=2, bins=[10,10,20,15,15,15,15]):
    '''
    Calls BOTH get_binned_MEANCp_timeseries_for_seq_TOPSIDE and 
        get_binned_MEANCp_timeseries_for_seq_BOTTOMSIDE.
    
    Gets pressure coefficient vs time in body segments ("bins") for a video sequence, 
        then downsamples to 15 time points and applies a savgol filter, and averages
        across left and right side of body.
        
    Left and right side of the body should have identical dynamics, but mirrored.
        To get an average of what happens on a side of the body, we can mirror the
        data from the left, and offset it by 50% of a tailbeat cycle.  This should
        align data with the right side, allowing for an average.
        
        
    Functions similarly to get_MEANpressure_vs_time_binned_all_TOP & _BOTTOM.
    
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCps - dataframe with pressure coefficents.  Columns are body segments, rows are time.
    
    '''
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]


    period = trialMetaData['period'][row]

    ncol = int(round(period*100))


    seqPath = individPath + '/' + seq

    pressures = pd.read_csv(seqPath + '/pressures.dat',delimiter='\t',header=None)/(0.5*1000*vel*vel)
    pressures = pressures.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
                
    #for right (0) and left (1) sides,          
    for side in [0,1]:
        #for the right ("top") side, we want the data starting at max tail tip
        #excursion to the left ("minimum" amplitude), as normal.
        if side == 0:
            mIdx = trialMetaData['tbMinIdx'][row]
        #for the left ("bottom") side, we want the data at a 50% offset, so
        #starting at max tail tip excursion to the RIGHT ("maximum" amplitude)
        elif side == 1:
            mIdx = trialMetaData['tbMaxIdx'][row]

        #reorder data
        oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
        newIdx = range(startFr,startFr+ncol)


        newcols = {}

        for i in range(len(oldIdx)):
            newcols[oldIdx[i]] = newIdx[i]


        pressuresShift = pressures.rename(columns=newcols, inplace=False)
        distValsShift = distVals.rename(columns=newcols, inplace=False)
        boundxShift = boundx.rename(columns=newcols, inplace=False)

        #get the mean pressure coefficient vs time on body segments for the
        #appropriate side of the body
        if side == 0:
            CpsTOP = get_binned_MEANCp_timeseries_for_seq_TOPSIDE(pressuresShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
        elif side == 1:
            CpsBOTTOM = get_binned_MEANCp_timeseries_for_seq_BOTTOMSIDE(pressuresShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
    
    #downsample each to 15 time points
    for col in range(len(bins)):
        if col == 0:
            trialCpsTOP = pd.DataFrame(interp_to_15_time_points(period,CpsTOP[col]),columns=[col])
            trialCpsBOTTOM = pd.DataFrame(interp_to_15_time_points(period,CpsBOTTOM[col]),columns=[col])
        else:
            trialCpsTOP[col] = interp_to_15_time_points(period,CpsTOP[col])
            trialCpsBOTTOM[col] = interp_to_15_time_points(period,CpsBOTTOM[col])
    trialCpsTOP['% Period'] = np.divide(range(15),15.)
    trialCpsBOTTOM['% Period'] = np.divide(range(15),15.)
    
    #apply a savgol filter to each        
    for num in range(len(bins)):
        trialCpsTOP[num] = savgol_filter(trialCpsTOP[num].values,11,6)
        trialCpsBOTTOM[num] = savgol_filter(trialCpsBOTTOM[num].values,11,6)
    
    #get the mean of the two sides
    Cps = (trialCpsTOP+trialCpsBOTTOM)/2.
    
    
    return Cps


def get_MEANCps_peaks(df, peakType):
    
    '''
    Finds the magnitude and timing of peak pressure coefficents
    
    Input:
    -df - dataframe with mean pressure coefficient vs time on body segments,
        such as output of get_MEANpressure_vs_time_binned_all_TOP & _BOTTOM.
        Columns are body segments, rows are time.
    -peakType - 'max' or 'min' - indicates if searching for a max or min value
    
    Returns:
    -peakvals - list of peak magnitude on each body segment
    -peakTimes - list of times of peak magnitude for each body segment
    
    '''
    
    #set up storage for peak magnitudes and times
    peakvals = np.zeros((7,))
    peakvals.fill(np.nan)
    peakTimes = np.zeros((7,))
    peakTimes.fill(np.nan)
    
    #go thru each bin and fill out lists
    if peakType == 'max':
        
        for i in range(0,7):
            peakvals[i] = df.max()[i]
            peakTimes[i] = df['% Period'][df.idxmax()[i]]
        
        #specify offsets of data for sync-ing up times across trials (always 
        #adding offsets of a whole tailbeat cycle)
        #offsets are used when the peak straddles the end of one tailbeat cycle
        #and start of next
        t = peakTimes[0]
        if t > 0.59:
            t -= 1
        peakTimes[0] = t
        
        t = peakTimes[1]
        if t > 0.5:
            t -= 1
        peakTimes[1] = t
        
        t = peakTimes[2]
        if t < 0.25:
            t += 1
        peakTimes[2] = t
        
        t = peakTimes[5]
        if t < 0.5:
            t += 1
        peakTimes[5] = t
        
        t = peakTimes[6]
        if t > 0.6:
            t -= 1
        peakTimes[6] = t
        
        
    elif peakType == 'min':
        
        for i in range(0,7):
            peakvals[i] = df.min()[i]
            peakTimes[i] = df['% Period'][df.idxmin()[i]]
        
        t = peakTimes[2]
        if t > 0.75:
            t -= 1
        peakTimes[2] = t
        
        t = peakTimes[3]
        if t > 0.6:
            t -= 1
        peakTimes[3] = t
        
        t = peakTimes[4]
        if t > 0.6:
            t -= 1
        peakTimes[4] = t        
    

    
    return peakvals, peakTimes




    
    
def get_binned_NETCF_timeseries_for_seq_TOPSIDE(forces,distVals,boundx,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets NET coefficient of force on body segments ("bins") vs time 
        on the right ("top") side of the body for a video sequence
        
    Functions similarly to get_binned_MEANCp_timeseries_for_seq_TOPSIDE.
    
    Input:
    -forces - dataframe of force coefficients at query points around the body.  
        Columns are successive video frames, rows are values at different query
        points
    -distVals - dataframe of distance along midline, measured in lengths from the
        snout.  Same structure as forces.
    -boundx - dataframe of x-coordinates of query points; same structure as forces
    
    forces, distVals, & boundx must be cut to period of interest
    
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -CFs - dataframe where columns are body segments and rows are net force
        coefficient over time    
    
    '''
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    CFs = np.zeros([len(bins),len(forces.columns)])
    
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    for col in forces.columns:
        
        halfForces = forces[col][0:maxes[col-startFr]]
        halfDist = distVals[col][0:maxes[col-startFr]]
    
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
            
            #notice difference from get_binned_MEANCp_timeseries_for_seq_TOPSIDE - net (sum), not mean
            CFs[segm,col-startFr] = segmCFs.sum()
 
    CFs = pd.DataFrame(CFs, columns=range(startFr,startFr+len(forces.columns)))
        
    CFs = CFs.T
    CFs['% Period'] = np.linspace(0.,1-1./len(forces.columns),len(forces.columns))

    
    return CFs
    
def get_binned_NETCF_timeseries_for_seq_BOTTOMSIDE(forces,distVals,boundx,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets NET coefficient of force on body segments ("bins") vs time 
        on the left ("bottom") side of the body for a video sequence
        
    Functions identically to get_binned_NETCF_timeseries_for_seq_TOPSIDE, but 
        for left side of body
        
    '''
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    CFs = np.zeros([len(bins),len(forces.columns)])
    
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    for col in forces.columns:
        
        halfForces = forces[col][maxes[col-startFr]:]
        halfDist = distVals[col][maxes[col-startFr]:]
        
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
            
            CFs[segm,col-startFr] = segmCFs.sum()


    CFs = pd.DataFrame(CFs, columns=range(startFr,startFr+len(forces.columns)))
        
    CFs = CFs.T
    CFs['% Period'] = np.linspace(0.,1-1./len(forces.columns),len(forces.columns))

    
    return CFs
    
    
def get_binned_NETCF_timeseries_for_seq_SIDESADDED(forces,distVals,boundx,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets NET coefficient of force on body segments ("bins") vs time summed across
        left & right sides of the body for a video sequence
        
    Functions identically to get_binned_NETCF_timeseries_for_seq_TOPSIDE & _BOTTOMSIDE,
        but doesn't divide data up by side
        
    '''
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    CFs = np.zeros([len(bins),len(forces.columns)])
    
    
    for col in forces.columns:
        
        allForces = forces[col]
        allDist = distVals[col]
    
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCFs = allForces[(allDist>=anteEdges[segm]) & (allDist<postEdges[segm])]
            else:
                segmCFs = allForces[(allDist>=anteEdges[segm]) & (allDist<=postEdges[segm])]
                
            CFs[segm,col-startFr] = segmCFs.sum()
 
    CFs = pd.DataFrame(CFs, columns=range(startFr,startFr+len(forces.columns)))
        
    CFs = CFs.T
    CFs['% Period'] = np.linspace(0.,1-1./len(forces.columns),len(forces.columns))

    
    return CFs


def get_NETforce_vs_time_wholebody_all_SIDESADDED(pathBase, seq, forceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls get_binned_NETCF_timeseries_for_seq_SIDESADDED to get force coefficient
        vs time in body segments ("bins") for a video sequence, then downsamples
        to 15 time points and applies a savgol filter
        
    Functions similarly to get_MEANpressure_vs_time_binned_all_TOP
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -forceType - 'Fx' or 'Fy' - specify if looking at axial or lateral forces
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    #copied from kinematics data - lateral surface areas of different individuals
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
    
   

    #read in the kinematics data
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    #read in the fish's swimming velocity (determined as mean flow velocity upstream of the fish)
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    #look up the fish's individual ID
    fishID = seq.split('_')[0]
    #build file path for individual's data
    individPath = pathBase + '/' + fishID
    #look up the lateral surface area for the individual
    latArea = latAreasm2[fishID]
    
    #find the row where the video sequence's kinematics data are
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    #look up the swimming speed in m/s for the video
    vel = upstream_vel['upstream_vel'][row]

    #look up the tailbeat period
    period = trialMetaData['period'][row]
    
    #there's one column of data per 0.01s, so use period to figure number of 
    #columns of data we expect
    ncol = int(round(period*100))

    #build file path for sequence's data
    seqPath = individPath + '/' + seq
    
    #read in force data.  Because x-axis has positive to the right, and
    #fish swims to the left side of screen, flip sign of x-forces so thrust is 
    #positive
    try:
        if forceType == 'Fx':
            forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
            forces = forces*-1
        elif forceType == 'Fy':
            forces = pd.read_csv(seqPath + '/yForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
    except:
        print('forceType must be Fx or Fy')
                    
    #trim force data down to the period of interest
    forces = forces.iloc[:,startFr:startFr+ncol]

    #read in distance values and trim to period of interest
    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    #repeat for x-coordinates of query points (and convert from m to mm)
    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
    #reorder data so sequence starts at left-most tail tip excursion
    mIdx = trialMetaData['tbMinIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    forcesShift = forces.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)
    
    #get force coefficients added across the two sides of the body
    CFs = get_binned_NETCF_timeseries_for_seq_SIDESADDED(forcesShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
           
    #downsample to 15 time points
    for col in range(len(bins)):
        if col == 0:
            trialCFs = pd.DataFrame(interp_to_15_time_points(period,CFs[col]),columns=[col])
        else:
            trialCFs[col] = interp_to_15_time_points(period,CFs[col])
    trialCFs['% Period'] = np.divide(range(15),15.)
    
    #apply a savgol filter
    for num in range(len(bins)):
        trialCFs[num] = savgol_filter(trialCFs[num].values,11,4)
    
    #add a time stamp
    totalCFs = pd.DataFrame(trialCFs['% Period'].values, columns = ['% Period'])
    #get the net force coefficient across all of the body segments
    totalCFs['net'] = trialCFs.drop(['% Period'], axis=1).sum(axis=1)
    
    
    return totalCFs

    
    
def get_NETforce_vs_time_binned_all_TOP(pathBase, seq, forceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls get_binned_NETCF_timeseries_for_seq_TOPSIDE to get force coefficient
        vs time in body segments ("bins") for a video sequence, then downsamples
        to 15 time points and applies a savgol filter, for data on the right
        ("top") side of the body
        
    Functions similarly to get_MEANpressure_vs_time_binned_all_TOP & 
        get_NETforce_vs_time_wholebody_all_SIDESADDED
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -forceType - 'Fx' or 'Fy' - specify if looking at axial or lateral forces
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
    
    
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
    latArea = latAreasm2[fishID]
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))



    seqPath = individPath + '/' + seq
    
    #must multiply x-forces by -1 because fish swims toward left but x-axis runs toward right
    #so flipping sign makes thrust forces have a (+) sign
    try:
        if forceType == 'Fx':
            forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
            forces = forces*-1
        elif forceType == 'Fy':
            forces = pd.read_csv(seqPath + '/yForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
    except:
        print('forceType must be Fx or Fy')
                    
                
    forces = forces.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
                
                
    mIdx = trialMetaData['tbMinIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    forcesShift = forces.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)



    CFs = get_binned_NETCF_timeseries_for_seq_TOPSIDE(forcesShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
                
    for col in range(len(bins)):
        if col == 0:
            trialCFs = pd.DataFrame(interp_to_15_time_points(period,CFs[col]),columns=[col])
        else:
            trialCFs[col] = interp_to_15_time_points(period,CFs[col])
    trialCFs['% Period'] = np.divide(range(15),15.)
    
            
    for num in range(len(bins)):
        trialCFs[num] = savgol_filter(trialCFs[num].values,11,6)
    
    
    return trialCFs
    
    
    
def get_NETforce_vs_time_binned_all_BOTTOM(pathBase, seq, forceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls get_binned_NETCF_timeseries_for_seq_BOTTOMSIDE to get force coefficient
        vs time in body segments ("bins") for a video sequence, then downsamples
        to 15 time points and applies a savgol filter, for data on the left
        ("bottom") side of the body
        
    Functions similarly to get_NETforce_vs_time_binned_all_TOP
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -forceType - 'Fx' or 'Fy' - specify if looking at axial or lateral forces
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
    
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
    latArea = latAreasm2[fishID]
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))



    seqPath = individPath + '/' + seq
                
    try:
        if forceType == 'Fx':
            forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
            forces = forces*-1
        elif forceType == 'Fy':
            forces = pd.read_csv(seqPath + '/yForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
            forces = forces*-1
    except:
        print('forceType must be Fx or Fy')
                    
                
    forces = forces.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
                
    #notice different index for reordering - reorder so tail tip is at maximum 
    #excursion to the right ("maximum" amplitude), rather than left.  Done so
    #that data can be mirrored and synced with data from right side of the body
    #later to get an average
    mIdx = trialMetaData['tbMaxIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    forcesShift = forces.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)


                    
    CFs = get_binned_NETCF_timeseries_for_seq_BOTTOMSIDE(forcesShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
    
    for col in range(len(bins)):
        if col == 0:
            trialCFs = pd.DataFrame(interp_to_15_time_points(period,CFs[col]),columns=[col])
        else:
            trialCFs[col] = interp_to_15_time_points(period,CFs[col])
    trialCFs['% Period'] = np.divide(range(15),15.)
    
            
    for num in range(len(bins)):
        trialCFs[num] = savgol_filter(trialCFs[num].values,11,6)
    
    
    return trialCFs
    
    
    
    
def get_NETforce_vs_time_binned_all_BOTH(pathBase, seq, forceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls BOTH get_binned_NETCF_timeseries_for_seq_TOPSIDE and 
        get_binned_NETCF_timeseries_for_seq_BOTTOMSIDE.
    
    Gets force coefficient vs time in body segments ("bins") for a video sequence, 
        then downsamples to 15 time points and applies a savgol filter, and averages
        across left and right side of body.
        
    Left and right side of the body should have identical dynamics, but mirrored.
        To get an average of what happens on a side of the body, we can mirror the
        data from the left, and offset it by 50% of a tailbeat cycle.  This should
        align data with the right side, allowing for an average.
        
        
    Functions similarly to get_MEANpressure_vs_time_binned_all_BOTH &
        get_NETforce_vs_time_binned_all_TOP & _BOTTOM.
    
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -forceType - 'Fx' or 'Fy' - specify if looking at axial or lateral forces
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
    

    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
    latArea = latAreasm2[fishID]
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))



    seqPath = individPath + '/' + seq

    try:
        if forceType == 'Fx':
            forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
            forces = forces*-1
        elif forceType == 'Fy':
            forces = pd.read_csv(seqPath + '/yForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
    except:
        print('forceType must be Fx or Fy')
                    
                
    forces = forces.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
            
            
    for side in [0,1]:
        if side == 0:
            mIdx = trialMetaData['tbMinIdx'][row]
        elif side == 1:
            mIdx = trialMetaData['tbMaxIdx'][row]
        if forceType == 'Fy':
            forces = forces*-1


        oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
        newIdx = range(startFr,startFr+ncol)


        newcols = {}

        for i in range(len(oldIdx)):
            newcols[oldIdx[i]] = newIdx[i]


        forcesShift = forces.rename(columns=newcols, inplace=False)
        distValsShift = distVals.rename(columns=newcols, inplace=False)
        boundxShift = boundx.rename(columns=newcols, inplace=False)


        if side == 0:
            CFsTOP = get_binned_NETCF_timeseries_for_seq_TOPSIDE(forcesShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
        elif side == 1:
            CFsBOTTOM = get_binned_NETCF_timeseries_for_seq_BOTTOMSIDE(forcesShift,distValsShift,boundxShift,startFr=startFr,bins=bins)
            
    for col in range(len(bins)):
        if col == 0:
            trialCFsTOP = pd.DataFrame(interp_to_15_time_points(period,CFsTOP[col]),columns=[col])
            trialCFsBOTTOM = pd.DataFrame(interp_to_15_time_points(period,CFsBOTTOM[col]),columns=[col])
        else:
            trialCFsTOP[col] = interp_to_15_time_points(period,CFsTOP[col])
            trialCFsBOTTOM[col] = interp_to_15_time_points(period,CFsBOTTOM[col])
    trialCFsTOP['% Period'] = np.divide(range(15),15.)
    trialCFsBOTTOM['% Period'] = np.divide(range(15),15.)
    
            
    for num in range(len(bins)):
        trialCFsTOP[num] = savgol_filter(trialCFsTOP[num].values,11,6)
        trialCFsBOTTOM[num] = savgol_filter(trialCFsBOTTOM[num].values,11,6)
        
    CFs = (trialCFsTOP+trialCFsBOTTOM)/2.
    
    return CFs
    
    
def get_binned_CFx_bysource_for_seq_TOPSIDE(sourceType,forces,distVals,boundx,xNorms,pressures,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets NET coefficient of force on body segments ("bins") vs time 
        on the right ("top") side of the body for a video sequence.
    Works to get a specific type of axial force coefficients:
        -fpush - positive pressure thrust (forward push)
        -fpull - negative pressure thrust (forward pull)
        -rpush - positive pressure drag (reverse push)
        -rpull - negative pressure drag (reverse pull)
        
    Functions similarly to get_binned_MEANCp_timeseries_for_seq_TOPSIDE (etc).
    
    Input:
    -sourceType - specifies what type of axial forces to tally
        options: fpush, fpull, rpush, rpull
    -forces - dataframe of force coefficients at query points around the body.  
        Columns are successive video frames, rows are values at different query
        points
    -distVals - dataframe of distance along midline, measured in lengths from the
        snout.  Same structure as forces.
    -boundx - dataframe of x-coordinates of query points; same structure as forces
    -xNorms - x-component of the outward-facing normal vector at query points;
        Same structure as forces
    -pressures - pressure coefficients at query points.  Same structure as forces
    
    forces, distVals, boundx, xNorms, & pressures must be cut to period of interest
    
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -CFs - dataframe where columns are body segments and rows are net force
        coefficient over time
    
    '''
    
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    CFs = np.zeros([len(bins),len(forces.columns)])
    
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    for col in forces.columns:
        
        #start like similar functions and get all values in a body segment
        
        halfForces = forces[col][0:maxes[col-startFr]]
        halfDist = distVals[col][0:maxes[col-startFr]]
        halfxNorm = xNorms[col][0:maxes[col-startFr]]
        halfPress = pressures[col][0:maxes[col-startFr]]
        
        
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
                segmxNorm = halfxNorm[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
                segmPress = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                segmxNorm = halfxNorm[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                segmPress = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                
            #use pressure and normal vector to identify axial force type we're interested;
            #only tally forces of this type
            
            #fpull: x-component of normal vector points left (direction of swimming) 
            #and pressure is negative (directing force out of surface)
            if sourceType == 'fpull':
                segmCFs = segmCFs[(segmxNorm<0)&(segmPress<0)]
            #rpush: x-component of normal vector points left (direction of swimming) 
            #and pressure is positive (directing force into surface)
            elif sourceType == 'rpush':
                segmCFs = segmCFs[(segmxNorm<0)&(segmPress>0)]
            #fpush: x-component of normal vector points right (opposite of swimming direction)
            #and pressure is positive (directing force into surface)
            elif sourceType == 'fpush':
                segmCFs = segmCFs[(segmxNorm>0)&(segmPress>0)]
            #fpush: x-component of normal vector points right (opposite of swimming direction) 
            #and pressure is negative (directing force out of surface)
            elif sourceType == 'rpull':
                segmCFs = segmCFs[(segmxNorm>0)&(segmPress<0)]
                
            CFs[segm,col-startFr] = segmCFs.sum()
 
    CFs = pd.DataFrame(CFs, columns=range(startFr,startFr+len(forces.columns)))
        
    CFs = CFs.T
    CFs['% Period'] = np.linspace(0.,1-1./len(forces.columns),len(forces.columns))
    
    return CFs
    
    
    
    
def get_binned_CFx_bysource_for_seq_BOTTOMSIDE(sourceType,forces,distVals,boundx,xNorms,pressures,startFr=2,bins=[10,10,20,15,15,15,15]):
    
    '''
    Gets NET coefficient of force on body segments ("bins") vs time 
        on the left ("bottom") side of the body for a video sequence.
    Works to get a specific type of axial force coefficients:
        -fpush - positive pressure thrust (forward push)
        -fpull - negative pressure thrust (forward pull)
        -rpush - positive pressure drag (reverse push)
        -rpull - negative pressure drag (reverse pull)
        
    Functions identically to get_binned_CFx_bysource_for_seq_TOPSIDE, but for 
        left side of body.
    
    Input:
    -sourceType - specifies what type of axial forces to tally
        options: fpush, fpull, rpush, rpull
    -forces - dataframe of force coefficients at query points around the body.  
        Columns are successive video frames, rows are values at different query
        points
    -distVals - dataframe of distance along midline, measured in lengths from the
        snout.  Same structure as forces.
    -boundx - dataframe of x-coordinates of query points; same structure as forces
    -xNorms - x-component of the outward-facing normal vector at query points;
        Same structure as forces
    -pressures - pressure coefficients at query points.  Same structure as forces
    
    forces, distVals, boundx, xNorms, & pressures must be cut to period of interest
    
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
        
    Returns:
    -CFs - dataframe where columns are body segments and rows are net force
        coefficient over time
    
    '''
    
    anteEdges = [0]
    postEdges = [bins[0]]
    for i in range(1,len(bins)):
        anteEdges.append(bins[i-1]+anteEdges[i-1])
        postEdges.append(bins[i]+postEdges[i-1])
    anteEdges = np.divide(anteEdges,100.)
    postEdges = np.divide(postEdges,100.)
        
        
    CFs = np.zeros([len(bins),len(forces.columns)])
    
    maxes = np.add(boundx.idxmax(axis=0).values,1)
    
    for col in forces.columns:
        
        halfForces = forces[col][maxes[col-startFr]:]
        halfDist = distVals[col][maxes[col-startFr]:]
        halfxNorm = xNorms[col][maxes[col-startFr]:]
        halfPress = pressures[col][maxes[col-startFr]:]
        
        for segm in range(len(postEdges)):
            if not(segm == len(postEdges)-1):
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
                segmxNorm = halfxNorm[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
                segmPress = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<postEdges[segm])]
            else:
                segmCFs = halfForces[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                segmxNorm = halfxNorm[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                segmPress = halfPress[(halfDist>=anteEdges[segm]) & (halfDist<=postEdges[segm])]
                
                
            if sourceType == 'fpull':
                segmCFs = segmCFs[(segmxNorm<0)&(segmPress<0)]
            elif sourceType == 'rpush':
                segmCFs = segmCFs[(segmxNorm<0)&(segmPress>0)]
            elif sourceType == 'fpush':
                segmCFs = segmCFs[(segmxNorm>0)&(segmPress>0)]
            elif sourceType == 'rpull':
                segmCFs = segmCFs[(segmxNorm>0)&(segmPress<0)]
            
            CFs[segm,col-startFr] = segmCFs.sum()


    CFs = pd.DataFrame(CFs, columns=range(startFr,startFr+len(forces.columns)))
        
    CFs = CFs.T
    CFs['% Period'] = np.linspace(0.,1-1./len(forces.columns),len(forces.columns))
    
    return CFs



def get_NETCFx_sourced_vs_time_binned_all_TOP(pathBase, seq, sourceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls get_binned_CFx_bysource_for_seq_TOPSIDE to get force coefficient vs 
        time in body segments ("bins") for a video sequence, then downsamples 
        to 15 time points and applies a savgol filter.  Collects force data for
        right ("top") side of body
    Works to get a specific type of axial force coefficients:
        -fpush - positive pressure thrust (forward push)
        -fpull - negative pressure thrust (forward pull)
        -rpush - positive pressure drag (reverse push)
        -rpull - negative pressure drag (reverse pull)
        
    Functions similarly to get_MEANpressure_vs_time_binned_all_TOP &
        get_NETforce_vs_time_binned_all_TOP & _BOTTOM.
    
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -sourceType - specifies what type of axial forces to tally
        options: fpush, fpull, rpush, rpull
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
    
            
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
    latArea = latAreasm2[fishID]
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))

    seqPath = individPath + '/' + seq
                
                
    forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
    forces = forces*-1
                
                
    forces = forces.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
    
    #xUnitNormals.dat gives inward-facing normal vector.  Multiply by -1 to have
    #direction of outward-facing normal vector.
    xNorms = pd.read_csv(seqPath + '/xUnitNormals.dat',delimiter='\t',header=None)*-1
    xNorms = xNorms.iloc[:,startFr:startFr+ncol]
                
    press = pd.read_csv(seqPath + '/pressures.dat',delimiter='\t',header=None)
    press = press.iloc[:,startFr:startFr+ncol]
            

    mIdx = trialMetaData['tbMinIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    forcesShift = forces.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)
    xNormsShift = xNorms.rename(columns=newcols, inplace=False)
    pressShift = press.rename(columns=newcols, inplace=False)


                
    CFs = get_binned_CFx_bysource_for_seq_TOPSIDE(sourceType,forcesShift,distValsShift,boundxShift,xNormsShift,pressShift,startFr=startFr,bins=bins)
              
    for col in range(len(bins)):
        if col == 0:
            trialCFs = pd.DataFrame(interp_to_15_time_points(period,CFs[col]),columns=[col])
        else:
            trialCFs[col] = interp_to_15_time_points(period,CFs[col])
    trialCFs['% Period'] = np.divide(range(15),15.)
    
            
    for num in range(len(bins)):
        trialCFs[num] = savgol_filter(trialCFs[num].values,11,6)
        
    return trialCFs
    

def get_NETCFx_sourced_vs_time_binned_all_BOTTOM(pathBase, seq, sourceType, startFr=2, bins=[10,10,20,15,15,15,15]):
    
    '''
    Calls get_binned_CFx_bysource_for_seq_BOTTOMSIDE to get force coefficient vs 
        time in body segments ("bins") for a video sequence, then downsamples 
        to 15 time points and applies a savgol filter.  Collects force data for
        left ("bottom") side of body
    Works to get a specific type of axial force coefficients:
        -fpush - positive pressure thrust (forward push)
        -fpull - negative pressure thrust (forward pull)
        -rpush - positive pressure drag (reverse push)
        -rpull - negative pressure drag (reverse pull)
        
    Functions identically to get_NETCFx_sourced_vs_time_binned_all_TOP but for
        left ("bottom") side of body
    
    
    Input:
    -pathBase - folder where the kinematics data summary table is found
    -seq - name of the video sequence
    -sourceType - specifies what type of axial forces to tally
        options: fpush, fpull, rpush, rpull
    -startFr - the frame of video where the tailbeat cycle of interest starts 
        (default is 2 - note Python indexes from 0, so skips first 2 frames)
    -bins - how long the body segments are.  Default lengths are 10, 10, 20, 
        15, 15, 15, and 15%.
    
    Returns:
    -trialCFs - dataframe with force coefficents.  Columns are body segments, rows are time.
    
    '''
    
    latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,
                 'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
  
    
    trialMetaData = pd.read_excel(pathBase+'/all_kinematics_data_summary.xlsx',header=0,index_col=0,sheetname='Sheet1')        
    upstream_vel = pd.read_excel(pathBase+'/upstream_velocity_from_PIV.xlsx',header=0,index_col=0)        
    
    fishID = seq.split('_')[0]
    individPath = pathBase + '/' + fishID
    latArea = latAreasm2[fishID]
                
    row = trialMetaData['sequence'].loc[trialMetaData['sequence']==seq].index[0]
    vel = upstream_vel['upstream_vel'][row]

    period = trialMetaData['period'][row]

    ncol = int(round(period*100))

    seqPath = individPath + '/' + seq
                
    forces = pd.read_csv(seqPath + '/xForcesAll.dat',delimiter='\t',header=None)/(0.5*1000*latArea*vel*vel)
    forces = forces*-1
                
                
    forces = forces.iloc[:,startFr:startFr+ncol]

    distVals = pd.read_csv(seqPath + '/distOnMidline.dat',delimiter='\t',header=None)
    distVals = distVals.iloc[:,startFr:startFr+ncol]

    boundx = pd.read_csv(seqPath + '/xvalsBound.dat',delimiter='\t',header=None)*1000
    boundx = boundx.iloc[:,startFr:startFr+ncol]
                
    xNorms = pd.read_csv(seqPath + '/xUnitNormals.dat',delimiter='\t',header=None)*-1
    xNorms = xNorms.iloc[:,startFr:startFr+ncol]
                
    press = pd.read_csv(seqPath + '/pressures.dat',delimiter='\t',header=None)
    press = press.iloc[:,startFr:startFr+ncol]
    
                
    mIdx = trialMetaData['tbMaxIdx'][row]
                    
    oldIdx = np.roll(range(startFr,startFr+ncol),-1*mIdx)
    newIdx = range(startFr,startFr+ncol)


    newcols = {}

    for i in range(len(oldIdx)):
        newcols[oldIdx[i]] = newIdx[i]


    forcesShift = forces.rename(columns=newcols, inplace=False)
    distValsShift = distVals.rename(columns=newcols, inplace=False)
    boundxShift = boundx.rename(columns=newcols, inplace=False)
    xNormsShift = xNorms.rename(columns=newcols, inplace=False)
    pressShift = press.rename(columns=newcols, inplace=False)


    CFs = get_binned_CFx_bysource_for_seq_BOTTOMSIDE(sourceType,forcesShift,distValsShift,boundxShift,xNormsShift,pressShift,startFr=startFr,bins=bins)
    
    for col in range(len(bins)):
        if col == 0:
            trialCFs = pd.DataFrame(interp_to_15_time_points(period,CFs[col]),columns=[col])
        else:
            trialCFs[col] = interp_to_15_time_points(period,CFs[col])
    trialCFs['% Period'] = np.divide(range(15),15.)
    
            
    for num in range(len(bins)):
        trialCFs[num] = savgol_filter(trialCFs[num].values,11,6)
        
    return trialCFs
             


def get_NETCFx_source_peaks(df, sourceType):
    
    '''
    Finds the magnitude and timing of peak axial force coefficents of a particular type
    Works to get a specific type of axial force coefficients:
        -fpush - positive pressure thrust (forward push)
        -fpull - negative pressure thrust (forward pull)
        -rpush - positive pressure drag (reverse push)
        -rpull - negative pressure drag (reverse pull)
    
    Input:
    -df - dataframe with mean axial coefficient vs time on body segments,
        such as output of get_NETCFx_sourced_vs_time_binned_all_TOP & _BOTTOM.
        Columns are body segments, rows are time.
    -sourceType - specifies what type of axial forces to tally
        options: fpush, fpull, rpush, rpull
        
    Returns:
    -peakvals - list of peak magnitude on each body segment
    -peakTimes - list of times of peak magnitude on each body segment
    -peakIdx - list of row indices where peak occurs in df for each body segment
    
    '''
    
    peakvals = np.zeros((7,))
    peakvals.fill(np.nan)
    peakTimes = np.zeros((7,))
    peakTimes.fill(np.nan)
    peakIdx = np.zeros((7,))
    peakIdx.fill(np.nan)
    
    #go thru each bin and fill out lists
    if sourceType == 'fpush':
        
        for i in range(2,7):
            peakvals[i] = df.max()[i]
            peakIdx[i] = df.idxmax()[i]
        
        peakTimes[2] = df['% Period'][df.idxmax()[2]]
        peakTimes[3] = df['% Period'][df.idxmax()[3]]
        
        #offset data for sync-ing purposes (notice all offsets are a full tailbeat cycle)
        #offsets are used when the peak straddles the end of one tailbeat cycle
        #and start of next
        t = df['% Period'][df.idxmax()[4]]
        if t < df['% Period'][2]:
            t += 1
        peakTimes[4] = t
        
        t = df['% Period'][df.idxmax()[5]]
        if t < df['% Period'][6]:
            t += 1
        peakTimes[5] = t
        
        t = df['% Period'][df.idxmax()[6]]
        if t >= df['% Period'][9]:
            t -= 1
        peakTimes[6] = t
        
        
    elif sourceType == 'fpull':
        
        for i in range(2,7):
            peakvals[i] = df.max()[i]
            peakIdx[i] = df.idxmax()[i]
            
        
        t = df['% Period'][df.idxmax()[2]]
        if t < df['% Period'][6]:
            t += 1
        peakTimes[2] = t
        
        t = df['% Period'][df.idxmax()[3]]
        if t > df['% Period'][8]:
            t -= 1
        peakTimes[3] = t
        
        t = df['% Period'][df.idxmax()[4]]
        if t > df['% Period'][11]:
            t -= 1
        peakTimes[4] = t
        
        peakTimes[5] = df['% Period'][df.idxmax()[5]]
        peakTimes[6] = df['% Period'][df.idxmax()[6]]
        
        
    elif sourceType == 'rpush':
        
        peakvals[0] = df.min()[0]
        peakvals[1] = df.min()[1]
        peakvals[4] = df.min()[4]
        peakvals[5] = df.min()[5]
        peakvals[6] = df.min()[6]
        
        peakIdx[0] = df.idxmin()[0]
        peakIdx[1] = df.idxmin()[1]
        peakIdx[4] = df.idxmin()[4]
        peakIdx[5] = df.idxmin()[5]
        peakIdx[6] = df.idxmin()[6]
        
        t = df['% Period'][df.idxmin()[0]]
        if t > df['% Period'][9]:
            t -= 1
        peakTimes[0] = t
        
        t = df['% Period'][df.idxmin()[1]]
        if t > df['% Period'][8]:
            t -= 1
        peakTimes[1] = t
        
        peakTimes[4] = df['% Period'][df.idxmin()[4]]
        peakTimes[5] = df['% Period'][df.idxmin()[5]]
        peakTimes[6] = df['% Period'][df.idxmin()[6]]
        
        
    elif sourceType == 'rpull':
        
        peakvals[1] = df.min()[1]
        peakvals[2] = df.min()[2]
        
        peakIdx[1] = df.idxmin()[1]
        peakIdx[2] = df.idxmin()[2]

        for i in range(4,7):
            peakvals[i] = df.min()[i]
            peakIdx[i] = df.idxmin()[i]
            
        peakTimes[1] = df['% Period'][df.idxmin()[1]]
        peakTimes[2] = df['% Period'][df.idxmin()[2]]
        
        t = df['% Period'][df.idxmin()[4]]
        if t < df['% Period'][5]:
            t += 1
        peakTimes[4] = t
        
        t = df['% Period'][df.idxmin()[5]]
        if t < df['% Period'][6]:
            t += 1
        peakTimes[5] = t
        
        t = df['% Period'][df.idxmin()[6]]
        if t > df['% Period'][8]:
            t -= 1
        peakTimes[6] = t
        

    
    return peakvals, peakTimes, peakIdx