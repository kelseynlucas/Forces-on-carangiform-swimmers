# -*- coding: utf-8 -*-
"""
Created on Fri Aug 03 10:57:54 2018

@author: Kelsey N. Lucas


Script to extract fish kinematics, pressure, and forces, assemble to a long 
dataframe, and export for stats processing in R

Make sure to load InterX and the helper functions!
"""

import numpy as np
import pandas as pd


#Load some built-in file-reading functions
from os import listdir
from os.path import isfile, join

#set path for bluegill & trout data
bgpath = 'E:/fish pressure_v2/bluegill'
trpath = 'E:/fish pressure_v2/trout'


#import existing kinematics dataframes
bgKinem = pd.read_excel(bgpath + '/all_kinematics_data_summary_2.5BLs.xlsx',header=0,index_col=0,sheetname='Sheet1')
trKinem = pd.read_excel(trpath + '/all_kinematics_data_summary_2.5BLs.xlsx',header=0,index_col=0,sheetname='Sheet1')

#add species label to each dataframe
bgKinem['species'] = np.repeat('bluegill',len(bgKinem['period']))
trKinem['species'] = np.repeat('trout',len(trKinem['period']))

#assign lateral areas
#first, make a dictionary of the lateral area of each individual
latAreasm2 = {'klbg2': 0.00229634, 'klbg3': 0.002814104, 'klbg4': 0.002140139, 'klbg5': 0.002082938, 'klbg7': 0.001834378,'trout4': 0.001679883, 'trout7': 0.001750481, 'trout8': 0.001645867}
#assemble a list of the lateral area for each entry in dataframes
bgLatArea = [latAreasm2[bgKinem['individual'][i]] for i in range(len(bgKinem['individual']))]
trLatArea = [latAreasm2[trKinem['individual'][i]] for i in range(len(trKinem['individual']))]
#add lateral areas to the dataframes
bgKinem['latAreasm2'] = bgLatArea
trKinem['latAreasm2'] = trLatArea

#drop unneeded columns
bgdf = bgKinem.drop(['tbMaxIdx','tbMinIdx','maxCaudAngle','sourceNumFrames'],axis=1)
trdf = trKinem.drop(['tbMaxIdx','tbMinIdx','maxCaudAngle','sourceNumFrames'],axis=1)





#collect any other information that is 1 value per trial (whole body CF)

#Set up columns to stick data into

#CFx - whole body
bgdf['meanCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p1maxCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p1minCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p1maxTimeCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p1minTimeCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p2maxCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p2minCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p2maxTimeCFxWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['p2minTimeCFxWholeBody'] = np.zeros_like(bgdf['sequence'])

#CFy - whole body
bgdf['meanCFyWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFyWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['maxCFyWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['minCFyWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['maxTimeCFyWholeBody'] = np.zeros_like(bgdf['sequence'])
bgdf['minTimeCFyWholeBody'] = np.zeros_like(bgdf['sequence'])

#CFx - whole body
trdf['meanCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['netCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p1maxCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p1minCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p1maxTimeCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p1minTimeCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p2maxCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p2minCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p2maxTimeCFxWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['p2minTimeCFxWholeBody'] = np.zeros_like(trdf['sequence'])

#CFy - whole body
trdf['meanCFyWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['netCFyWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['maxCFyWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['minCFyWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['maxTimeCFyWholeBody'] = np.zeros_like(trdf['sequence'])
trdf['minTimeCFyWholeBody'] = np.zeros_like(trdf['sequence'])

#get data for each trial
for seq in bgdf['sequence'].unique():
    
    #get index of row to save into
    row = bgdf['sequence'].loc[bgdf['sequence']==seq].index[0]
    
    #get whole body CFx
    CFxWholeBody = get_NETforce_vs_time_wholebody_all_SIDESADDED(bgpath, seq, 'Fx')
    
    #get and save metrics of interest
    bgdf.loc[row,'meanCFxWholeBody'] = CFxWholeBody.mean()['net']
    bgdf.loc[row,'netCFxWholeBody'] = CFxWholeBody.sum()['net']
    bgdf.loc[row,'p1maxCFxWholeBody'] = CFxWholeBody['net'][0:8].max()
    bgdf.loc[row,'p1minCFxWholeBody'] = CFxWholeBody['net'][0:8].min()
    bgdf.loc[row,'p1maxTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmax(CFxWholeBody['net'][0:8])]
    bgdf.loc[row,'p1minTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmin(CFxWholeBody['net'][0:8])]
    bgdf.loc[row,'p2maxCFxWholeBody'] = CFxWholeBody['net'][8:15].max()
    bgdf.loc[row,'p2minCFxWholeBody'] = CFxWholeBody['net'][8:15].min()
    bgdf.loc[row,'p2maxTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmax(CFxWholeBody['net'][8:15])]
    bgdf.loc[row,'p2minTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmin(CFxWholeBody['net'][8:15])]

    #get whole body CFy
    CFyWholeBody = get_NETforce_vs_time_wholebody_all_SIDESADDED(bgpath, seq, 'Fy')
    
    #get and save metrics of interest
    bgdf.loc[row,'meanCFyWholeBody'] = CFyWholeBody.mean()['net']
    bgdf.loc[row,'netCFyWholeBody'] = CFyWholeBody.sum()['net']
    bgdf.loc[row,'maxCFyWholeBody'] = CFyWholeBody['net'].max()
    bgdf.loc[row,'minCFyWholeBody'] = CFyWholeBody['net'].min()
    bgdf.loc[row,'maxTimeCFyWholeBody'] = CFyWholeBody['% Period'][np.argmax(CFyWholeBody['net'])]
    bgdf.loc[row,'minTimeCFyWholeBody'] = CFyWholeBody['% Period'][np.argmin(CFyWholeBody['net'])]

for seq in trdf['sequence'].unique():
    
    #get index of row to save into
    row = trdf['sequence'].loc[trdf['sequence']==seq].index[0]
    
    #get whole body CFx
    CFxWholeBody = get_NETforce_vs_time_wholebody_all_SIDESADDED(trpath, seq, 'Fx')
    
    #get and save metrics of interest
    trdf.loc[row,'meanCFxWholeBody'] = CFxWholeBody.mean()['net']
    trdf.loc[row,'netCFxWholeBody'] = CFxWholeBody.sum()['net']
    trdf.loc[row,'p1maxCFxWholeBody'] = CFxWholeBody['net'][0:8].max()
    trdf.loc[row,'p1minCFxWholeBody'] = CFxWholeBody['net'][0:8].min()
    trdf.loc[row,'p1maxTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmax(CFxWholeBody['net'][0:8])]
    trdf.loc[row,'p1minTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmin(CFxWholeBody['net'][0:8])]
    trdf.loc[row,'p2maxCFxWholeBody'] = CFxWholeBody['net'][8:15].max()
    trdf.loc[row,'p2minCFxWholeBody'] = CFxWholeBody['net'][8:15].min()
    trdf.loc[row,'p2maxTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmax(CFxWholeBody['net'][8:15])]
    trdf.loc[row,'p2minTimeCFxWholeBody'] = CFxWholeBody['% Period'][np.argmin(CFxWholeBody['net'][8:15])]
    
    #get whole body CFy
    CFyWholeBody = get_NETforce_vs_time_wholebody_all_SIDESADDED(trpath, seq, 'Fy')
    
    #get and save metrics of interest
    trdf.loc[row,'meanCFyWholeBody'] = CFyWholeBody.mean()['net']
    trdf.loc[row,'netCFyWholeBody'] = CFyWholeBody.sum()['net']
    trdf.loc[row,'maxCFyWholeBody'] = CFyWholeBody['net'].max()
    trdf.loc[row,'minCFyWholeBody'] = CFyWholeBody['net'].min()
    trdf.loc[row,'maxTimeCFyWholeBody'] = CFyWholeBody['% Period'][np.argmax(CFyWholeBody['net'])]
    trdf.loc[row,'minTimeCFyWholeBody'] = CFyWholeBody['% Period'][np.argmin(CFyWholeBody['net'])]
    
    

#Will need to repeat dataframe rows multiple times to deal with entries that 
#share previous parameters (segment data, etc)

#repeat dataframes 56x (once for each segment, side, forcetype, and pressuretype combo)
bgdf = pd.concat([bgdf]*56)
bgdf = bgdf.sort_index()
bgdf = bgdf.reset_index()
trdf = pd.concat([trdf]*56)
trdf = trdf.sort_index()
trdf = trdf.reset_index()


segmentNames = [1,2,3,4,5,6,7]
#repeat segmNames #forcetype*#sides*#ptypes=8x for each individ
bgSegments = np.tile(np.repeat(segmentNames,8),15)
trSegments = np.tile(np.repeat(segmentNames,8),9)

#for each segment, repeat side designation #forcetype*#ptypes=4x
sides = np.concatenate((np.tile(['left'],4),np.tile(['right'],4)),axis=0)
#repeat these labels #individuals*#segm times
bgSides = np.tile(sides,15*7)
trSides = np.tile(sides,9*7)

#for each side, will need forceType #ptype times=2x
forceTypes=['thrust','thrust','drag','drag']
#repeat these labels #individuals*#segm*#sides=15*7*2x or 9*7*2
bgForceTypes = np.tile(forceTypes,15*7*2)
trForceTypes = np.tile(forceTypes,9*7*2)

#for each force time, will need pType 1x
pTypes = ['lo','hi']
#repeat these labels #individ*#segm*#sides*#forceTypes=15*7*2*2
bgPTypes = np.tile(pTypes,15*7*2*2)
trPTypes = np.tile(pTypes,9*7*2*2)


#assemble dataframes
bgdf['segment'] = bgSegments
bgdf['side'] = bgSides
bgdf['forceType'] = bgForceTypes
bgdf['pressType'] = bgPTypes
trdf['segment'] = trSegments
trdf['side'] = trSides
trdf['forceType'] = trForceTypes
trdf['pressType'] = trPTypes



#Make columns for data per segment

#amplitude
bgdf['ampInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['ampInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['ampMaxTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['ampMaxTimeInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['ampMinTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['ampMinTimeInSegm'] = np.zeros_like(trdf['sequence'])

#body angle
bgdf['bodyAngleInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['bodyAngleInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['bodyAngleMaxTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['bodyAngleMaxTimeInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['bodyAngleMinTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['bodyAngleMinTimeInSegm'] = np.zeros_like(trdf['sequence'])

#angle of attack
bgdf['AoAInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['AoAInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['AoAMaxTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['AoAMaxTimeInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['AoAMinTimeInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['AoAMinTimeInSegm'] = np.zeros_like(trdf['sequence'])

#Cp 
bgdf['meanCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['absMeanCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['maxCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['minCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['maxTimeCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['minTimeCpsInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['meanCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['absMeanCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['maxCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['minCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['maxTimeCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['minTimeCpsInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])

trdf['meanCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['absMeanCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['maxCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['minCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['maxTimeCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['minTimeCpsInSegm'] = np.zeros_like(trdf['sequence'])
trdf['meanCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['absMeanCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['maxCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['minCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['maxTimeCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['minTimeCpsInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])

#CFx (total)
bgdf['meanCFxInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFxInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['meanCFxInSegm'] = np.zeros_like(trdf['sequence'])
trdf['netCFxInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['meanCFxInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFxInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
trdf['meanCFxInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['netCFxInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])

#CFy (total)
bgdf['meanCFyInSegm'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFyInSegm'] = np.zeros_like(bgdf['sequence'])
trdf['meanCFyInSegm'] = np.zeros_like(trdf['sequence'])
trdf['netCFyInSegm'] = np.zeros_like(trdf['sequence'])
bgdf['meanCFyInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFyInSegmSidesAvg'] = np.zeros_like(bgdf['sequence'])
trdf['meanCFyInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['netCFyInSegmSidesAvg'] = np.zeros_like(trdf['sequence'])

#CFx by segment
bgdf['meanCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['peakCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['absMeanCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['absNetCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['absPeakCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['peakTimeCFxByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['areaUnderCurveByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['ampAtPeakByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['angleAtPeakByMech'] = np.zeros_like(bgdf['sequence'])
bgdf['meanCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['netCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['peakCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['absMeanCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['absNetCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['absPeakCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['peakTimeCFxByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['areaUnderCurveByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['ampAtPeakByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])
bgdf['angleAtPeakByMechSidesAvg'] = np.zeros_like(bgdf['sequence'])

trdf['meanCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['netCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['peakCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['absMeanCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['absNetCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['absPeakCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['peakTimeCFxByMech'] = np.zeros_like(trdf['sequence'])
trdf['areaUnderCurveByMech'] = np.zeros_like(trdf['sequence'])
trdf['ampAtPeakByMech'] = np.zeros_like(trdf['sequence'])
trdf['angleAtPeakByMech'] = np.zeros_like(trdf['sequence'])
trdf['meanCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['netCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['peakCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['absMeanCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['absNetCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['absPeakCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['peakTimeCFxByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['areaUnderCurveByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['ampAtPeakByMechSidesAvg'] = np.zeros_like(trdf['sequence'])
trdf['angleAtPeakByMechSidesAvg'] = np.zeros_like(trdf['sequence'])

#get all data that is per segment (bluegill)
for seq in bgdf['sequence'].unique():
    print(seq)
    
    placementIdx = [i for i in range(len(bgdf['sequence'])) if bgdf['sequence'][i] == seq]    
    
    #get segment amplitude 
    amp = get_amp_vs_time_binned(bgpath,seq)
    ampInSegm = (abs(amp.min())+amp.max())/2.
    ampInSegm = ampInSegm[0:7]
    ampMaxTimes = amp['% Period'][amp.idxmax()[0:7]].values
    ampMinTimes = amp['% Period'][amp.idxmin()[0:7]].values
    #make list of amps in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    bgdf.loc[placementIdx,'ampInSegm'] = np.repeat(ampInSegm.values,8)
    bgdf.loc[placementIdx,'ampMaxTimeInSegm'] = np.repeat(ampMaxTimes,8)
    bgdf.loc[placementIdx,'ampMinTimeInSegm'] = np.repeat(ampMinTimes,8)
    
    #get segment angle
    angle = get_angles_vs_time_binned(bgpath,seq)
    bodyAngleInSegm = (abs(angle.min())+angle.max())/2.
    bodyAngleInSegm = bodyAngleInSegm[0:7]
    angleMaxTimes = angle['% Period'][angle.idxmax()[0:7]].values
    angleMinTimes = angle['% Period'][angle.idxmin()[0:7]].values
    #make list of angles in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    bgdf.loc[placementIdx,'bodyAngleInSegm'] = np.repeat(bodyAngleInSegm.values,8)
    bgdf.loc[placementIdx,'bodyAngleMaxTimeInSegm'] = np.repeat(angleMaxTimes,8)
    bgdf.loc[placementIdx,'bodyAngleMinTimeInSegm'] = np.repeat(angleMinTimes,8)    
    
    
    #get angle of attack
    AoA = get_AoA_vs_time_binned(bgpath,seq)
    AoAInSegm = (abs(AoA.min())+AoA.max())/2.
    AoAInSegm = AoAInSegm[0:7]
    AoAMaxTimes = amp['% Period'][AoA.idxmax()[0:7]].values
    AoAMinTimes = amp['% Period'][AoA.idxmin()[0:7]].values
    #make list of angles in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    bgdf.loc[placementIdx,'AoAInSegm'] = np.repeat(AoAInSegm.values,8)
    bgdf.loc[placementIdx,'AoAMaxTimeInSegm'] = np.repeat(AoAMaxTimes,8)
    bgdf.loc[placementIdx,'AoAMinTimeInSegm'] = np.repeat(AoAMinTimes,8)
    
    
    #Get pressure information
    
    CpsLeft = get_MEANpressure_vs_time_binned_all_TOP(bgpath,seq)
    CpsRight = get_MEANpressure_vs_time_binned_all_BOTTOM(bgpath,seq)
    
    meanCpsLeft = CpsLeft.mean()[0:7]
    meanCpsRight = CpsRight.mean()[0:7]
    
    meanMeanCps = np.mean([meanCpsLeft.values, meanCpsRight.values], axis=0)
    
    
    maxCpsLeft, maxTimeCpsLeft = get_MEANCps_peaks(CpsLeft, 'max')
    maxCpsRight, maxTimeCpsRight = get_MEANCps_peaks(CpsRight, 'max')
    minCpsLeft, minTimeCpsLeft = get_MEANCps_peaks(CpsLeft, 'min')
    minCpsRight, minTimeCpsRight = get_MEANCps_peaks(CpsRight, 'min')
    
    meanMaxCps = np.mean([maxCpsLeft, maxCpsRight], axis=0)
    meanMaxTimeCps = np.mean([maxTimeCpsLeft, maxTimeCpsRight], axis=0)
    meanMinCps = np.mean([minCpsLeft, minCpsRight], axis=0)
    meanMinTimeCps = np.mean([minTimeCpsLeft, minTimeCpsRight], axis=0)

    #save pressure data to df
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left')]
    #save left side data
    bgdf.loc[placementIdx,'meanCpsInSegm'] = np.repeat(meanCpsLeft.values,4)
    bgdf.loc[placementIdx,'absMeanCpsInSegm'] = np.repeat(abs(meanCpsLeft.values),4)
    bgdf.loc[placementIdx,'maxCpsInSegm'] = np.repeat(maxCpsLeft,4)
    bgdf.loc[placementIdx,'minCpsInSegm'] = np.repeat(minCpsLeft,4)
    bgdf.loc[placementIdx,'maxTimeCpsInSegm'] = np.repeat(maxTimeCpsLeft,4)
    bgdf.loc[placementIdx,'minTimeCpsInSegm'] = np.repeat(minTimeCpsLeft,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right')]
    #save right side data
    bgdf.loc[placementIdx,'meanCpsInSegm'] = np.repeat(meanCpsRight.values,4)
    bgdf.loc[placementIdx,'absMeanCpsInSegm'] = np.repeat(abs(meanCpsRight.values),4)
    bgdf.loc[placementIdx,'maxCpsInSegm'] = np.repeat(maxCpsRight,4)
    bgdf.loc[placementIdx,'minCpsInSegm'] = np.repeat(minCpsRight,4)
    bgdf.loc[placementIdx,'maxTimeCpsInSegm'] = np.repeat(maxTimeCpsRight,4)
    bgdf.loc[placementIdx,'minTimeCpsInSegm'] = np.repeat(minTimeCpsRight,4)
    
    #make placementIdx for side averages
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq)]
    #save mean side data
    bgdf.loc[placementIdx,'meanCpsInSegmSidesAvg'] = np.repeat(meanMeanCps,8)
    bgdf.loc[placementIdx,'absMeanCpsInSegmSidesAvg'] = np.repeat(abs(meanMeanCps),8)
    bgdf.loc[placementIdx,'maxCpsInSegmSidesAvg'] = np.repeat(meanMaxCps,8)
    bgdf.loc[placementIdx,'minCpsInSegmSidesAvg'] = np.repeat(meanMinCps,8)
    bgdf.loc[placementIdx,'maxTimeCpsInSegmSidesAvg'] = np.repeat(meanMaxTimeCps,8)
    bgdf.loc[placementIdx,'minTimeCpsInSegmSidesAvg'] = np.repeat(meanMinTimeCps,8)
    
    
    
    #Get total CFx information
    
    CFxLeft = get_NETforce_vs_time_binned_all_TOP(bgpath, seq, 'Fx')
    CFxRight = get_NETforce_vs_time_binned_all_BOTTOM(bgpath, seq, 'Fx')
    
    meanCFxLeft = CFxLeft.mean()[0:7]
    meanCFxRight = CFxRight.mean()[0:7]
    
    netCFxLeft = CFxLeft.sum()[0:7]
    netCFxRight = CFxRight.sum()[0:7]
    
    meanMeanCFx = np.mean([meanCFxLeft, meanCFxRight], axis=0)
    meanNetCFx = np.mean([netCFxLeft, netCFxRight], axis=0)
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left')]
    #save left side data
    bgdf.loc[placementIdx,'meanCFxInSegm'] = np.repeat(meanCFxLeft.values,4)
    bgdf.loc[placementIdx,'netCFxInSegm'] = np.repeat(netCFxLeft.values,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right')]
    #save right side data
    bgdf.loc[placementIdx,'meanCFxInSegm'] = np.repeat(meanCFxRight.values,4)
    bgdf.loc[placementIdx,'netCFxInSegm'] = np.repeat(netCFxRight.values,4)
    
    #make placementIdx for average of sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq)]
    #save mean side data
    bgdf.loc[placementIdx,'meanCFxInSegmSidesAvg'] = np.repeat(meanMeanCFx,8)
    bgdf.loc[placementIdx,'netCFxInSegmSidesAvg'] = np.repeat(meanNetCFx,8)
    
    
    
    #Get total CFy information
    
    CFyLeft = get_NETforce_vs_time_binned_all_TOP(bgpath, seq, 'Fy')
    CFyRight = get_NETforce_vs_time_binned_all_BOTTOM(bgpath, seq, 'Fy')
    
    meanCFyLeft = CFyLeft.mean()[0:7]
    meanCFyRight = CFyRight.mean()[0:7]
    
    netCFyLeft = CFyLeft.sum()[0:7]
    netCFyRight = CFyRight.sum()[0:7]
    
    meanMeanCFy = np.mean([meanCFyLeft, meanCFyRight], axis=0)
    meanNetCFy = np.mean([netCFyLeft, netCFyRight], axis=0)
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left')]
    #save left side data
    bgdf.loc[placementIdx,'meanCFyInSegm'] = np.repeat(meanCFyLeft.values,4)
    bgdf.loc[placementIdx,'netCFyInSegm'] = np.repeat(netCFyLeft.values,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right')]
    #save right side data
    bgdf.loc[placementIdx,'meanCFyInSegm'] = np.repeat(meanCFyRight.values,4)
    bgdf.loc[placementIdx,'netCFyInSegm'] = np.repeat(netCFyRight.values,4)
    
    #make placementIdx for average of sides
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq)]
    #save mean side data
    bgdf.loc[placementIdx,'meanCFyInSegmSidesAvg'] = np.repeat(meanMeanCFy,8)
    bgdf.loc[placementIdx,'netCFyInSegmSidesAvg'] = np.repeat(meanNetCFy,8)

    
    
    #Get single force type info - Fx
    
    thrustHiLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(bgpath, seq, 'fpush')
    thrustHiRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(bgpath, seq, 'fpush') 
    
    thrustLoLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(bgpath, seq, 'fpull')
    thrustLoRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(bgpath, seq, 'fpull')
    
    dragHiLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(bgpath, seq, 'rpush')
    dragHiRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(bgpath, seq, 'rpush')
    
    dragLoLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(bgpath, seq, 'rpull')
    dragLoRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(bgpath, seq, 'rpull')
    
    #calculate values of interest
    
    #means
    meanThrustHiLeft = thrustHiLeft.mean()[0:7]
    meanThrustHiRight = thrustHiRight.mean()[0:7]
    
    meanThrustLoLeft = thrustLoLeft.mean()[0:7]
    meanThrustLoRight = thrustLoRight.mean()[0:7]
    
    meanDragHiLeft = dragHiLeft.mean()[0:7]
    meanDragHiRight = dragHiRight.mean()[0:7]
    
    meanDragLoLeft = dragLoLeft.mean()[0:7]
    meanDragLoRight = dragLoRight.mean()[0:7]
    
    meanThrustHi = np.mean([meanThrustHiLeft, meanThrustHiRight], axis=0)
    meanThrustLo = np.mean([meanThrustLoLeft, meanThrustLoRight], axis=0)
    meanDragHi = np.mean([meanDragHiLeft, meanDragHiRight], axis=0)
    meanDragLo = np.mean([meanDragLoLeft, meanDragLoRight], axis=0)
    
    #net
    netThrustHiLeft = thrustHiLeft.sum()[0:7]
    netThrustHiRight = thrustHiRight.sum()[0:7]
    
    netThrustLoLeft = thrustLoLeft.sum()[0:7]
    netThrustLoRight = thrustLoRight.sum()[0:7]
    
    netDragHiLeft = dragHiLeft.sum()[0:7]
    netDragHiRight = dragHiRight.sum()[0:7]
    
    netDragLoLeft = dragLoLeft.sum()[0:7]
    netDragLoRight = dragLoRight.sum()[0:7]
    
    netThrustHi = np.mean([netThrustHiLeft, netThrustHiRight], axis=0)
    netThrustLo = np.mean([netThrustLoLeft, netThrustLoRight], axis=0)
    netDragHi = np.mean([netDragHiLeft, netDragHiRight], axis=0)
    netDragLo = np.mean([netDragLoLeft, netDragLoRight], axis=0)
    
    #peak values
    peakThrustHiLeft, peakTimeThrustHiLeft, peakIdxThrustHiLeft = get_NETCFx_source_peaks(thrustHiLeft,'fpush')
    peakThrustHiRight, peakTimeThrustHiRight, peakIdxThrustHiRight = get_NETCFx_source_peaks(thrustHiRight,'fpush')
    
    peakThrustLoLeft, peakTimeThrustLoLeft, peakIdxThrustLoLeft = get_NETCFx_source_peaks(thrustLoLeft,'fpull')
    peakThrustLoRight, peakTimeThrustLoRight, peakIdxThrustLoRight = get_NETCFx_source_peaks(thrustLoRight,'fpull')
    
    peakDragHiLeft, peakTimeDragHiLeft, peakIdxDragHiLeft = get_NETCFx_source_peaks(dragHiLeft,'rpush')
    peakDragHiRight, peakTimeDragHiRight, peakIdxDragHiRight = get_NETCFx_source_peaks(dragHiRight,'rpush')
    
    peakDragLoLeft, peakTimeDragLoLeft, peakIdxDragLoLeft = get_NETCFx_source_peaks(dragLoLeft,'rpull')
    peakDragLoRight, peakTimeDragLoRight, peakIdxDragLoRight = get_NETCFx_source_peaks(dragLoRight,'rpull')
    
    meanPeakThrustHi = np.mean([peakThrustHiLeft, peakThrustHiRight], axis=0)
    meanPeakThrustLo = np.mean([peakThrustLoLeft, peakThrustLoRight], axis=0)
    meanPeakDragHi = np.mean([peakDragHiLeft, peakDragHiRight], axis=0)
    meanPeakDragLo = np.mean([peakDragLoLeft, peakDragLoRight], axis=0)
    
    meanPeakTimeThrustHi = np.mean([peakTimeThrustHiLeft, peakTimeThrustHiRight], axis=0)
    meanPeakTimeThrustLo = np.mean([peakTimeThrustLoLeft, peakTimeThrustLoRight], axis=0)
    meanPeakTimeDragHi = np.mean([peakTimeDragHiLeft, peakTimeDragHiRight], axis=0)
    meanPeakTimeDragLo = np.mean([peakTimeDragLoLeft, peakTimeDragLoRight], axis=0)
    
    #amp and angle at peak values
    ampatPeakThrustHiLeft = np.zeros((7,))
    ampatPeakThrustHiRight = np.zeros((7,))
    ampatPeakThrustLoLeft = np.zeros((7,))
    ampatPeakThrustLoRight = np.zeros((7,))
    ampatPeakDragHiLeft = np.zeros((7,))
    ampatPeakDragHiRight = np.zeros((7,))
    ampatPeakDragLoLeft = np.zeros((7,))
    ampatPeakDragLoRight = np.zeros((7,))
    
    ampatPeakThrustHi = np.zeros((7,))
    ampatPeakThrustLo = np.zeros((7,))
    ampatPeakDragHi = np.zeros((7,))
    ampatPeakDragLo = np.zeros((7,))    
    
    angleatPeakThrustHiLeft = np.zeros((7,))
    angleatPeakThrustHiRight = np.zeros((7,))
    angleatPeakThrustLoLeft = np.zeros((7,))
    angleatPeakThrustLoRight = np.zeros((7,))
    angleatPeakDragHiLeft = np.zeros((7,))
    angleatPeakDragHiRight = np.zeros((7,))
    angleatPeakDragLoLeft = np.zeros((7,))
    angleatPeakDragLoRight = np.zeros((7,))
    
    angleatPeakThrustHi = np.zeros((7,))
    angleatPeakThrustLo = np.zeros((7,))
    angleatPeakDragHi = np.zeros((7,))
    angleatPeakDragLo = np.zeros((7,))

    
    
    
    for i in range(7):
        
        try:
            ampatPeakThrustHiLeft[i] = amp.iloc[int(peakIdxThrustHiLeft[i]), i]
        except:
            ampatPeakThrustHiLeft[i] = np.nan
        
        try:
            ampatPeakThrustHiRight[i] = amp.iloc[int(peakIdxThrustHiRight[i]), i]
        except:
            ampatPeakThrustHiRight[i] = np.nan
        
        ampatPeakThrustHi[i] = np.mean([abs(ampatPeakThrustHiLeft[i]),abs(ampatPeakThrustHiRight[i])])
        
        try:
            ampatPeakThrustLoLeft[i] = amp.iloc[int(peakIdxThrustLoLeft[i]), i]
        except:
            ampatPeakThrustLoLeft[i] = np.nan
        
        try:
            ampatPeakThrustLoRight[i] = amp.iloc[int(peakIdxThrustLoRight[i]), i]
        except:
            ampatPeakThrustLoRight[i] = np.nan
            
        ampatPeakThrustLo[i] = np.mean([abs(ampatPeakThrustLoLeft[i]),abs(ampatPeakThrustLoRight[i])])
            
        try:
            ampatPeakDragHiLeft[i] = amp.iloc[int(peakIdxDragHiLeft[i]), i]
        except:
            ampatPeakDragHiLeft[i] = np.nan
            
        try:
            ampatPeakDragHiRight[i] = amp.iloc[int(peakIdxDragHiRight[i]), i]
        except:
            ampatPeakDragHiRight[i] = np.nan
            
        ampatPeakDragHi[i] = np.mean([abs(ampatPeakDragHiLeft[i]),abs(ampatPeakDragHiRight[i])])
        
        try:
            ampatPeakDragLoLeft[i] = amp.iloc[int(peakIdxDragLoLeft[i]), i]
        except:
            ampatPeakDragLoLeft[i] = np.nan
          
        try:
            ampatPeakDragLoRight[i] = amp.iloc[int(peakIdxDragLoRight[i]), i]
        except:
            ampatPeakDragLoRight[i] = np.nan
        
        ampatPeakDragLo[i] = np.mean([abs(ampatPeakDragLoLeft[i]),abs(ampatPeakDragLoRight[i])])
    
    
    for i in range(7):
        try:
            angleatPeakThrustHiLeft[i] = angle.iloc[int(peakIdxThrustHiLeft[i]), i]
        except:
            angleatPeakThrustHiLeft[i] = np.nan
        
        try:
            angleatPeakThrustHiRight[i] = angle.iloc[int(peakIdxThrustHiRight[i]), i]
        except:
            angleatPeakThrustHiRight[i] = np.nan
            
        angleatPeakThrustHi[i] = np.mean([abs(angleatPeakThrustHiLeft[i]),abs(angleatPeakThrustHiRight[i])])
        
        try:
            angleatPeakThrustLoLeft[i] = angle.iloc[int(peakIdxThrustLoLeft[i]), i]
        except:
            angleatPeakThrustLoLeft[i] = np.nan
        
        try:
            angleatPeakThrustLoRight[i] = angle.iloc[int(peakIdxThrustLoRight[i]), i]
        except:
            angleatPeakThrustLoRight[i] = np.nan
            
        angleatPeakThrustLo[i] = np.mean([abs(angleatPeakThrustLoLeft[i]),abs(angleatPeakThrustLoRight[i])])
        
        try:
            angleatPeakDragHiLeft[i] = angle.iloc[int(peakIdxDragHiLeft[i]), i]
        except:
            angleatPeakDragHiLeft[i] = np.nan
            
        try:
            angleatPeakDragHiRight[i] = angle.iloc[int(peakIdxDragHiRight[i]), i]
        except:
            angleatPeakDragHiRight[i] = np.nan
        
        angleatPeakDragHi[i] = np.mean([abs(angleatPeakDragHiLeft[i]),abs(angleatPeakDragHiRight[i])])        
        
        try:
            angleatPeakDragLoLeft[i] = angle.iloc[int(peakIdxDragLoLeft[i]), i]
        except:
            angleatPeakDragLoLeft[i] = np.nan
            
        try:
            angleatPeakDragLoRight[i] = angle.iloc[int(peakIdxDragLoRight[i]), i]
        except:
            angleatPeakDragLoRight[i] = np.nan
            
        angleatPeakDragLo[i] = np.mean([abs(angleatPeakDragLoLeft[i]),abs(angleatPeakDragLoRight[i])])

    
    #area under curve
    areaThrustHiLeft = []
    areaThrustHiRight = []
    areaThrustLoLeft = []
    areaThrustLoRight = []
    areaDragHiLeft = []
    areaDragHiRight = []
    areaDragLoLeft = []
    areaDragLoRight = []
    
    meanAreaThrustHi = []
    meanAreaThrustLo = []
    meanAreaDragHi = []
    meanAreaDragLo = []
    
    for col in range(7):
        
        areaThrustHiLeft.append(np.trapz(thrustHiLeft[col],thrustHiLeft['% Period']))
        areaThrustHiRight.append(np.trapz(thrustHiRight[col],thrustHiRight['% Period']))
        meanAreaThrustHi.append(np.mean([areaThrustHiLeft[-1], areaThrustHiRight[-1]]))
        
        areaThrustLoLeft.append(np.trapz(thrustLoLeft[col],thrustLoLeft['% Period']))
        areaThrustLoRight.append(np.trapz(thrustLoRight[col],thrustLoRight['% Period']))
        meanAreaThrustLo.append(np.mean([areaThrustLoLeft[-1], areaThrustLoRight[-1]]))
        
        areaDragHiLeft.append(np.trapz(dragHiLeft[col],dragHiLeft['% Period']))
        areaDragHiRight.append(np.trapz(dragHiRight[col],dragHiRight['% Period']))
        meanAreaDragHi.append(np.mean([areaDragHiLeft[-1], areaDragHiRight[-1]]))
        
        areaDragLoLeft.append(np.trapz(dragLoLeft[col],dragLoLeft['% Period']))
        areaDragLoRight.append(np.trapz(dragLoRight[col],dragLoRight['% Period']))
        meanAreaDragLo.append(np.mean([areaDragLoLeft[-1], areaDragLoRight[-1]]))
    
    
    #save all values
    
    #thrust hi
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left') & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanThrustHiLeft.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netThrustHiLeft.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakThrustHiLeft
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustHiLeft.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustHiLeft.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustHiLeft)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustHiLeft
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustHiLeft
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustHiLeft
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustHiLeft
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right') & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanThrustHiRight.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netThrustHiRight.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakThrustHiRight
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustHiRight.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustHiRight.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustHiRight)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustHiRight
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustHiRight
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustHiRight
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustHiRight
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanThrustHi,2)
    bgdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netThrustHi,2)
    bgdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakThrustHi,2)
    bgdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanThrustHi),2)
    bgdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netThrustHi),2)
    bgdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakThrustHi),2)
    bgdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeThrustHi,2)
    bgdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaThrustHi,2)
    bgdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakThrustHi,2)
    bgdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakThrustHi,2)



    #thrust lo
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left') & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanThrustLoLeft.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netThrustLoLeft.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakThrustLoLeft
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustLoLeft.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustLoLeft.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustLoLeft)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustLoLeft
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustLoLeft
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustLoLeft
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustLoLeft
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right') & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanThrustLoRight.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netThrustLoRight.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakThrustLoRight
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustLoRight.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustLoRight.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustLoRight)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustLoRight
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustLoRight
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustLoRight
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustLoRight
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['forceType'][i] == 'thrust') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanThrustLo,2)
    bgdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netThrustLo,2)
    bgdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakThrustLo,2)
    bgdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanThrustLo),2)
    bgdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netThrustLo),2)
    bgdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakThrustLo),2)
    bgdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeThrustLo,2)
    bgdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaThrustLo,2)
    bgdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakThrustLo,2)
    bgdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakThrustLo,2)
    
    #drag hi
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left') & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanDragHiLeft.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netDragHiLeft.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakDragHiLeft
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragHiLeft.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragHiLeft.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragHiLeft)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragHiLeft
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragHiLeft
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragHiLeft
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragHiLeft
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right') & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanDragHiRight.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netDragHiRight.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakDragHiRight
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragHiRight.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragHiRight.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragHiRight)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragHiRight
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragHiRight
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragHiRight
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragHiRight
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'hi')]
    bgdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanDragHi,2)
    bgdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netDragHi,2)
    bgdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakDragHi,2)
    bgdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanDragHi),2)
    bgdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netDragHi),2)
    bgdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakDragHi),2)
    bgdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeDragHi,2)
    bgdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaDragHi,2)
    bgdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakDragHi,2)
    bgdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakDragHi,2)
    
    
    #drag lo
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'left') & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanDragLoLeft.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netDragLoLeft.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakDragLoLeft
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragLoLeft.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragLoLeft.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragLoLeft)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragLoLeft
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragLoLeft
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragLoLeft
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragLoLeft
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['side'][i] == 'right') & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMech'] = meanDragLoRight.values
    bgdf.loc[placementIdx,'netCFxByMech'] = netDragLoRight.values
    bgdf.loc[placementIdx,'peakCFxByMech'] = peakDragLoRight
    bgdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragLoRight.values)
    bgdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragLoRight.values)
    bgdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragLoRight)
    bgdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragLoRight
    bgdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragLoRight
    bgdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragLoRight
    bgdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragLoRight
    placementIdx = [i for i in range(len(bgdf['sequence'])) if (bgdf['sequence'][i] == seq) & (bgdf['forceType'][i] == 'drag') & (bgdf['pressType'][i] == 'lo')]
    bgdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanDragLo,2)
    bgdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netDragLo,2)
    bgdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakDragLo,2)
    bgdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanDragLo),2)
    bgdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netDragLo),2)
    bgdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakDragLo),2)
    bgdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeDragLo,2)
    bgdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaDragLo,2)
    bgdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakDragLo,2)
    bgdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakDragLo,2)




#get all data that is per segment (trout)
for seq in trdf['sequence'].unique():
    print(seq)
    
    placementIdx = [i for i in range(len(trdf['sequence'])) if trdf['sequence'][i] == seq]    
    
    #get segment amplitude 
    amp = get_amp_vs_time_binned(trpath,seq)
    ampInSegm = (abs(amp.min())+amp.max())/2.
    ampInSegm = ampInSegm[0:7]
    ampMaxTimes = amp['% Period'][amp.idxmax()[0:7]].values
    ampMinTimes = amp['% Period'][amp.idxmin()[0:7]].values
    #make list of amps in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    trdf.loc[placementIdx,'ampInSegm'] = np.repeat(ampInSegm.values,8)
    trdf.loc[placementIdx,'ampMaxTimeInSegm'] = np.repeat(ampMaxTimes,8)
    trdf.loc[placementIdx,'ampMinTimeInSegm'] = np.repeat(ampMinTimes,8)
    
    #get segment angle
    angle = get_angles_vs_time_binned(trpath,seq)
    bodyAngleInSegm = (abs(angle.min())+angle.max())/2.
    bodyAngleInSegm = bodyAngleInSegm[0:7]
    angleMaxTimes = angle['% Period'][angle.idxmax()[0:7]].values
    angleMinTimes = angle['% Period'][angle.idxmin()[0:7]].values
    #make list of angles in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    trdf.loc[placementIdx,'bodyAngleInSegm'] = np.repeat(bodyAngleInSegm.values,8)
    trdf.loc[placementIdx,'bodyAngleMaxTimeInSegm'] = np.repeat(angleMaxTimes,8)
    trdf.loc[placementIdx,'bodyAngleMinTimeInSegm'] = np.repeat(angleMinTimes,8)    
    
    
    #get angle of attack
    AoA = get_AoA_vs_time_binned(trpath,seq)
    AoAInSegm = (abs(AoA.min())+AoA.max())/2.
    AoAInSegm = AoAInSegm[0:7]
    AoAMaxTimes = amp['% Period'][AoA.idxmax()[0:7]].values
    AoAMinTimes = amp['% Period'][AoA.idxmin()[0:7]].values
    #make list of angles in dimensions to match df.  Need each value to repeat #sides*#forceType*#PType = 2*2*2 = 8x
    trdf.loc[placementIdx,'AoAInSegm'] = np.repeat(AoAInSegm.values,8)
    trdf.loc[placementIdx,'AoAMaxTimeInSegm'] = np.repeat(AoAMaxTimes,8)
    trdf.loc[placementIdx,'AoAMinTimeInSegm'] = np.repeat(AoAMinTimes,8)
    
    
    #Get pressure information
    
    CpsLeft = get_MEANpressure_vs_time_binned_all_TOP(trpath,seq)
    CpsRight = get_MEANpressure_vs_time_binned_all_BOTTOM(trpath,seq)
    
    meanCpsLeft = CpsLeft.mean()[0:7]
    meanCpsRight = CpsRight.mean()[0:7]
    
    meanMeanCps = np.mean([meanCpsLeft.values, meanCpsRight.values], axis=0)
    
    
    maxCpsLeft, maxTimeCpsLeft = get_MEANCps_peaks(CpsLeft, 'max')
    maxCpsRight, maxTimeCpsRight = get_MEANCps_peaks(CpsRight, 'max')
    minCpsLeft, minTimeCpsLeft = get_MEANCps_peaks(CpsLeft, 'min')
    minCpsRight, minTimeCpsRight = get_MEANCps_peaks(CpsRight, 'min')
    
    meanMaxCps = np.mean([maxCpsLeft, maxCpsRight], axis=0)
    meanMaxTimeCps = np.mean([maxTimeCpsLeft, maxTimeCpsRight], axis=0)
    meanMinCps = np.mean([minCpsLeft, minCpsRight], axis=0)
    meanMinTimeCps = np.mean([minTimeCpsLeft, minTimeCpsRight], axis=0)
    
    #save pressure data to df
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left')]
    #save left side data
    trdf.loc[placementIdx,'meanCpsInSegm'] = np.repeat(meanCpsLeft.values,4)
    trdf.loc[placementIdx,'absMeanCpsInSegm'] = np.repeat(abs(meanCpsLeft.values),4)
    trdf.loc[placementIdx,'maxCpsInSegm'] = np.repeat(maxCpsLeft,4)
    trdf.loc[placementIdx,'minCpsInSegm'] = np.repeat(minCpsLeft,4)
    trdf.loc[placementIdx,'maxTimeCpsInSegm'] = np.repeat(maxTimeCpsLeft,4)
    trdf.loc[placementIdx,'minTimeCpsInSegm'] = np.repeat(minTimeCpsLeft,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right')]
    #save right side data
    trdf.loc[placementIdx,'meanCpsInSegm'] = np.repeat(meanCpsRight.values,4)
    trdf.loc[placementIdx,'absMeanCpsInSegm'] = np.repeat(abs(meanCpsRight.values),4)
    trdf.loc[placementIdx,'maxCpsInSegm'] = np.repeat(maxCpsRight,4)
    trdf.loc[placementIdx,'minCpsInSegm'] = np.repeat(minCpsRight,4)
    trdf.loc[placementIdx,'maxTimeCpsInSegm'] = np.repeat(maxTimeCpsRight,4)
    trdf.loc[placementIdx,'minTimeCpsInSegm'] = np.repeat(minTimeCpsRight,4)
    
    #make placementIdx for side averages
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq)]
    #save mean side data
    trdf.loc[placementIdx,'meanCpsInSegmSidesAvg'] = np.repeat(meanMeanCps,8)
    trdf.loc[placementIdx,'absMeanCpsInSegmSidesAvg'] = np.repeat(abs(meanMeanCps),8)
    trdf.loc[placementIdx,'maxCpsInSegmSidesAvg'] = np.repeat(meanMaxCps,8)
    trdf.loc[placementIdx,'minCpsInSegmSidesAvg'] = np.repeat(meanMinCps,8)
    trdf.loc[placementIdx,'maxTimeCpsInSegmSidesAvg'] = np.repeat(meanMaxTimeCps,8)
    trdf.loc[placementIdx,'minTimeCpsInSegmSidesAvg'] = np.repeat(meanMinTimeCps,8)
    
    
    #Get total CFx information
    
    CFxLeft = get_NETforce_vs_time_binned_all_TOP(trpath, seq, 'Fx')
    CFxRight = get_NETforce_vs_time_binned_all_BOTTOM(trpath, seq, 'Fx')
    
    meanCFxLeft = CFxLeft.mean()[0:7]
    meanCFxRight = CFxRight.mean()[0:7]
    
    netCFxLeft = CFxLeft.sum()[0:7]
    netCFxRight = CFxRight.sum()[0:7]
    
    meanMeanCFx = np.mean([meanCFxLeft, meanCFxRight], axis=0)
    meanNetCFx = np.mean([netCFxLeft, netCFxRight], axis=0)
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left')]
    #save left side data
    trdf.loc[placementIdx,'meanCFxInSegm'] = np.repeat(meanCFxLeft.values,4)
    trdf.loc[placementIdx,'netCFxInSegm'] = np.repeat(netCFxLeft.values,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right')]
    #save right side data
    trdf.loc[placementIdx,'meanCFxInSegm'] = np.repeat(meanCFxRight.values,4)
    trdf.loc[placementIdx,'netCFxInSegm'] = np.repeat(netCFxRight.values,4)

    #make placementIdx for average of sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq)]
    #save mean side data
    trdf.loc[placementIdx,'meanCFxInSegmSidesAvg'] = np.repeat(meanMeanCFx,8)
    trdf.loc[placementIdx,'netCFxInSegmSidesAvg'] = np.repeat(meanNetCFx,8)
    
    
    #Get total CFy information
    
    CFyLeft = get_NETforce_vs_time_binned_all_TOP(trpath, seq, 'Fy')
    CFyRight = get_NETforce_vs_time_binned_all_BOTTOM(trpath, seq, 'Fy')
    
    meanCFyLeft = CFyLeft.mean()[0:7]
    meanCFyRight = CFyRight.mean()[0:7]
    
    netCFyLeft = CFyLeft.sum()[0:7]
    netCFyRight = CFyRight.sum()[0:7]
    
    meanMeanCFy = np.mean([meanCFyLeft, meanCFyRight], axis=0)
    meanNetCFy = np.mean([netCFyLeft, netCFyRight], axis=0)
    
    #make placementIdx for left sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left')]
    #save left side data
    trdf.loc[placementIdx,'meanCFyInSegm'] = np.repeat(meanCFyLeft.values,4)
    trdf.loc[placementIdx,'netCFyInSegm'] = np.repeat(netCFyLeft.values,4)
    
    #make placementIdx for right sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right')]
    #save right side data
    trdf.loc[placementIdx,'meanCFyInSegm'] = np.repeat(meanCFyRight.values,4)
    trdf.loc[placementIdx,'netCFyInSegm'] = np.repeat(netCFyRight.values,4)
    
    #make placementIdx for average of sides
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq)]
    #save mean side data
    trdf.loc[placementIdx,'meanCFyInSegmSidesAvg'] = np.repeat(meanMeanCFy,8)
    trdf.loc[placementIdx,'netCFyInSegmSidesAvg'] = np.repeat(meanNetCFy,8)
    
    
    
    
    
    #Get single force type info - Fx
    
    thrustHiLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(trpath, seq, 'fpush')
    thrustHiRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(trpath, seq, 'fpush') 
    
    thrustLoLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(trpath, seq, 'fpull')
    thrustLoRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(trpath, seq, 'fpull')
    
    dragHiLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(trpath, seq, 'rpush')
    dragHiRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(trpath, seq, 'rpush')
    
    dragLoLeft = get_NETCFx_sourced_vs_time_binned_all_TOP(trpath, seq, 'rpull')
    dragLoRight = get_NETCFx_sourced_vs_time_binned_all_BOTTOM(trpath, seq, 'rpull')
    
    #calculate values of interest
    
    #means
    meanThrustHiLeft = thrustHiLeft.mean()[0:7]
    meanThrustHiRight = thrustHiRight.mean()[0:7]
    
    meanThrustLoLeft = thrustLoLeft.mean()[0:7]
    meanThrustLoRight = thrustLoRight.mean()[0:7]
    
    meanDragHiLeft = dragHiLeft.mean()[0:7]
    meanDragHiRight = dragHiRight.mean()[0:7]
    
    meanDragLoLeft = dragLoLeft.mean()[0:7]
    meanDragLoRight = dragLoRight.mean()[0:7]
    
    meanThrustHi = np.mean([meanThrustHiLeft, meanThrustHiRight], axis=0)
    meanThrustLo = np.mean([meanThrustLoLeft, meanThrustLoRight], axis=0)
    meanDragHi = np.mean([meanDragHiLeft, meanDragHiRight], axis=0)
    meanDragLo = np.mean([meanDragLoLeft, meanDragLoRight], axis=0)
    
    #net
    netThrustHiLeft = thrustHiLeft.sum()[0:7]
    netThrustHiRight = thrustHiRight.sum()[0:7]
    
    netThrustLoLeft = thrustLoLeft.sum()[0:7]
    netThrustLoRight = thrustLoRight.sum()[0:7]
    
    netDragHiLeft = dragHiLeft.sum()[0:7]
    netDragHiRight = dragHiRight.sum()[0:7]
    
    netDragLoLeft = dragLoLeft.sum()[0:7]
    netDragLoRight = dragLoRight.sum()[0:7]
    
    netThrustHi = np.mean([netThrustHiLeft, netThrustHiRight], axis=0)
    netThrustLo = np.mean([netThrustLoLeft, netThrustLoRight], axis=0)
    netDragHi = np.mean([netDragHiLeft, netDragHiRight], axis=0)
    netDragLo = np.mean([netDragLoLeft, netDragLoRight], axis=0)
    
    #peak values
    peakThrustHiLeft, peakTimeThrustHiLeft, peakIdxThrustHiLeft = get_NETCFx_source_peaks(thrustHiLeft,'fpush')
    peakThrustHiRight, peakTimeThrustHiRight, peakIdxThrustHiRight = get_NETCFx_source_peaks(thrustHiRight,'fpush')
    
    peakThrustLoLeft, peakTimeThrustLoLeft, peakIdxThrustLoLeft = get_NETCFx_source_peaks(thrustLoLeft,'fpull')
    peakThrustLoRight, peakTimeThrustLoRight, peakIdxThrustLoRight = get_NETCFx_source_peaks(thrustLoRight,'fpull')
    
    peakDragHiLeft, peakTimeDragHiLeft, peakIdxDragHiLeft = get_NETCFx_source_peaks(dragHiLeft,'rpush')
    peakDragHiRight, peakTimeDragHiRight, peakIdxDragHiRight = get_NETCFx_source_peaks(dragHiRight,'rpush')
    
    peakDragLoLeft, peakTimeDragLoLeft, peakIdxDragLoLeft = get_NETCFx_source_peaks(dragLoLeft,'rpull')
    peakDragLoRight, peakTimeDragLoRight, peakIdxDragLoRight = get_NETCFx_source_peaks(dragLoRight,'rpull')


    meanPeakThrustHi = np.mean([peakThrustHiLeft, peakThrustHiRight], axis=0)
    meanPeakThrustLo = np.mean([peakThrustLoLeft, peakThrustLoRight], axis=0)
    meanPeakDragHi = np.mean([peakDragHiLeft, peakDragHiRight], axis=0)
    meanPeakDragLo = np.mean([peakDragLoLeft, peakDragLoRight], axis=0)
    
    meanPeakTimeThrustHi = np.mean([peakTimeThrustHiLeft, peakTimeThrustHiRight], axis=0)
    meanPeakTimeThrustLo = np.mean([peakTimeThrustLoLeft, peakTimeThrustLoRight], axis=0)
    meanPeakTimeDragHi = np.mean([peakTimeDragHiLeft, peakTimeDragHiRight], axis=0)
    meanPeakTimeDragLo = np.mean([peakTimeDragLoLeft, peakTimeDragLoRight], axis=0)
    
    #amp and angle at peak values
    ampatPeakThrustHiLeft = np.zeros((7,))
    ampatPeakThrustHiRight = np.zeros((7,))
    ampatPeakThrustLoLeft = np.zeros((7,))
    ampatPeakThrustLoRight = np.zeros((7,))
    ampatPeakDragHiLeft = np.zeros((7,))
    ampatPeakDragHiRight = np.zeros((7,))
    ampatPeakDragLoLeft = np.zeros((7,))
    ampatPeakDragLoRight = np.zeros((7,))
    
    ampatPeakThrustHi = np.zeros((7,))
    ampatPeakThrustLo = np.zeros((7,))
    ampatPeakDragHi = np.zeros((7,))
    ampatPeakDragLo = np.zeros((7,))    
    
    angleatPeakThrustHiLeft = np.zeros((7,))
    angleatPeakThrustHiRight = np.zeros((7,))
    angleatPeakThrustLoLeft = np.zeros((7,))
    angleatPeakThrustLoRight = np.zeros((7,))
    angleatPeakDragHiLeft = np.zeros((7,))
    angleatPeakDragHiRight = np.zeros((7,))
    angleatPeakDragLoLeft = np.zeros((7,))
    angleatPeakDragLoRight = np.zeros((7,))
    
    angleatPeakThrustHi = np.zeros((7,))
    angleatPeakThrustLo = np.zeros((7,))
    angleatPeakDragHi = np.zeros((7,))
    angleatPeakDragLo = np.zeros((7,))

    
    
    
    for i in range(7):
        
        try:
            ampatPeakThrustHiLeft[i] = amp.iloc[int(peakIdxThrustHiLeft[i]), i]
        except:
            ampatPeakThrustHiLeft[i] = np.nan
        
        try:
            ampatPeakThrustHiRight[i] = amp.iloc[int(peakIdxThrustHiRight[i]), i]
        except:
            ampatPeakThrustHiRight[i] = np.nan
        
        ampatPeakThrustHi[i] = np.mean([abs(ampatPeakThrustHiLeft[i]),abs(ampatPeakThrustHiRight[i])])
        
        try:
            ampatPeakThrustLoLeft[i] = amp.iloc[int(peakIdxThrustLoLeft[i]), i]
        except:
            ampatPeakThrustLoLeft[i] = np.nan
        
        try:
            ampatPeakThrustLoRight[i] = amp.iloc[int(peakIdxThrustLoRight[i]), i]
        except:
            ampatPeakThrustLoRight[i] = np.nan
            
        ampatPeakThrustLo[i] = np.mean([abs(ampatPeakThrustLoLeft[i]),abs(ampatPeakThrustLoRight[i])])
            
        try:
            ampatPeakDragHiLeft[i] = amp.iloc[int(peakIdxDragHiLeft[i]), i]
        except:
            ampatPeakDragHiLeft[i] = np.nan
            
        try:
            ampatPeakDragHiRight[i] = amp.iloc[int(peakIdxDragHiRight[i]), i]
        except:
            ampatPeakDragHiRight[i] = np.nan
            
        ampatPeakDragHi[i] = np.mean([abs(ampatPeakDragHiLeft[i]),abs(ampatPeakDragHiRight[i])])
        
        try:
            ampatPeakDragLoLeft[i] = amp.iloc[int(peakIdxDragLoLeft[i]), i]
        except:
            ampatPeakDragLoLeft[i] = np.nan
          
        try:
            ampatPeakDragLoRight[i] = amp.iloc[int(peakIdxDragLoRight[i]), i]
        except:
            ampatPeakDragLoRight[i] = np.nan
        
        ampatPeakDragLo[i] = np.mean([abs(ampatPeakDragLoLeft[i]),abs(ampatPeakDragLoRight[i])])
    
    
    for i in range(7):
        try:
            angleatPeakThrustHiLeft[i] = angle.iloc[int(peakIdxThrustHiLeft[i]), i]
        except:
            angleatPeakThrustHiLeft[i] = np.nan
        
        try:
            angleatPeakThrustHiRight[i] = angle.iloc[int(peakIdxThrustHiRight[i]), i]
        except:
            angleatPeakThrustHiRight[i] = np.nan
            
        angleatPeakThrustHi[i] = np.mean([abs(angleatPeakThrustHiLeft[i]),abs(angleatPeakThrustHiRight[i])])
        
        try:
            angleatPeakThrustLoLeft[i] = angle.iloc[int(peakIdxThrustLoLeft[i]), i]
        except:
            angleatPeakThrustLoLeft[i] = np.nan
        
        try:
            angleatPeakThrustLoRight[i] = angle.iloc[int(peakIdxThrustLoRight[i]), i]
        except:
            angleatPeakThrustLoRight[i] = np.nan
            
        angleatPeakThrustLo[i] = np.mean([abs(angleatPeakThrustLoLeft[i]),abs(angleatPeakThrustLoRight[i])])
        
        try:
            angleatPeakDragHiLeft[i] = angle.iloc[int(peakIdxDragHiLeft[i]), i]
        except:
            angleatPeakDragHiLeft[i] = np.nan
            
        try:
            angleatPeakDragHiRight[i] = angle.iloc[int(peakIdxDragHiRight[i]), i]
        except:
            angleatPeakDragHiRight[i] = np.nan
        
        angleatPeakDragHi[i] = np.mean([abs(angleatPeakDragHiLeft[i]),abs(angleatPeakDragHiRight[i])])        
        
        try:
            angleatPeakDragLoLeft[i] = angle.iloc[int(peakIdxDragLoLeft[i]), i]
        except:
            angleatPeakDragLoLeft[i] = np.nan
            
        try:
            angleatPeakDragLoRight[i] = angle.iloc[int(peakIdxDragLoRight[i]), i]
        except:
            angleatPeakDragLoRight[i] = np.nan
            
        angleatPeakDragLo[i] = np.mean([abs(angleatPeakDragLoLeft[i]),abs(angleatPeakDragLoRight[i])])

    
    #area under curve
    areaThrustHiLeft = []
    areaThrustHiRight = []
    areaThrustLoLeft = []
    areaThrustLoRight = []
    areaDragHiLeft = []
    areaDragHiRight = []
    areaDragLoLeft = []
    areaDragLoRight = []
    
    meanAreaThrustHi = []
    meanAreaThrustLo = []
    meanAreaDragHi = []
    meanAreaDragLo = []
    
    for col in range(7):
        
        areaThrustHiLeft.append(np.trapz(thrustHiLeft[col],thrustHiLeft['% Period']))
        areaThrustHiRight.append(np.trapz(thrustHiRight[col],thrustHiRight['% Period']))
        meanAreaThrustHi.append(np.mean([areaThrustHiLeft[-1], areaThrustHiRight[-1]]))
        
        areaThrustLoLeft.append(np.trapz(thrustLoLeft[col],thrustLoLeft['% Period']))
        areaThrustLoRight.append(np.trapz(thrustLoRight[col],thrustLoRight['% Period']))
        meanAreaThrustLo.append(np.mean([areaThrustLoLeft[-1], areaThrustLoRight[-1]]))
        
        areaDragHiLeft.append(np.trapz(dragHiLeft[col],dragHiLeft['% Period']))
        areaDragHiRight.append(np.trapz(dragHiRight[col],dragHiRight['% Period']))
        meanAreaDragHi.append(np.mean([areaDragHiLeft[-1], areaDragHiRight[-1]]))
        
        areaDragLoLeft.append(np.trapz(dragLoLeft[col],dragLoLeft['% Period']))
        areaDragLoRight.append(np.trapz(dragLoRight[col],dragLoRight['% Period']))
        meanAreaDragLo.append(np.mean([areaDragLoLeft[-1], areaDragLoRight[-1]]))
    
    
    #save all values
    
    #thrust hi
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left') & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanThrustHiLeft.values
    trdf.loc[placementIdx,'netCFxByMech'] = netThrustHiLeft.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakThrustHiLeft
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustHiLeft.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustHiLeft.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustHiLeft)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustHiLeft
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustHiLeft
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustHiLeft
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustHiLeft
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right') & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanThrustHiRight.values
    trdf.loc[placementIdx,'netCFxByMech'] = netThrustHiRight.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakThrustHiRight
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustHiRight.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustHiRight.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustHiRight)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustHiRight
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustHiRight
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustHiRight
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustHiRight
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanThrustHi,2)
    trdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netThrustHi,2)
    trdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakThrustHi,2)
    trdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanThrustHi),2)
    trdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netThrustHi),2)
    trdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakThrustHi),2)
    trdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeThrustHi,2)
    trdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaThrustHi,2)
    trdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakThrustHi,2)
    trdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakThrustHi,2)

    #thrust lo
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left') & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanThrustLoLeft.values
    trdf.loc[placementIdx,'netCFxByMech'] = netThrustLoLeft.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakThrustLoLeft
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustLoLeft.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustLoLeft.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustLoLeft)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustLoLeft
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustLoLeft
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustLoLeft
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustLoLeft
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right') & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanThrustLoRight.values
    trdf.loc[placementIdx,'netCFxByMech'] = netThrustLoRight.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakThrustLoRight
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanThrustLoRight.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netThrustLoRight.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakThrustLoRight)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeThrustLoRight
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaThrustLoRight
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakThrustLoRight
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakThrustLoRight
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['forceType'][i] == 'thrust') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanThrustLo,2)
    trdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netThrustLo,2)
    trdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakThrustLo,2)
    trdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanThrustLo),2)
    trdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netThrustLo),2)
    trdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakThrustLo),2)
    trdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeThrustLo,2)
    trdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaThrustLo,2)
    trdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakThrustLo,2)
    trdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakThrustLo,2)
    
    #drag hi
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left') & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanDragHiLeft.values
    trdf.loc[placementIdx,'netCFxByMech'] = netDragHiLeft.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakDragHiLeft
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragHiLeft.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragHiLeft.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragHiLeft)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragHiLeft
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragHiLeft
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragHiLeft
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragHiLeft
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right') & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanDragHiRight.values
    trdf.loc[placementIdx,'netCFxByMech'] = netDragHiRight.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakDragHiRight
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragHiRight.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragHiRight.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragHiRight)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragHiRight
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragHiRight
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragHiRight
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragHiRight
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'hi')]
    trdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanDragHi,2)
    trdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netDragHi,2)
    trdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakDragHi,2)
    trdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanDragHi),2)
    trdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netDragHi),2)
    trdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakDragHi),2)
    trdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeDragHi,2)
    trdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaDragHi,2)
    trdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakDragHi,2)
    trdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakDragHi,2)
    
    
    #drag lo
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'left') & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanDragLoLeft.values
    trdf.loc[placementIdx,'netCFxByMech'] = netDragLoLeft.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakDragLoLeft
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragLoLeft.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragLoLeft.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragLoLeft)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragLoLeft
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragLoLeft
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragLoLeft
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragLoLeft
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['side'][i] == 'right') & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMech'] = meanDragLoRight.values
    trdf.loc[placementIdx,'netCFxByMech'] = netDragLoRight.values
    trdf.loc[placementIdx,'peakCFxByMech'] = peakDragLoRight
    trdf.loc[placementIdx,'absMeanCFxByMech'] = abs(meanDragLoRight.values)
    trdf.loc[placementIdx,'absNetCFxByMech'] = abs(netDragLoRight.values)
    trdf.loc[placementIdx,'absPeakCFxByMech'] = abs(peakDragLoRight)
    trdf.loc[placementIdx,'peakTimeCFxByMech'] = peakTimeDragLoRight
    trdf.loc[placementIdx, 'areaUnderCurveByMech'] = areaDragLoRight
    trdf.loc[placementIdx, 'ampAtPeakByMech'] = ampatPeakDragLoRight
    trdf.loc[placementIdx, 'angleAtPeakByMech'] = angleatPeakDragLoRight
    placementIdx = [i for i in range(len(trdf['sequence'])) if (trdf['sequence'][i] == seq) & (trdf['forceType'][i] == 'drag') & (trdf['pressType'][i] == 'lo')]
    trdf.loc[placementIdx,'meanCFxByMechSidesAvg'] = np.repeat(meanDragLo,2)
    trdf.loc[placementIdx,'netCFxByMechSidesAvg'] = np.repeat(netDragLo,2)
    trdf.loc[placementIdx,'peakCFxByMechSidesAvg'] = np.repeat(meanPeakDragLo,2)
    trdf.loc[placementIdx,'absMeanCFxByMechSidesAvg'] = np.repeat(abs(meanDragLo),2)
    trdf.loc[placementIdx,'absNetCFxByMechSidesAvg'] = np.repeat(abs(netDragLo),2)
    trdf.loc[placementIdx,'absPeakCFxByMechSidesAvg'] = np.repeat(abs(meanPeakDragLo),2)
    trdf.loc[placementIdx,'peakTimeCFxByMechSidesAvg'] = np.repeat(meanPeakTimeDragLo,2)
    trdf.loc[placementIdx, 'areaUnderCurveByMechSidesAvg'] = np.repeat(meanAreaDragLo,2)
    trdf.loc[placementIdx, 'ampAtPeakByMechSidesAvg'] = np.repeat(ampatPeakDragLo,2)
    trdf.loc[placementIdx, 'angleAtPeakByMechSidesAvg'] = np.repeat(angleatPeakDragLo,2)
    
    
    
#combine to one df
allfishdf = pd.concat([bgdf,trdf],ignore_index=True)
allfishdf = allfishdf.drop(['index'],axis=1)
'''
allfishdf.to_excel('D:/fish pressure_v2/fishData_compiledForStats_v6.xlsx')
allfishdf.to_csv('D:/fish pressure_v2/fishData_compiledForStats_v6.csv')
'''
