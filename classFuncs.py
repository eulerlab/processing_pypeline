# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 13:57:59 2016

@author: andre
"""
#necessary libraries
import matplotlib.pyplot as plt
#from matplotlib import colors as clrs
import numpy as np
from scipy import signal
import math
import seaborn as sns
import os
import tkinter as tk
import tkinter.filedialog as fd
import pandas as pd
import time
import pycircstat as circ
import peakutils as pk
from configparser import ConfigParser

sns.set_style("white")

os.chdir("E:\\github\\processing_pypeline\\" )
import readScanM as rsm


def avg_matrix(matrix=None,grouping=None):
    """function to get the average response to each variant of a certain stimulus
    (e.g. diferent movement directions of moving bars)
    An example on how the grouping list should be arranged:
    grouping = [[0,8,16],[1,9,17],[2,10,18],[3,11,19],[4,12,20],[5,13,21],[6,14,22],[7,15,23]]
    where each integer corresponds to the trial number, if no grouping is 
    provided, than average the entire matrix
    inputs: 
        matrix: timepointsXtrials
        grouping: list of indices. see description above
    outputs:
        dsMatrix: new matrix with average of responses for each group
    """
    if grouping is not None:
        indices = grouping
        dsRows = np.size(matrix,axis=0)
        dsMatrix=np.zeros((dsRows,len(indices)))
    
        for j in range(0,len(indices)):          
            meanDatum = np.mean(matrix[:,indices[j]],axis=1)
            dsMatrix[:,j] = meanDatum
    
        
    else:
        dsMatrix = np.mean(matrix,axis=0)
    
    return dsMatrix

def calc_area(rois,cellnum=1):
    """calculate the area of the rois, based on the number of pixels
    inputs: 
        rois: numpy array containing all the ROIS - filled with zeros where there is no ROI
        and a number for every pixel corresponding to a specific ROI
        cellnum: the id of the cell/ROI to be measured
    outputs:
        row: row indexes for all values found
        col: column indexes for all values found
        pix: the number of pixels of a particular roi.
    """
    
    indx = np.where(np.array(rois)==cellnum)
    row = indx[0]
    col = indx[1]    
    pix=len(row)
    
    
    return pix,row,col
    
def calc_delay(stimDelay=17.0,sampRate=7.8):
    """function to calculate the delay in between the trigger and stimulator
    it requires samprate(hz) and stimdelay(ms). It outputs the delay in bins,
    this number can be used to move the stimulus and/or response arrays so that
    they are aligned."""
    #stim = np.empty(shape=(1,int(p["chirpDur_s"]/p["durFr_s"])))
    delay = stimDelay/1000. #convert delay to seconds
    delay = delay / (1/sampRate) # now converto to bins
    delay = np.floor(delay)+1
    return delay
    
def calcGeoC(rois,cellnum=1):
    """calculate the geometric centre of the rois
    input:
        rois: 2D matrix containing numbers on the location of each roi
        cellNum: the number representing a ROI
    return:
        rowMed: the roi median value at rows
        colMed: the roi median value at columns"""
    indx = np.where(rois==cellnum)
    rowMed = np.median(indx[0])
    colMed = np.median(indx[1])
    return rowMed, colMed

def clean_field(filePath,minQual):
    cleanList=list()
    qualAll = list()
    hdfStore = pd.HDFStore(filePath)
    keys = hdfStore.keys()
    for cell in keys:
#        print (cell + ": ")
        data = pd.read_hdf(hdfStore,cell)
        data = data.transpose()
        qual = data["qualIndex"].dropna()
        qualAll.append(qual)
#        print(qual)
        if qual.values[0] > minQual:
            cleanList.append(cell)
    
    return cleanList,qualAll

    
def cluster_distribution(hdFilePath, keys ):
    """function to go throught the chirp hdf5 file 
    (all cell responses to chirp stimulus) and get the information
    about which cluster each cell belongs to.
    The second input is a list with the keys inside the hdf file"""
    
    #oldL = None
    oldL = dir()
    clusterList = list()
    for cell in keys:
#        hdfStore = pd.HDFStore(hdFilePath)
        data = pd.read_hdf(hdFilePath,cell)
        data = data.transpose()
        cluster = int(max(data["clusIndx"]["chirp"]))
        if "c"+str(cluster) not in locals():
            exec("c"+str(cluster)+"=pd.DataFrame()")
            clusterList.append("c"+str(cluster))
        series = pd.Series(data["medianTrace"]["chirp"],name=(cell[1:]))
        exec("c"+str(cluster)+"=c"+str(cluster)+".append(series.dropna())")
        #exec("output['c"+str(cluster)+"']=output['c"+str(cluster)+"'].append(series.dropna())")
        
    output=dict()
    for c in clusterList:
        exec("output[c]="+c)
  
    return output

def coor2series(headerDictionary,key,sufix=""):
    coor = float(headerDictionary[key])
    output = pd.Series(coor,name=(key,sufix))
    return output   
    
    
def correlate_traces(trace1,trace2):
    """cross correlation of two traces. 
    If trace1=trace2, than this is autocorrelation.
    inputs:
        trace1
        trace2
    returns:
        maxCorr: the maximum correlation found
        maxLag: the shift in between traces for
        the point of maximum correlation
        zeroCorr: the correlation when traces are aligned (middle)
        zeroLag: lag at zeroCorr (0)
            """
    correlated = plt.xcorr(x=trace1,y=trace2,
                              normed=True, usevlines=False,
                              maxlags=None)    
    
    maxCorr = np.max(correlated[1])
    maxLag = correlated[0][correlated[1] == maxCorr]
    maxLag = maxLag[0]
    
    zeroLag = np.where(correlated[0] == 0)    
    zeroLag = int(zeroLag[0])
    zeroCorr = correlated[1][zeroLag]    
    return maxCorr,maxLag,zeroCorr,zeroLag
    


def create_trigger_trace(trace,traceTime,triggerTime):
    """create a trigger trace with the same length as the stimulus trace,
        It has 1s as trigger points and 0s otherwise. The function also outputs
        the indexes where the triggers ocurred"""

    trigger = np.zeros(len(trace))
    triggerInd=list()
    for time1 in triggerTime:
        #print(time1)
        dummie = np.where(time1>=traceTime)
        try:
            dummie = dummie[0][-1]
            triggerInd.append(dummie)
            trigger[dummie]=1
        except IndexError:
            print ("trigger before traceTime")
    
    

    return trigger,triggerInd

def create_2d_gaussian_window(rows,columns,sd):
    """create a 2d gaussian window that can be used to filter matrices
    returns:
        gauss: rowsXcolumns matrix containing the gaussian distrib.
    """
    
    gauss = signal.gaussian(M=columns,std=sd,sym=True)
    gauss = np.tile(gauss,(rows,1)) 
    gauss = np.dot(np.transpose(gauss),gauss)
    return gauss
    
def create_bg_stim(sampFreq=7.8, greenFirst = 0):
    """create a trace representing the green/blue full field
    stimulation. the green part of the stimulation is represented 
    as double as high as the blue part.
    returns:
        stim: stimulus trace
        tStim: time vector in sec.
    
    """
    p = {"nTrials"         : 3, 
     "dxStim_um"       : 1000,   # Stimulus size
     "durFr_s"         : 1/60.0, # Frame duration
     "nFrPerMarker"    : 3,
     "RGB_green"       : (0,255,0),
     "RGB_blue"        : (255,0,0),
     "begOff"          : 3.0,
     "begOn"           : 3.0}
     
    begOff = np.zeros(shape=(int(np.ceil(p["begOff"]/p["durFr_s"])),1)) 
    begOn = np.ones(shape=(int(np.ceil(p["begOn"]/p["durFr_s"])),1)) 
    green = begOn*2
    
    if greenFirst==1:
        stim = np.concatenate((green,begOff,begOn,begOff),axis=0)
    else:
        stim = np.concatenate((begOn,begOff,green,begOff),axis=0)
    
#    tStim = np.linspace(0,int(len(stim)/sampFreq),len(stim))
    tStim = np.linspace(0,len(stim),len(stim)) 
    tStim = tStim*(p["durFr_s"])
    return stim, tStim, greenFirst
    
def create_chirp_stim(sampFreq=7.8,stimDelay=0):
    """This function recreates one trial of 
    the chirp stimulus. it has a hardcoded dictionary taken from the script
    to generate the chirp in QDSpy. It needs to be changed if the user used
    QDSstim to stimulate the retina. The function takes sampling frequency(Hz) 
    and stimulator delay(bins) as inputs"""
    
    p = {"nTrials"         : 4, 
     "chirpDur_s"      : 8.0,    # Rising chirp phase
     "chirpMaxFreq_Hz" : 8.0,    # Peak frequency of chirp (Hz)
     "ContrastFreq_Hz" : 2.0,    # Freqency at which contrast 
                                 # is modulated
     "tSteadyOFF_s"    : 3.0,    # Light OFF at beginning ...
     "tSteadyOFF2_s"   : 2.0,    # ... and end of stimulus
     "tSteadyON_s"     : 3.0,    # Light 100% ON before and 
                                 # after chirp
     "tSteadyMID_s"    : 2.0,    # Light at 50% for steps
     "IHalf"           : 127,
     "IFull"           : 254,
     "dxStim_um"       : 1000,   # Stimulus size
     "StimType"        : 2,      # 1 = Box, 2 = Circle/
     "durFr_s"         : 1/60.0, # Frame duration
     "nFrPerMarker"    : 3}
     
    #stim = np.empty(shape=(1,int(p["chirpDur_s"]/p["durFr_s"])))
    delay = stimDelay
    
    begOff = np.zeros(shape=(int(np.ceil(p["tSteadyOFF2_s"]/p["durFr_s"])),1))
    endOff = np.zeros(shape=(int(np.ceil(p["tSteadyOFF_s"]/p["durFr_s"])),1))
         
    begOn = np.ones(shape=(int(p["tSteadyON_s"]/p["durFr_s"]),1))*p["IFull"]
    midOn = np.ones(shape=(int(p["tSteadyMID_s"]/p["durFr_s"]),1))*p["IHalf"]
    
    nPntChirp       = int(p["chirpDur_s"] /p["durFr_s"])
    K_HzPerSec      = p["chirpMaxFreq_Hz"] /p["chirpDur_s"]  # acceleration in Hz/s
    Intensity=list()
    for iT in range(nPntChirp):
        t_s       = iT*p["durFr_s"] # in ms
        Intensity.append(math.sin(math.pi *K_HzPerSec *t_s**2) *p["IHalf"] +p["IHalf"])
    
    Intensity = np.reshape(Intensity,[len(Intensity),1])    
    
    Intensity2=list()
    for iPnt in range(nPntChirp):
        t_s       = iPnt*p["durFr_s"]
        IRamp     = int(p["IHalf"] *t_s /p["chirpDur_s"])
        Intensity2.append(math.sin(2*math.pi *p["ContrastFreq_Hz"] *t_s) *IRamp +p["IHalf"])
        
    Intensity2 = np.reshape(Intensity2,[len(Intensity2),1])    
    
    stim=np.concatenate((begOff,begOn,endOff,midOn,Intensity,midOn,Intensity2,midOn,endOff),
                    axis=0)
    
    stim = stim/np.max(stim)

    tStim = np.linspace(0,len(stim),len(stim)) 
    tStim = tStim*(p["durFr_s"])
#    tStim = np.linspace(0,int(len(stim)/sampFreq),len(stim))
#    #convert time vector to secs    
#    tStim = tStim*p["durFr_s"]
    #remove delay
    stim = stim[delay:]
    tStim = tStim[0:len(tStim)-delay]
    return stim, tStim

def create_dict(traceMat = None,
                medianTrace = None, sampFreq = None,
                stimulusTrace = None,taxis = None,
                t_stim=None,obs=None):
   
    """create a dictionary. This is a legacy function used before the dataFrames with
    pandas were implemented"""
    
    tempDict = {"traces" : traceMat,
               "median" : medianTrace,
               "sampfreq" : sampFreq,
               "stimulus" : stimulusTrace,
               "taxis" : taxis,
               "t_stim" : t_stim,
               "obs" : obs}
    return tempDict
                    
def create_ds_stim(sampFreq=7.8,
                   p = {"nTrials"         : 3, 
                        "DirList"         : [0,180, 45,225, 90,270, 135,315],       
                        "vel_umSec"       : 1000.0, # speed of moving bar in um/sec
                        "tMoveDur_s"      : 4.0,    # duration of movement (defines distance
                                                    # the bar travels, not its speed)
                        "barDx_um"        : 300.0,  # bar dimensions in um
                        "barDy_um"        : 1000.0,
                        "bkgColor"        : (0,0,0),       # background color
                        "barColor"        : (255,255,0), # bar color
                        "durFr_s"         : 1/60.0, # Frame duration
                        "nFrPerMarker"    : 3,
                        "begOff"          : 3*(1/60.),
                        #     "endOff"          : 3.0
                        }):
    """creates a trace indicating moment of stimulation for each trial for 
    moving bar stimuli. The time vector is given seconds, it takes delay in 
    bins as input as well as recording samp freq"""
    
    begOff = np.zeros(shape=(int(np.ceil(p["begOff"]/p["durFr_s"])),1))
    begOn =  np.ones(shape=(int(np.ceil(p["tMoveDur_s"]/p["durFr_s"])),1))
#    endOff = np.zeros(shape=(int(np.ceil(p["endOff"]/p["durFr_s"])),1))
    
#    stim = np.concatenate((begOn,begOff,endOff),axis = 0)
    stim = np.concatenate((begOff,begOn),axis = 0)
#    tStim = np.linspace(0,int(len(stim)/sampFreq),len(stim))
    tStim = np.linspace(0,len(stim),len(stim)) 
    tStim = tStim*(p["durFr_s"])
    directions = p["DirList"]
    screenDur = p["tMoveDur_s"]
    return stim,tStim, directions, screenDur

def create_field_plot(fh=None,ax=None):
    if fh==None:
        fh = plt.figure()
    if ax==None:
        ax=plt.gca()
            
    plt.hlines(y=0,xmin=-5,xmax=5,linestyle="--")
    plt.vlines(x=0,ymin=-5,ymax=5,linestyle="--")      
    return fh,ax
            
def create_step(sampFreq=7.8,sizes=1,onTime=2.0,offTime=2.0):
    """create a trace representing the step (ON/OFF, flicker) stimulus
    returns:
        stim: stimulus trace.
        tStim: time vector in seconds"""
    p = {"durFr_s" : 1/60.0}
    
    begOff = np.zeros(shape=(int(np.ceil(offTime/p["durFr_s"])),1))
    begOn =  np.ones(shape=(int(np.ceil(onTime/p["durFr_s"])),1))
    
    stim=np.empty(shape=(1,1))
    for i in range(0,sizes):
        stim=np.concatenate((stim,begOn,begOff),axis=0)
#    endOff = np.zeros(shape=(int(np.ceil(p["endOff"]/p["durFr_s"])),1))
    tStim = np.linspace(0,len(stim),len(stim)) 
    tStim = tStim*(p["durFr_s"])
    return stim[1:], tStim[0:-1]

def delay_ds(roix=[1],framex=64,framexmicro=110,timeToCross=4.0,sampRate=7.8):
    """calculate how many frames it takes for the front part of the moving
    bar to reach a certain roi, this can be used for the on/off index 
    calculation, as ROIs on the opposite edge of the stimulus would only be
    stimulated later"""
    
    #number of frames to cross based on recording sampRate
    fCross = round(timeToCross*sampRate)
    
    maxx = max(roix)        
    minx = min(roix)
    

    
    #time to cross one pixel
    pCross = timeToCross/framex
    #time to cross the Roix
    rCross = pCross*(maxx-minx)
    
    firstX = int(round((minx*fCross)/framex))
    lastX  = round((maxx*fCross)/framex)
    lastX  = int(lastX+round(rCross*sampRate))
    
    return firstX, lastX
    

    
def detrend(sampRate=8,cutoff=0.1,length=125,trace=None):
    """detrend traces using a high cutoff filter - this is different from
    igor detrending, as there the method is applied to the data before trace
    extraction"""


    # Configuration.
    fS = sampRate  # Sampling rate.
    fH = cutoff  # Cutoff frequency.
    N = length  # Filter length, must be odd.

    # Compute sinc filter.
    h = np.sinc(2 * fH / fS * (np.arange(N) - (N - 1) / 2.))

    # Apply window.
    h *= np.hamming(N)

    # Normalize to get unity gain.
    h /= np.sum(h)

    # Create a high-pass filter from the low-pass filter through spectral inversion.
    h = -h
    h[(N - 1) / 2] += 1

#    print(h)

    # Applying the filter to a signal s can be as simple as writing
    s = np.convolve(trace, h)
    return s

def direction_selectivity(matrix=None):
    """check the direction selectivity of RGC cells that were stimulated
    with moving bar stimulus.
    use single value decomposition of the average response
    matrix (timeXcondition).
    inputs:
        matrix: time X condition(moving bar angle) 
    returns:
        normResp: 
            normalized response 
        dsVector:
            directive selection vector
        tc:
           time component  
    
    """    
    
    #SVD analysis to determine direction and orientation selectivity
    U,S,V = np.linalg.svd(matrix)
    sv =  np.sign(np.mean(np.sign(V[:,0])))
                
    if (np.mean((-1*U[:,0]-np.mean(matrix[:,0],axis=0))**2) < \
        np.mean((U[:,0]-np.mean(matrix[:,0],axis=0))**2)):
        
        su = -1
    else:
        su = 1
                    
    if sv==1 and su==1:
        s = 1
    elif sv==-1 and su==-1:
        s = -1
    elif sv==1 and su==-1:
        s = 1
    elif sv==0:
        s = su
    else:
        s = 1
                    
    dsVector = s * V[:,0]
    dsVector = dsVector - np.min(dsVector)
    dsVector = dsVector/max(dsVector);
                
    tc = s * U[:,0]
    tc = tc - np.mean(tc[0:7])
    tc = tc/np.max(abs(tc))
                
    normResp = np.zeros(shape = np.shape(matrix))
    for j in range(0,np.size(matrix,axis=1)):
        normResp[:,j] = np.transpose(tc)*matrix[:,j]

    return normResp,dsVector,tc
    
def field_loc(ODseries,filePath="",pattern=pd.Series({"-x":"dorsal",
                                                      "+x":"ventral",
                                                      "-y":"nasal",
                                                      "+y":"temporal"})):
    
    if filePath=="":
        root = tk.Tk()
        root.withdraw()
                 
        filePath = fd.askopenfilenames(initialdir= filePath,
                                    filetypes=[('header files', '.smh')],
                                    title = ["choose the file containing the field recording"])
        filePath=filePath[0]


        
    fieldHeader = rsm.read_in_header(filePath = filePath)
    
    xfield = float(fieldHeader["XCoord_um"]) - ODseries["x"]
    yfield = float(fieldHeader["YCoord_um"]) - ODseries["y"]
    #zfield = float(fieldHeader["ZCoord_um"]) - ODseries["z"]

    if xfield >= 0:
        fx = pattern["+x"]
    else:
        fx = pattern["-x"]

    if yfield >= 0:
        fy = pattern["+y"]
    else:
        fy = pattern["-y"]

#    if zfield >= 0:
#        fz = "aboveOD"
#    else:
#        fz = "belowOD"
#    ind1 = pd.MultiIndex.from_tuples([(xfield,yfield)], names=[fx, fy])
    xSeries = pd.Series(xfield,index=["x"],name=fx)
    ySeries = pd.Series(yfield,index=["y"],name=fy)
    fieldOut=pd.DataFrame([xSeries,ySeries])
#    fieldOut = pd.DataFrame(data = [xfield,yfield],#,[zfield]], 
#                            index = ["x","y"],#["x","y"],#,"z"],
#                            columns = [[fx],[fy] ])#+" "+fz])
#    
    return fieldOut
   

def get_folder_tree(folder):
    treeList = list()
    
    for path, subdirs, files in os.walk(folder):
#        print(path)
#        print(subdirs)
#        print (files)
        for name in files:
            #print (os.path.join(path, name))
            treeList.append(os.path.join(path, name))
        
    return treeList

    
def get_labels(label="trial",dataFrame=None):
    outLabels = list()
    for key in dataFrame.keys():
        if label in key[0]:
            outLabels.append(key)
    return outLabels
    
def get_color_resp_max_value(respMatrix=None,interval=[5,10]):
    """"""    
    respMatrix=respMatrix[interval[0]:interval[1],:]
    maxVal=np.max(respMatrix**2)
    return maxVal
  
def get_cluster_histogram(dictionary = None,plot=0,fh=None,export=0,fn=None):
    """this function uses the "cluster_distribution" output
    to get the number of cells in each cluster,the cluster indentification
    and the identification of the cells in each cluster
    input:
        dictionary: a dictionary containing cluster id as the key and the number of cells as items
        plot: flag to generate a plot of the distribution
        fh: figure handle
        sh: subplot handle
        export: flag to indicatte whether the figure should be saved to disk or not
        fn: filename for the figure. If None, figure is saved on desktop with date_time as filename
    outputs:
        labels: the cluster ids
        numbers: the number of cells in each cluster
        cells: The cells belonging to each cluster
        if plot is 1 than:
            fh: Figure handle
            sh: subplot handle
            
            """
    labels = list()
    numbers = list()
    cells = list()
    for key in dictionary.keys():
        labels.append(key)
        numbers.append(len(dictionary[key]))
        cells.append(list(dictionary[key].transpose().keys()))
            

    if plot==1:
        sns.set_style("white")
        if fn == None:
            userhome = os.path.expanduser('~')
            desktop = userhome + '\\Desktop\\'
            currTime = time.localtime()
            currTime = str(currTime[0])+"_"+str(currTime[1])+"_"+ \
                       str(currTime[2])+"_"+str(currTime[3])+"_"+ \
                       str(currTime[4])+"_"+str(currTime[5])
            fn = desktop + currTime
            
        if fh == None:
            fh=sns.plt.figure()
            
        labels = [s.strip('c') for s in labels]
        plt.bar(left = np.subtract(np.array(labels,dtype="float"),[0.4]*len(labels)),height = numbers)
        sns.despine(top=True,right=True)
        plt.gca().set(xlabel='Cluster id', ylabel='Observations')
        #plt.gca().set(xticks=range(1,max(map(int,labels))))
        plt.gca().set(xticks=range(1,51))
#        plt.legend(labels)
        if export==1:
            fh.savefig(filename=fn, dpi=600)
        return labels, numbers, cells,fh
        
    else:
        fh=None  
        return labels, numbers, cells

def get_spikes(traceTime,triggertime):
    pass

def get_traces_n_triggers(allData):
    trace = allData.transpose()["normTrace"].dropna()
    trace = np.array(trace)
    trace = trace.flatten()
    triggerInd = allData.transpose()["trigInd"].dropna()
    trigger = allData.transpose()["triggerTrace"].dropna()

    return  trace,triggerInd,trigger

def get_vel_trace(trace):
    velTrace=trace[:]
                    
    #get first derivative
    velTrace=np.diff(velTrace,n=1,axis=0)
    #convolve vel trace with calcium signals
    velTrace = np.multiply(trace[1:],velTrace)
    #get normalized trace [0 to 1] to calculate peak
    normTrace = velTrace/max(abs(velTrace))
    #calculate SD
    sd = np.median(abs(normTrace),axis=0)/0.6745
    return velTrace,normTrace,sd

def get_peaks_above_sd(trace,sd,onlypos=1):
    #get indexes of peaks above 1 SD (findpeaks in matlab)
    indexes=pk.indexes(y=abs(trace),thres=float(sd),min_dist=1)
    if onlypos==1:
        #get only the positive spikes
        indexes = indexes[trace[indexes]>sd]
#        indexes=pk.indexes(y=trace,thres=float(sd),min_dist=1)
#        val=trace[indexes]>=sd
#        indexes=indexes[val]#.values]#indexes[val.values]
#    else:
        
#        indexes=pk.indexes(y=abs(trace),thres=float(sd),min_dist=1)
    return indexes
                        
def grab_baseline(trace=None,traceTime=None,triggerTime=None):
    """grab the baseline trace, aka the signal recorded before the 
    first trigger happens
    inputs:
        trace: the recorded trace - all trials in one trace
        traceTime: the time vector for the recorded trace
        triggerTime: list or array with the trigger timepoints
    outputs:
        baseline: snip of trace before the first trigger
        """
    
    #create binary array to indicate trigger points
    trigger,triggerInd = create_trigger_trace(trace,traceTime,triggerTime)    
    #catch the baseline    
    baseline = trace[0:triggerInd[0]]
    if len(baseline) is 0:
        baseline=0
        print("there was no baseline before 1st trigger")
    return baseline


def meanperc(datum=None,preStimDur=250,stimDur=500,
               postStimDur=250,zThresh=2.5,
               activThresh=0.5,areaThresh=0.7,bins = 4):
                   
    """function get mean cell response and percentiles.Also to divide the response
    period into equally sized intervals and get the area below the curve in those
    intervals"""
    
    resp = list()
    for i in range(0,len(datum)):
        # the data is already plotted as z-scores, so we look for variations in
        #between the maximum and minimum point that are above 2.5, which means 2.5SDs    
        ##if the difference between baseline and stimulus phase is bigger than 2.5
        ##at any point, than this is a cell response.
        ##if np.any((abs(stim[i])+abs(minima[i]))>1):
        if np.any((np.max(datum[i])-np.min(datum[i])) >= zThresh):
            resp.append(1)
        else:
            resp.append(0)
        
    #get the percentage of cell reponse give the number of trials.        
    respScore=np.sum(resp)/float(i+1)

    #Get the mean value of all traces
    meanDatum = np.mean(datum,axis=0)
    
    #get percentiles
    percentiles = np.percentile(datum,[5.,25.,50.,75.,95.],axis=0)
    
    #normalize data so that minimum point is 0
    normMeanDatum = meanDatum-np.min(meanDatum)
    
    #normalize so that maximum point is 1 --> this way the maximum response area is 
    #equal to the number of bins. This will be useful later, to know cell activity
    #without having to visually inspect cell by cell
    normMeanDatum = normMeanDatum/np.max(normMeanDatum)

    if respScore >= activThresh:
        totalTime = preStimDur+stimDur+postStimDur
        #divide the triggered window into 4 equally sized subwindows 
        intervals = np.linspace(start = 0, stop= totalTime, num=bins+1, dtype=int)
                     #[0,250,500,750,1000]

        #loop to go trhough all intervals and get responses from cells
        subresponses=list()
        for i in range (1,len(intervals)):
            dummie = np.trapz(normMeanDatum[intervals[i-1]:intervals[i]])    
            subresponses.append(np.divide(dummie,(intervals[i]-intervals[i-1])))        
        
        return (meanDatum, normMeanDatum, respScore, subresponses,percentiles) #,classStr)
    else:
        #classStr="no classification"
        subresponses=None
        return(meanDatum,normMeanDatum,respScore,
               subresponses, percentiles)#,classStr)    

def normalize(respMatrix=None):
    """normalize the data so max value of data  is 1.
    Matrix should be timeXrepetitions"""
    #median = np.median(respMatrix,axis=0)
    #respMatrix1=respMatrix-abs(median)
    maxval=np.max(abs(respMatrix),axis=0)
    respMatrix = np.divide(respMatrix,np.transpose(maxval))
    return respMatrix

def n_max_correlations(trace,dictionary,n_max):
    """function to grab the a specific number of correlation values, 
    in decreasing order. If n_max=3 the function will grab the values of the 3
    highest correlations at time lag 0"""
    allCorr = list()
    allClus = list()
    allClusTrace = list()
    corr=0
    corrClus=""
    clusTrace = 0
    keyUsed=list()
    for i in range (0,n_max+1):
        #        #set things back to 0
#        corr=0
#        corrClus=""
        newItem=0
        for key in dictionary.keys():
            #grab cluster trace

            if key[0] =="c" and key not in keyUsed:
                dummie2 = dictionary[key]["chirpMean"][0][0][0]
                #normalize cluster trace so they max(|trace|)=1
                dummie2 = normalize(dummie2)                

                #get the correlation at zero lag
                tempCorr = np.correlate(trace,dummie2,"same")
                zeroCorr = tempCorr[int(len(tempCorr)/2)]
                
                #if the current correlation if bigger than the previous one, 
                #set it as the new correlation.
                if zeroCorr > corr:
                    #print(zeroCorr)
                    #print(key)
                    newItem=1
                    corr = zeroCorr
                    corrClus = key
                    #raw cluster trace 
                    clusTrace = dummie2
                    keyUsed.append(key)
        if newItem==1:
            #once the loop is done, append the max correlation and 
            #the cluster name to lists
        
            allCorr.append(corr)
            allClus.append(corrClus)
            allClusTrace.append(clusTrace)
        
        #remove the key that originated the best correlation in order to 
        #get the 2nd best one (and os on for the number of loops)
#        try:
#            #print(corrClus)
#            dictionary.pop(corrClus)
#        except KeyError:
#            pass
#        #set things back to 0
        corr=0
        corrClus=""
        
    return allCorr, allClus, allClusTrace

 
def plot2cells_chirp(absPath1="",absPath2="",cell1="",cell2="",
                     clusterFlag=1,fh=None,sbHandle=None,sb1Handle=None):
    """compare median responses of two cells responding to chirp.
        inputs:
            absPath1: absolute path to the processed2 data (after raw2panda function)    
            cell1: the group to be read in the hdf5 file in absPath1. eg. "/cell18"
            absPath2: absolute path to the processed2 data (after raw2panda function)    
            cell2: the group to be read in the hdf5 file in absPath2. eg. "/cell35"
            clusterFlag: if == 1, function will plot the cluster traces that best correlated
            with each cell.
            fh: figure handle, can be used to add the plots to an existing figure
            sbhandle: subplot handles, can be used to add the plots to existing axes
            sb1handle: subplot handle cand be used to add the cluster traces to an existing axes
        returns:
            fh: figure handle
            [sh,sh1]: trace subplot handles
            [sh2,sh3]: cluster trace subplot handles
    """
    if absPath1 !="":
        data1 = pd.read_hdf(absPath1,cell1)
    
    else:
        print("no path given. aborting")
        return
    if absPath2 !="":
        data2 = pd.read_hdf(absPath2,cell2)
    else:
        print("no path given. aborting")
        return
    sns.set_style("white")
    data1=data1.transpose()
    tVector1 = data1["timeVector"].dropna()
    labels1 = get_labels(label="trial",dataFrame=data1)
    trials1 = np.array(data1[labels1].dropna())
    trials1 = normalize(trials1)
    trials1 = np.median(trials1,1)
    
    data2=data2.transpose()
    tVector2 = data2["timeVector"].dropna()
    labels2 = get_labels(label="trial",dataFrame=data2)
    trials2 = np.array(data2[labels2].dropna())
    trials2 = normalize(trials2)
    trials2 = np.median(trials2,1)
    
    if fh==None:
        fh=plt.figure()
    
    if sbHandle==None:
        sh=fh.add_subplot(211)
    else:
        sh = sbHandle
    if sb1Handle==None:
        sh2=fh.add_subplot(212)
        
    sh.axes.plot(tVector1,trials1,"k")
    
    
    sh.axes.plot(tVector2,trials2,"r")
    sh.set_ylabel("norm. trace")
    
    if clusterFlag==1:
        #trace of the first cluster
        clu1Label = int(data1["clusIndx"]["chirp"][0])
        clu2Label = int(data2["clusIndx"]["chirp"][0])
        
        clu1 = data1["c"+str(clu1Label)].dropna()
        clu2 = data2["c"+str(clu2Label)].dropna() 
        
        sh2.axes.plot(tVector1[0:len(clu1)],clu1,"k") 
        sh2.axes.plot(tVector2[0:len(clu2)],clu2,"r")
        sh.legend([cell1,cell1+"_c"+str(clu1Label)])      
        sh2.legend([cell2,cell2+"_c"+str(clu2Label)])
        return fh,sh,sh2
    else:
        sh.legend([cell1])
        sh2.legend([cell2])
        return fh,sh
    
def plot_cluster(dataFrame = None,fh=None,sh=None):
    sns.set_style("whitegrid")    
    if fh is None:
        fh=plt.figure()
    if sh is None:
        sh=fh.add_subplot(111)
        
    sh.axes.plot()
    return

def plot_cluster_one_cell(dataFrame = None,fh=None,sh=None):
    """take the dataFrame for the chirp stimulus and plot the 
    cell correlation to the N first best clusters - N is set by the
    'n_max_correlations' function"""
    sns.set_style("whitegrid")
    
    if fh is None:
        fh=plt.figure()
    if sh is None:
        sh=fh.add_subplot(111)
        
    #tickNum = len(dictionary["obs"]["corr_at_0"])
    
    sh.axes.plot(dataFrame["clusCorrs"].dropna(),color=(0.85,0.85,0.85))
    sh.axes.plot(dataFrame["clusCorrs"].dropna(),"bo")
    
    sh.axes.xaxis.set_ticks(range(dataFrame["clusCorrs"].dropna()))
    sh.axes.xaxis.set_ticklabels(dataFrame["clusIndx"].dropna())
    sh.margins(0.01)
    sh.axes.yaxis.set_view_interval(0.0,1.0)
#    sh.axes.xaxis.set_view_interval(0.0,5)   
 
    sh.set_xlabel("clusters")
    sh.set_ylabel("correlation")

    return fh,sh
  
    
def plot_trials_n_median(dataFrame = None,plotStim=1,traceColor="r",fh=None,sh=None,sh1=None):
    """plot all trials as gray traces and median as tracecolor. Input is the data frame created
    by 'raw2panda'. If handles are give, the plot will be done in a specific figure and 
    subplot. 'plotstim' flag plots the respective stimulus trace.
    returns:
        fh: figure handle.
        sh: subplot handle.
        sh1: stimulus subplot handle (if it was plotted)."""
    
    sns.set_style("white")    
    if fh is None:
        fh = plt.figure()
    if sh is None:
        sh = fh.add_subplot(211)
#    tAxis=dataFrame["timeVector"].dropna()
    labels = get_labels(label="trial",dataFrame=dataFrame)
    trials = np.array(dataFrame[labels].dropna())
    trials = normalize(trials)
    #trials = trials.transpose()
    if plotStim == 1:
#        tStim = dataFrame["stimVector"].dropna()
        if sh1 is None:
            sh1 = fh.add_subplot(211+1)
        try:
            sh1.axes.plot(dataFrame["stimVector"].dropna(),dataFrame["stimTrace"].dropna(),color="k")
        except KeyError:
            print("processed data incomplete")
            
    sh.axes.plot(dataFrame["timeVector"].dropna(),trials,color=(0.85,0.85,0.85))
    #sh.axes.plot(dataFrame["timeVector"],dataFrame[labels],color=(0.85,0.85,0.85))
    
    sh.axes.plot(dataFrame["timeVector"].dropna(),np.median(trials,axis=1),color=traceColor)



        
    return fh,sh,sh1

def get_files_to_load(folder,string):
    if folder is not None:
        contentList = os.listdir(folder)
        files2load = list()
        for item in contentList:
            if string in item:
                files2load.append(item)
    
        return files2load
    else:
        print("folder path invalid")
        return
        
def play_noise(noise,fh=None,ax=None):
    if fh == None:
        fh = plt.figure()
    if ax == None:
        ax = plt.gca()
    
    maxVal = np.max(noise)
    minVal = np.min(noise)
    
    plt.imshow(X=noise[0],vmin=minVal,vmax=maxVal)
    plt.colorbar(boundaries=[minVal,maxVal])
    plt.pause(0.5)
    
    for i in range(1,len(noise)):
       
        plt.imshow(X=noise[i],vmin=minVal,vmax=maxVal)
        
        plt.pause(0.5)
        plt.suptitle(str(i))
    return fh,ax

def plot_resp_all_stim(dataPath=None,fileid="",cell="cell1",fh=None):
    """function to plot a cell response to all stimuli used in the experiment"""
    if fh==None:
        fh = plt.figure()

    files2load = get_files_to_load(dataPath,fileid)
    print(files2load)
    for index,item in enumerate(files2load):
        hdfStore = pd.HDFStore(dataPath+item)
        keys = hdfStore.keys()
        
        sh=fh.add_subplot(3,3,index+1)
        if "/"+cell in keys:
            ind = keys.index("/"+cell)
       
            data = pd.read_hdf(dataPath+item,keys[ind])
            data = data.transpose()
            print(data["timeVector"].columns[0])
            if "noise" not in data["timeVector"].columns[0]:
                plot_trials_n_median(dataFrame=data,plotStim=1,sh=sh,sh1=sh,fh=fh)
                sh.set_title(data["timeVector"].columns[0])
#            else:
#                plt.imshow()
            if index==0:
                sh.set_xlabel("time")
                sh.set_ylabel("normalized trace")
        hdfStore.close()
            
    return fh
    
def raw2panda(rawTrace,traceTime, triggerTime,trigMode,
              sampRate=7.8,stimName="stim",trialFlag=1):
    """function to convert raw data into panda table.
    it requires the raw ROI trace, the trace time vector, the trigger time
    array, the trigger mode(1 for taking every trigger, 2 for skipping every
    other, and so on..), and sampling rate."""
    #rawTrace = roiDict["ROI"+indx]
    rawTrace = np.array(rawTrace)
    
    baseline = grab_baseline(trace=rawTrace,
                             traceTime=traceTime,
                             triggerTime=triggerTime)
                                                                                                                
    resTrace = subtract_baseline(baseline=baseline, 
                                     trace=rawTrace)
    
    if trialFlag==1:
        resMatrix,triggerTrace = trace2trial(trace=resTrace,
                                traceTime=traceTime,
                                triggerTime=triggerTime,
                                triggerMode=trigMode) 
                
        trialLabel = ['trial{0}'.format(i) for i in range(1,np.shape(resMatrix)[1]+1)]
        stimName = [stimName]*len(trialLabel)    
        array = [trialLabel,stimName]
        ind = list(zip(*array))
    
        ind1 = pd.MultiIndex.from_tuples(ind, names=['trace', 'stimType'])
           
        allData = pd.DataFrame(np.transpose(resMatrix),
                           index = ind1)
                
        medianTrace = np.median(resMatrix,axis=1)
        medianTrace = normalize(medianTrace)
                
        medianTrace = pd.Series(medianTrace,name=("medianTrace",stimName[0]))
                
        allData = allData.append(medianTrace)
               
        taxis = np.linspace(0,int(len(resMatrix)/sampRate),len(resMatrix))
#        taxis = pd.Series(taxis,name=())
        
        qi,minIdx = response_quality_index(stimMatrix=resMatrix)
        allData = allData.append(pd.Series(qi,name = ("qualIndex",stimName[0])))
        allData = allData.append(pd.Series(minIdx,name = ("minIndex",stimName[0])))
    else:
        stimName=[stimName,stimName]
#        resTrace = pd.Series(resTrace,name=("normTrace",stimName))
#        ind = ["normTrace",stimName]
#        ind = list(zip(*ind))
#                
#        ind = pd.MultiIndex.from_tuples(("trial","next"))   
        resTrace=pd.Series(resTrace,name=("normTrace",stimName[0]))
        allData = pd.DataFrame(resTrace)
        allData = allData.transpose()
        
        triggerTrace,triggerInd = create_trigger_trace(rawTrace,traceTime,triggerTime)
        triggerInd = pd.Series(triggerInd,name=("trigInd",stimName[0]))
        allData = allData.append(triggerInd)
        taxis = np.linspace(0,int(len(rawTrace)/sampRate),len(rawTrace))
        
        
        
    taxis = pd.Series(taxis,name=("timeVector",stimName[0])) 
    allData=allData.append(taxis)
    
    triggerTrace = pd.Series(triggerTrace,name=("triggerTrace",stimName[0]))
    allData = allData.append(triggerTrace)

    return allData
   
def response_quality_index(stimMatrix=None):
    """function to calculate quality index of a given cell under a certain
    stimulation. 2D Matrix has to be arranged by (timeXstimulation repetitions)
    The calculation is the variance of the mean divided by the mean of the
    variance"""
    import numpy as np
    varMean = np.var(np.mean(stimMatrix,axis=1))
    meanVar = np.mean(np.var(stimMatrix,axis=0))
    index = varMean/meanVar
    minIndex = 1./np.size(stimMatrix,axis=1)
    return index,minIndex

def reconstruct_noise(noiseList=None,greenInd=2,blueInd=3,whiteInd=1):
    """use the file used to contruct the dense noise in 
    order to reconstruct the stimulus sequence in python"""
    xPix = int(noiseList[0][0]) 
    yPix = int(noiseList[0][1])
    frames = int(noiseList[0][2])
    #noiseMat = np.zeros(shape = (frames,xPix,yPix),dtype=int)
    #greenChan = np.zeros(shape = (frames,xPix,yPix),dtype=int)
    #blueChan = np.zeros(shape = (frames,xPix,yPix),dtype=int)
    
    
    noiseMat=np.zeros(shape=(frames,xPix,yPix,3))
    i=0
    for line in noiseList:
        if len(line)>10:  
            #print(len(line))
            line=list(map(int,line))
            square=np.reshape(line,(xPix,yPix),order="C")

         
            bSquare = (square==blueInd)
            #bSquare = np.array(bSquare,dtype=float)
            gSquare = (square==greenInd)
            wSquare = (square==whiteInd)
            bSquare = bSquare+wSquare
            gSquare = gSquare+wSquare
            
            noiseMat[i,:,:,1] = np.array(gSquare,dtype=float)
            noiseMat[i,:,:,2] = np.array(bSquare,dtype=float)
            #greenChan[i,:,:] = np.array(gSquare,dtype=float)
            #blueChan[i,:,:] = np.array(bSquare,dtype=float)

            i=i+1

    
        
    return noiseMat
    
        
def read_stimulus(filePath=None):
    """read files, line by line and store them to a list"""
    fid = open(filePath,'r')
    stim=list()
    for line in fid.readlines():
        line=line[0:-1]
        line = line.split(",")
        stim.append(line)    
    fid.close()
    return stim

    
def retina_edges(rawPath="",fileList=None,pattern = pd.Series({"-x":"dorsal",
                                                 "+x":"ventral",
                                                 "-y":"nasal",
                                                 "+y":"temporal"})):
    """read in the header of the raw data files that were recorded at the edges of the retina,
    and at the optic disk (OD). Use this information to calculate the relative distance from
    the edges to the optic disk.
    inputs: 
        rawPath, optional input. Specificies from which directory to start the user dialogue
        pattern - Input required to establish how the retina was positioned on the chamber. It
        should contain a dictionary with the following keys: "+x","-x","+y","-y" and the correspondent
        values (as strings) defined by the user.
    returns:
        edgesOut : a panda data frame containing each edge location, i.e. temporal ventral,
        and the relative distance to the OD for all three coordinates (x,y,z). 
        ODout :  a panda series containing the absolute coordinate of the OD. these values can
        be used to calculate back what the absolute positioning of the recording was"""
    if fileList is None:
        root = tk.Tk()
        root.withdraw()
    
    
        #choose files to open header
        edges = fd.askopenfilenames(initialdir= rawPath,
                                filetypes=[('header files', '.smh')],
                                title = ["choose the files containing the edges recordings"])
        od =    fd.askopenfilenames(initialdir= rawPath,
                                filetypes=[('header files', '.smh')],
                                title = ["choose the files containing the OD recording"])
    
    else:
        od = [s for s in fileList if "od" in s]
        od = od[0]
        fileList.remove(od)
        edges = fileList[:]
        
    ODheader = rsm.read_in_header(filePath = od)
    odx = float(ODheader["XCoord_um"])
    ody = float(ODheader["YCoord_um"])
    #odz = float(ODheader["ZCoord_um"])
    
    ODout = pd.Series([odx,ody],["x","y"],name="OD")#odz],["x","y","z"])
    positions = list()
    coord = list()
    

    for item in edges:
                
        header = rsm.read_in_header(filePath = item)
        
        xtemp = float(header["XCoord_um"]) - odx

        if xtemp >= 0:
            locx = pattern["+x"]
        else:
            locx = pattern["-x"]
        
        ytemp = float(header["YCoord_um"]) - ody

        
        if ytemp >= 0:
            locy = pattern["+y"]
        else:
            locy = pattern["-y"]
        
        #ztemp = float(header["ZCoord_um"]) - odz
        
        #if ztemp >= 0:
        #    locz = "aboveOD"
        #else:
        #    locz = "belowOD"
        
        coord.append([xtemp,ytemp])#,ztemp])
        positions.append((locx,locy))#+" "+locz)
        
    #test = pd.DataFrame(data = [5,20],columns=[["nasal"],["dorsal"]],index=["x","y"])
    edgesOut = pd.DataFrame(data=coord,
                          index=tuple(positions),#,"z"],
                          columns=["x","y"])

    return edgesOut,ODout

def cluster_all_cells(dataFrame=None):
    
    pass


def STA(spkInd,triggerInd,stimMatrix,
        responseTrace,timeDelay=-1,
        gaussianFilter=np.ones(shape=(1,1))):
    
    #timeDelay=abs(timeDelay)
    sta = list()
    for spk in spkInd:
        
        #points where the spike happened after at least one trigger minus the timedelay
        dummie,_ = np.where(triggerInd[:] <= np.array(spk+timeDelay))
        
        if len(dummie) > 0:#if array is not empty
            
            #here only the last index is necessary
            dummie = dummie[-1]
            
            #grab the respective frame from the stimulus matrix,
            #and multiply it by the spike value
            matrix = stimMatrix[dummie]*responseTrace[spk]
            
            #filter matrix with gaussian (standard value doesn't affect the matrix)
            matrix=signal.convolve2d(in1=matrix,in2=gaussianFilter,
                                         mode="same",boundary='symm')
                                
            sta.append(matrix)
    
    #sta = np.sum(sta,axis=0)/len(sta)
            
    return sta

def subtract_baseline_matrix(baseline=None,resMatrix=None):
    """function to subtract the median of the baseline from a matrix of responses
    (timeXrepetitions)"""
    baseline=np.mean(baseline)
    resMatrix=resMatrix-baseline
    return resMatrix

def subtract_baseline(baseline=None,trace=None):
    """function to subtract the median of the baseline from the recorded trace"""
    baseline=np.median(baseline)
    trace=trace-baseline
    return trace

def trace2trial(trace=None,traceTime=None,triggerTime=None,triggerMode=1):
    """transform recorded traces and recorded 
    triggers into a timeXstimulus matrix"""
    #create binary array to indicate trigger points
    trigger,triggerInd = create_trigger_trace(trace,traceTime,triggerTime)
    #select which triggers to keep. Chirp stimululs has a wrong trigger every 
    #other trigger
    triggerInd=triggerInd[0::triggerMode]
    
    #create empty array to store data
    rows = len(triggerInd)
    cols = np.max(np.diff(triggerInd))
    #check if recording ended long enough after last trigger:
    last = len(trace[triggerInd[-1]:])
    if last < cols:
        rows=rows-1
    #add one last value to the triggerInd to cycle through the trace easier
    add = triggerInd[-1]+cols
    triggerInd.append(add)
    trials = np.empty(shape=(rows,cols))
    
    #run through all triggers.
    for i in range(1,rows+1):
        if triggerInd[i]-triggerInd[i-1]==cols:
            temp = trace[triggerInd[i-1]:triggerInd[i]]
        else:
            temp = trace[triggerInd[i-1]:triggerInd[i]+1]
        trials[i-1,0:len(temp)] = temp[:]
    trials=np.transpose(trials)
    return trials, trigger



def testTuningpy(dirs,counts,per):
    """
    Created on Fri Aug 19 10:23:32 2016
    testTuningpy is a translation to python for the matlab function "testTuning",
    originally developed in 2014 by Alexander Ecker and Philipp Berens

    @author: Andre M Chagas
    """
    itera = 1000
    counts1 = counts[:]
    
    c = np.shape(counts1)
    k = dirs[:]
    v = np.exp(per*np.multiply(np.complex(0,1),k))
    v = v/np.sqrt(c[0])
    q = np.abs(np.mean(counts1*v))
    #q = q[0]
    qdistr = np.zeros(shape=(itera,1))
    for j in range(itera):
        #r = np.random.randint(low=0,high=(c[0]*c[1])+1,size=(c[0],c[1]))
        #counts=[counts[z] for z in r]  

        
        np.random.shuffle(counts1)
        #counts1 = counts1.reshape(c)
        

        temp = np.abs(np.mean(counts1*v))
        qdistr[j] =temp
    p = np.mean(qdistr>q)
    
    return p,q,qdistr

  

def plot_field_location(field,
                        ax=None,fh=None,
                        #horizontal = "ventral",
                        #vertical = "temporal",
                        colour="b"):
    """function to plot the location of a recording respective to the retina.
    input:
        field:      dataFrame containing normalized coordinates (recording location/retina edge)
                    x and y should be indexes and the quadrant column (e.g. "ventral nasal")
        ax:         axis handle. necessary if more than one field is going to be plotted on the same 
                    figure
        fh:        figure/subplot handle. necessary if the user wants to add this plot to a figure 
                    containing more than one subplot
        #not used horizontal: the recording location respective to the horizontal axis - ventral or dorsal
        #not used vertical:   the recording location respective to the vertical axis - nasal or temporal
        colour:     the color to plot on the figure. useful if several recordings are plotted on
                    the same figure.
    output:
        fig: figure handle
        ax : Axis handle
    """
                    
    figxedge=2
    figyedge=2
    bufferzone=0.5
    

    if fh==None:
        fh = plt.figure()
    if ax==None:

        ax = fh.gca()
        retArea = plt.Circle((0,0),1,color='k',
                     linestyle="--",linewidth=2, fill=False)
        ax.add_artist(retArea)
        plt.hlines(y=0,xmax=figxedge-bufferzone,xmin=-figxedge+bufferzone)
        plt.vlines(x=0,ymax=figyedge-bufferzone,ymin=-figyedge+bufferzone)
        # change default range so that new circles will work
        ax.set_xlim((-figxedge, figxedge))
        ax.set_ylim((-figyedge, figyedge))

        ax.text(figxedge-1.5*bufferzone, figyedge-bufferzone, 'temp dorsal', style='italic',
                bbox={'facecolor':'white', 'alpha':0.5, 'pad':0})
        ax.text(figxedge-1.5*bufferzone, -figyedge+0.5*bufferzone, 'temp ventral', style='italic',
                bbox={'facecolor':'white', 'alpha':0.5, 'pad':0})
        ax.text(-figxedge+0.5*bufferzone, figyedge-bufferzone, 'nasal dorsal', style='italic',
                bbox={'facecolor':'white', 'alpha':0.5, 'pad':0})
        ax.text(-figxedge+0.5*bufferzone, -figyedge+0.5*bufferzone, 'nasal ventral', style='italic',
                bbox={'facecolor':'white', 'alpha':0.5, 'pad':0})

    plt.plot(field.loc["x"].dropna(),field.loc["y"].dropna(),'o',color=colour)
        
    return fh,ax

def remove_nans(array=None):
    array = array[-np.isnan(array)]
    return array


def process_bg(allData):
    sampRate = allData.transpose()["sampRate"].dropna().values[0][0]
    trials = ['trial{0}'.format(i) for i in range(1,4)]                           
    resMatrix = np.array(allData.transpose()[trials])
    notNans = ~np.isnan(resMatrix)                
    indx,indy = np.where(notNans==True)
    resMatrix = resMatrix[0:max(indx)+1,:]
    stim,tStim,_ = create_bg_stim(sampFreq=sampRate, greenFirst = 1)
    stim = pd.Series(stim.flatten(),name=("stimTrace","bg"))
    tStim = pd.Series(tStim.flatten(),name = ("stimVector","bg"))
    allData = allData.append(stim)
    allData = allData.append(tStim)
    green=resMatrix[0:int(len(resMatrix)/2),:]
    blue =resMatrix[int(len(resMatrix)/2)+1:,:]
                
                    #midPoint = 5 # 3 sec stimulation at 8Hz
                
    greenOn=get_color_resp_max_value(respMatrix=green,
                                                         interval=[0,12])
    blueOn =get_color_resp_max_value(respMatrix=blue,
                                                   interval=[0,12])
                
    gbIon = (greenOn-blueOn)/(greenOn+blueOn)
    gbIon = pd.Series(gbIon,name=("colorOnInd","bg"))
                
    greenOff=get_color_resp_max_value(respMatrix=green,
                                                   interval=[13,26])
    blueOff =get_color_resp_max_value(respMatrix=blue,
                                                   interval=[13,26])
                
    gbIoff = (greenOff-blueOff)/(greenOff+blueOff)
    gbIoff = pd.Series(gbIoff,name=("colorOffInd","bg"))
                
    allData = allData.append(gbIon)
    allData = allData.append(gbIoff)
                
    green = pd.Series(np.mean(green,axis=1),name=("avgGreen","bg"))
    allData = allData.append(green)
                    
    blue  = pd.Series(np.mean(blue,axis=1),name=("avgBlues","bg"))
    allData = allData.append(blue)
    
    return allData
    
def process_ds(allData, sufix):
    sampRate = allData.transpose()["sampRate"].dropna().values[0][0]
    stim,tStim,directions,screendur = create_ds_stim(sampFreq=sampRate)
    indices =[[0,8,16],[1,9,17],[2,10,18],[3,11,19],
              [4,12,20],[5,13,21],[6,14,22],[7,15,23]]
    trials = ['trial{0}'.format(i) for i in range(1,25)]

                           
    resMatrix = allData.transpose()[trials]
    resMatrix = resMatrix.dropna()
                
    resMatrix = np.array(resMatrix)
    dsMatrix = avg_matrix(matrix=resMatrix,grouping=indices)
                
    trials = ['avgTrial{0}'.format(i) for i in range(1,9)]
    name = [sufix]*len(trials) 
    ind = [trials,name]
    ind = list(zip(*ind))
                
    ind = pd.MultiIndex.from_tuples(ind)                
                 
    tempDS = pd.DataFrame(dsMatrix,columns=ind)
    allData = allData.append(tempDS.transpose())
                    
    #normalize matrix mean matrix:
    dsMatrix=(dsMatrix+np.max(np.abs(dsMatrix)))/np.max(np.abs(dsMatrix))
                
    #SVD analysis to determine direction and orientation selectivity
    normTrace,dsVector,tc = direction_selectivity(matrix=dsMatrix)
                
                
    #circular shifted trace
    #idx=np.argmax(x)
    #ctrace = np.roll(meanDatum)
                    
    # make indices                    
    # DS/OS indices
    #convert bar angles to radians
    dirRad=np.divide(directions,360.0)
    dirRad=np.multiply(dirRad,(2*np.pi))
    
                
    p,q,qdist = testTuningpy(dirs=dirRad, counts=dsVector, per=1) 
                
    allData = allData.append(pd.Series(p,name=("ds_stat_signif",sufix)))
    allData = allData.append(pd.Series(q,name=("projected_index",sufix)))
                
#    ind = [sufix]*len(qdist)    
#    ind = [qdist,ind]
#    ind = list(zip(*ind))
#    
#    ind = pd.MultiIndex.from_tuples(ind)
                
    allData = allData.append(pd.Series(qdist.flatten(),name=("ds_shuff_projected_dist",sufix)))
    #get the vector size on direction selectivity
    dsIndex = circ.resultant_vector_length(alpha=dirRad,w=dsVector,d=np.diff(dirRad[0:2]))
                                               
    dsIndex = dsIndex[0]
    dsIndex = pd.Series(dsIndex,name=("dirSelec",sufix))
    dirRad = pd.Series(dirRad,name=("rad",sufix))
    allData = allData.append(dirRad)
    allData = allData.append(dsIndex)
                    
                
    ## needs finishing for orientation selectivity
#    dsP = testTuning(dirRad,xx',1);
#    pref_dir = circ_mean(dir,x);                
#    os_index = circ_r(2*dir,x,2*diff(dir(1:2)));
#    os_p = testTuning(dir,xx',2);
#    pref_ori = circ_mean(2*dir,x);                
                


    ########ON OFF INDEX --> OOI ######
    onPix=dsMatrix[3:7,:]
    offPix=dsMatrix[28:32,:]
                
    onResp=onPix #first 250ms
    offResp=offPix#last 250ms

    ooi = (onResp.mean(axis=0)-offResp.mean(axis=0))/   \
            (onResp.mean(axis=0)+offResp.mean(axis=0))
                
    ooi = pd.Series(ooi.mean(),name=("ooi",sufix))
                
    allData = allData.append(ooi)
                
    #stimulus
    stim,tStim,directions,screenDur= create_ds_stim(sampFreq=sampRate)
                
    stim = pd.Series(stim.flatten(),name=("stimTrace",sufix))
    tStim = pd.Series(tStim,name=("stimVector",sufix))
                
    allData = allData.append(stim)
    allData = allData.append(tStim)
    
    return allData
    
def process_chirp(allData,natureData):
    sufix="chirp"
    sampRate = allData.transpose()["sampRate"].dropna().values[0][0]
    tempNatData  = natureData.copy()

    medianTrace = np.array(allData.transpose()["medianTrace"])
    medianTrace = medianTrace[~np.isnan(medianTrace)]
                                                  
    corr,corrClus,clusTrace = n_max_correlations(medianTrace,tempNatData,5)
    
    ind = [sufix]*len(corrClus)    
    ind = [corrClus,ind]
    ind = list(zip(*ind))
    
    ind = pd.MultiIndex.from_tuples(ind)
                
    corrClus = [w.replace('c', '') for w in corrClus]
    corrClus = list(map(int,corrClus))
                
    #corrClus
    allData = allData.append(pd.Series(corr,name = ("clusCorrs",sufix)))
    allData = allData.append(pd.Series(corrClus,name = ("clusIndx",sufix)))
    allData=allData.append(pd.DataFrame(clusTrace,index=ind))
                
    stim,tStim = create_chirp_stim(sampFreq=sampRate)
    stim = pd.Series(list(stim.flatten()),name=("stimTrace",sufix))
    tStim = pd.Series(list(tStim.flatten()),name=("stimVector",sufix))
                
    allData = allData.append(stim)
    allData = allData.append(tStim)                

    #chirp=1

    return allData
    
    
def process_field_location (iniFile,edgesFolder,fieldPath,pattern=None,fileList=None):
    
    if pattern == None:
        # Read .ini file cf: https://wiki.python.org/moin/ConfigParserExamples
        parser = ConfigParser()
        parser.read(iniFile)
    
        eye = parser["Animal"]["string_eye"]
    
        if eye == "left":
            pattern = pd.Series({"-x":"dorsal",
                         "+x":"ventral",
                         "-y":"nasal",
                         "+y":"temporal"})
        else:
            pattern = pd.Series({"-x":"dorsal",
                         "+x":"ventral",
                         "-y":"temporal",
                         "+y":"nasal"})

    
    edgesOut,odOut = retina_edges(rawPath=edgesFolder,pattern = pattern,fileList = fileList[:])
    print (edgesOut)
    
    fieldOut = field_loc(ODseries=odOut,
                         filePath = fieldPath,pattern = pattern)

    fieldLoc = fieldOut.index#fieldOut.columns[0].split(" ")

    for val in edgesOut.index:
        if fieldLoc[0] in val and fieldLoc[1] in val:
            index=val
#            ind1 = pd.MultiIndex.from_tuples(ind, names=['trace', 'stimType'])

    #index = tuple(index)
    fieldOut.loc[fieldLoc[0],"x"] = fieldOut.loc[fieldLoc[0],"x"]/edgesOut["x"][index]
    fieldOut.loc[fieldLoc[1],"y"] = fieldOut.loc[fieldLoc[1],"y"]/edgesOut["y"][index]

    return fieldOut
    
###########original MATLAB FUNCTION###############
#function [p, q, qdistr] = testTuning(dirs, counts, per)
#% Test significance of orientation tuning by permutation test.
#%   [p, q, qdistr] = testTuning(dirs, counts) computes a p-value for
#%   orientation tuning by running a permutation test on the per-th Fourier
#%   component.
#%
#%   Inputs:
#%       counts      matrix of spike counts as returned by getSpikeCounts.
#%       dirs        vector of directions (#directions x 1)
#%       per         fourier component to test (1 = direction, 2 =
#%                   orientation)
#%
#%   Outputs:
#%       p           p-value
#%       q           magnitude of second Fourier component
#%       qdistr      sampling distribution of |q| under the null hypothesis
#%
#% 2014 Alexander Ecker and Philipp Berens, University of Tübingen

#iter = 1000;
#[M, N] = size(counts);
#k = dirs(:);% / 180 * pi;
#v = exp(per* 1i * k) / sqrt(N);
#q = abs(mean(counts, 1) * v);
#qdistr = zeros(iter, 1);
#for j = 1 : iter
#    r = randperm(M * N);
#    counts(:) = counts(r);
#    qdistr(j) = abs(mean(counts, 1) * v);
#end
#p = mean(qdistr > q);

##############################################################################       




    
