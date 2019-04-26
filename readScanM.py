# -*- coding: utf-8 -*-
"""
Created on Wed May 18 14:45:10 2016

@author: andre
"""

import numpy as np
import struct
import sys
#def get_key(key=None,dictionary=None):
#    """small function to get values from keys in the dictionary of the header.
#    necessary since the header output has a funny format. To be removed in 
#    later revisions."""
#    return dictionary[str(key)][0:-1]
    
    
def read_in_data(filePath=None, header = None, 
                 readChan1=False,readChan2=False,readChan3=False,readChan4=False):
    """function to read the binary data (the actual data coming from the
    Analog Inputs of the NI cards), as recorded by scanM. It requires the
    dictionary provided by "read_in_header" function to properly process data"""
    
    # grab some variables from header dictionary:
    frameWidth = int(header["FrameWidth"])
    frameHeight = int(header["FrameHeight"])
    
    # Stimulus Buffer per Frame (aka how many frames are stored in one chunk)
    StimBufPerFr = int(header["StimBufPerFr"])

    #it seems that the framer counter, counts backwards, so to get the number of 
    #frames, one needs to take the total number of frames and subtract from the counter
    nFrames = (int(header["NumberOfFrames"])-int(header["FrameCounter"]))*StimBufPerFr
    

    #recording buffer (aka how much of each channel is saved before the next one 
    #starts - necessary so that later one can sort the binary data according to the 
    #channels)
    pixBuffer = int(header["PixelBuffer_#0_Length"])
    
   

    ##########___________READ IN PIXELS___________________

    
    #open binary - right now the whole file is read, which could lead to problems for large files

#    
    #read data in x steps
    steps=4
    #open binary data 
    with open(filePath,"rb") as fid:
        #move to end of file and get the total number of bytes
        numBytes = fid.seek(0,2)
        #move back to beg of file
        fid.seek(0,0)
        #the data is composed of 16bit integers, 
        #meaning each number is represented with 2 bytes
        values=np.array([],dtype="int16")
        for i in range(steps):
            #print (fid.tell())
            #read the first fifth of the file
            length2read=int(numBytes/steps)
            temp = fid.read(length2read)
        
            numList = ["h"]*int(length2read/2) # each number is represented by two bytes
            numList = "".join(numList)
        
            values = np.concatenate((values,struct.unpack_from(numList,temp,offset=0)))
            
        

        

    #number of channels recorded is given by the data lenght divided by result 
    # of frameWidthXframeHeightXnFrames 
    #- couldn't find this information in the header
    nChannels = int(len(values)/(nFrames*frameWidth*frameHeight))
    
    
    #empty arrays to store data
    ###to do: preallocate array the size of each should be (nFrames*frameWidth*frameHeight)
    if readChan1 is True:
        data1=np.zeros(shape = ((nFrames+1)*frameWidth*frameHeight),dtype=int)
#        data1=np.array([],dtype="int32")
    if readChan2 is True:
        data2=np.zeros(shape = ((nFrames+1)*frameWidth*frameHeight),dtype=int)
#        data2=np.array([],dtype="int32")
    if readChan3 is True:
        data3=np.zeros(shape = ((nFrames+1)*frameWidth*frameHeight),dtype=int)
#        data3=np.array([],dtype="int32")
    if readChan4 is True:
        data4=np.zeros(shape = ((nFrames+1)*frameWidth*frameHeight),dtype=int)
        
    beg1=0
    beg2=0
    beg3=0
    beg4=0
    for i in range(0,len(values),nChannels*pixBuffer):
        if readChan1 is True:
            end1 = i+int(pixBuffer)
            data1[beg1:beg1+len(values[i:end1])]=values[i:end1]
            beg1=beg1+len(values[i:end1])
        
        if readChan2 is True and nChannels>=2:
            end2 = i+int(2*pixBuffer)
            chanInd2 = i+pixBuffer
            data2[beg2:beg2+len(values[chanInd2:end2])]=values[chanInd2:end2]
            beg2=beg2+len(values[chanInd2:end2])
        
        if readChan3 is True and nChannels>=3:
            end3 = i+int(3*pixBuffer)
            chanInd3 = i+2*pixBuffer
            data3[beg3:beg3+len(values[chanInd3:end3])]=values[chanInd3:end3]
            beg3=beg3+len(values[chanInd3:end3])
        if readChan4 is True and nChannels>=4:
            end4 = i+int(4*pixBuffer)
            chanInd4 = i+3*pixBuffer
            data4[beg4:beg4+len(values[chanInd4:end4])]=values[chanInd4:end4]
            beg4=beg4+len(values[chanInd4:end4])
                
    #run through data array to sort into the different channels
    #x=0
    
#    for i in range(0,len(values),nChannels*int(pixBuffer)):
#        if readChan1 is True:
#            channel1Indx = i 
#            data1 = np.concatenate((data1,values[channel1Indx:channel1Indx+int(pixBuffer)]))
#            
#        if nChannels > 1 and readChan2 is True:
#            channel2Indx = i+pixBuffer
#            data2 = np.concatenate((data2,values[channel2Indx:channel2Indx+int(pixBuffer)]))
#            
#        if nChannels > 2 and readChan3 is True:
#            channel3Indx = i+(2*pixBuffer)
#            data3 = np.concatenate((data3,values[channel3Indx:channel3Indx+int(pixBuffer)]))
#    
    
    output = dict()
    if readChan1 is True:
        output["chan1"] = data1
    if readChan2 is True and nChannels>=2:
        output["chan2"] = data2
    if readChan3 is True and nChannels>=3:
        output["chan3"] = data3
    if readChan4 is True and nChannels>=4:
        output["chan4"] = data4
    return output


def read_in_header(filePath=None):
    """function to read the header file recorded with scanM. 
    it stores the header data into a dictionary"""
    ###########_________READ IN HEADER FILE____________
        
    #open header
    fid = open(filePath,encoding="latin-1")

    dicHead =dict()
    for line in fid.readlines():
        #for some reason (couldn't open unicode) the lines contain a mixture of binary and string
        #skip every other one to get only the string values and
        #use .split(",") to separate the data type description from the data
        #temp=line.split(",")
        temp=line[1:-1:2].split(",")
        #print(temp)
        #print(temp)
        #store only the data, since in python the data type is defined in a different way
        data =  temp[1:]    
        if data:     #means if there is something stored in the "data" variable
            #now use the "=" sign to split the value description from the value itself
            dicInput=data[0].split("=")
            
            #remove first empty space, if it exists
            if dicInput[0][0] == " ":
                dicInput[0] = dicInput[0][1:]
            
            #remove last empty space, if it exists
            if dicInput[0][-1] == " ":
                dicInput[0] = dicInput[0][0:-1]
            
            try:
                dicHead[dicInput[0]]=dicInput[1][0:-1]
            except IndexError:
                print("read more than necessary")


    fid.close()
    return dicHead



def to_frame(dataArray=[],frameTotal=2,frameCounter=1,frameBuffer=1,frameHeight=512,frameWidth=652):
    """function to reshape the dataArray into frame format. Currently it only
    works with the direct scan mode (s shaped).Note that this function does not 
    cut off retrace periods.
    
    
    frameTotal is the total number of frames.\n
    frameCounter counts backwards from frameTotal \n
    frameBuffer is the number of frames stored in one chunck. \n
    frameHeight is the number of pixels in the y axis.\n
    frameWidth is the number of pixels in the x axis\n"""
    
    nFrames=(frameTotal-frameCounter)*frameBuffer
    
    c1=np.reshape(dataArray[0:nFrames*frameHeight*frameWidth],
                 (nFrames,frameHeight,frameWidth),
                 order="C")
    return c1


def trigger_detection(frameData,triggerLevel=220,triggerMode=1):
    """detect triggers from one of the recorded channels.
    in our setup normally channel3 contains trigger data"""
    #create zeros array. later "ones" will be placed to indicate trigger point.
    #the array length is equal to the number of frames on the channel used.
    
    trigArray = np.zeros(shape=(len(frameData)),dtype=int)
    indexes=list()
    for i,frame in enumerate(frameData):
    
        if np.any(frame>=triggerLevel):
            indexes.append(i)
    
    ########TODO#########
    ### figure out if killing triggers that happen next to one another is desirable
    #####################
    
    #drop triggers depending on trigger mode.

    indexes=indexes[::triggerMode]
    #populate triggerArray with ones
    trigArray[indexes]=1
    return indexes,trigArray
