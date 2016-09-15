# -*- coding: utf-8 -*-
"""
Created on Wed May 18 14:45:10 2016

@author: andre
"""

import numpy as np

#def get_key(key=None,dictionary=None):
#    """small function to get values from keys in the dictionary of the header.
#    necessary since the header output has a funny format. To be removed in 
#    later revisions."""
#    return dictionary[str(key)][0:-1]
    
    
def read_in_data(filePath=None, header = None):
    """function to read the binary data (the actual data coming from the
    Analog Inputs of the NI cards), as recorded by scanM. It requires the
    dictionary provided by "read_in_header" function. to properly process data"""
    
    # grab some variables from header dictionary:
    frameWidth = int(header["FrameWidth"])
    frameHeight = int(header["FrameHeight"])
    #it seems that the framer counter, counts backwards, so to get the number of 
    #frames, one needs to take the total number of frames and subtract from the counter
    nFrames = int(header["NumberOfFrames"])-int(header["FrameCounter"])
    


    #recording buffer (aka how much of each channel is saved before the next one 
    #starts - necessary so that later one can sort the binary data according to the 
    #channels)
    pixBuffer = int(header["PixelBuffer_#0_Length"])

    ##########___________READ IN PIXELS___________________

    #open binary
    fid = open(filePath,"rb")
    dummie1=fid.read(-1)
    fid.close()
    dummie = list(dummie1)

    #I suspect that only every other value contains information from the pmts
    #although it could also be that the data stored in "index" is also indicating 
    #something else - sign for instance?
#        index = np.array(dummie[1::2],dtype=int)
    values = np.array(dummie[0::2],dtype=int)

    #number of channels recorded is given by the data lenght divided by result 
    # of frameWidthXframeHeightXnFrames 
    #- couldn't find this information in the header
    nChannels = int(len(values)/(nFrames*frameWidth*frameHeight))

    #empty arrays to store data
    ###to do: preallocate array the size of each should be (nFrames*frameWidth*frameHeight)
    data1=np.array([],dtype="int32")
    data2=np.array([],dtype="int32") 
    data3=np.array([],dtype="int32")

    #run through data array to sort into the different channels
    for i in range(0,len(values),nChannels*int(pixBuffer)):
        channel1Indx = i  
        data1 = np.concatenate((data1,values[channel1Indx:channel1Indx+int(pixBuffer)]))
            
        if nChannels > 1:
            channel2Indx = i+pixBuffer
            data2 = np.concatenate((data2,values[channel2Indx:channel2Indx+int(pixBuffer)]))
            
        if nChannels > 2:
            channel3Indx = i+(2*pixBuffer)
            data3 = np.concatenate((data3,values[channel3Indx:channel3Indx+int(pixBuffer)]))
    
    return data1,data2,data3


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
                
            dicHead[dicInput[0]]=dicInput[1][0:-1]


    fid.close()
    return dicHead



def to_frame(dataArray=[],nFrames=1,frameHeight=512,frameWidth=652):
    """function to reshape the dataArray into frame format. Currently it only
    works with the direct scan mode (s shaped).Note that this function does not 
    cut off retrace periods.
    
    
    nFrames is the number of recorded frames.\n
    frameHeight is the number of pixels in the y axis.\n
    frameWidth is the number of pixels in the x axis\n"""

    
    c1=np.reshape(dataArray[0:nFrames*frameHeight*frameWidth],
                 (nFrames,frameHeight,frameWidth),
                    order="C")
    return c1
    