# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:08:02 2017

@author: andre
"""
#the next two lines are used so that Ipython automatically reloads updated libraries
%load_ext autoreload
%autoreload 2  

#import necessary libraries
import numpy as np
import os
import matplotlib.pyplot as plt
   
os.chdir("E:\\github\\processing_pypeline\\")
import readScanM as rsm

#data location
#filePath = "E:\\github\\processing_pypeline\\example_data\\"
filePath = "Z:\\User\\Chagas\\zzdata\\"
#header file name
#fileName =  "rgc_ogb1_regular_scan"
fileName = "M1_RRS1_1" 
#binary data file name
#binaryName =  "rgc_ogb1_regular_scan"

#grab header information
dicHeader = rsm.read_in_header(filePath = filePath+fileName+".smh")

#grab information from the header
frameN = int(dicHeader["NumberOfFrames"])
frameC = int(dicHeader["FrameCounter"])
frameB = int(dicHeader["StimBufPerFr"])
frameH = int(dicHeader["FrameHeight"])
frameW = int(dicHeader["FrameWidth"])

#sampRate = int(dicHeader[""])

#read in binary data, output is a dictionary, where each key is one channel.
#up to this point, the data is still serialized
output = rsm.read_in_data(filePath=filePath+fileName+".smp",header=dicHeader,
                          readChan1=True,readChan2=False,readChan3=False,readChan4=False)

#converting the data from serialized to frames. Only doing this for channel1
frame1 = rsm.to_frame(output["chan1"],frameTotal=frameN,frameCounter=frameC,frameBuffer=frameB,frameHeight=frameH,frameWidth=frameW)


#convert data from channel three to detect triggers
frame3 = rsm.to_frame(output["chan3"],frameTotal=frameN,frameCounter=frameC,frameBuffer=frameB,frameHeight=frameH,frameWidth=frameW)
indexes,trigArray1 = rsm.trigger_detection(frameData=frame3,triggerLevel=220,triggerMode=2)



