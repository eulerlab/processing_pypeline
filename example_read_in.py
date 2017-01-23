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
filePath = "E:\\github\\processing_pypeline\\example_data\\"
#header file name
headerName =  "rgc_ogb1_regular_scan.smh"
#binary data file name
binaryName =  "rgc_ogb1_regular_scan.smp"

#grab header information
dicHeader = rsm.read_in_header(filePath = filePath+headerName)

#grab information from the header
frameN = int(dicHeader["FrameCounter"])
frameH = int(dicHeader["FrameHeight"])
frameW = int(dicHeader["FrameWidth"])

#read in binary data, output is a dictionary, where each key is one channel.
#up to this point, the data is still serialized
output = rsm.read_in_data(filePath=filePath+binaryName,header=dicHeader,
                          readChan1=True,readChan2=True,readChan3=True)

#converting the data from serialized to frames. Only doing this for channel1
frame1 = rsm.to_frame(output["chan1"],nFrames=frameN,frameHeight=frameH,frameWidth=frameW)
