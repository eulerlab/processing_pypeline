# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:03:55 2016

@author: andre
"""

#script to assign data collected to a previously existing cluster
#use the next two lines to add reload to ipython. 
%load_ext autoreload

%autoreload 2  

#import copy
import pandas as pd
import os
import numpy as np
#import h5py
import re
import matplotlib.pyplot as plt
import scipy.signal as signal
from configparser import ConfigParser
import scipy.io.matlab as matlab

os.chdir("E:\\Dropbox\\code\\repositories\\retina_analysis\\PreprocExport")
import ImportPreprocessedData as ipd

#os.chdir("E:\\Dropbox\\code\\repositories\\retina_analysis\\")
import classFuncs as cfs

os.chdir("E:\\github\\processing_pypeline\\")
import readScanM as rs

#load data from nature paper
natureFile = "Z:\\User\\Chagas\\nature2016_RGCs\\2016_clusters.mat"
natureData = matlab.loadmat(natureFile)

dbPath = "Z:\\Data\\Chagas\\"


folderName = "20170119\\"
subFolder = "1\\"



filePath = dbPath+folderName+subFolder+"Pre\\"
headerPath = dbPath+folderName+subFolder+"Raw\\"

#not to comprimise the folder/file structure of the DB,
#use another location on the server to save the data
storePath = "Z:\\User\\Chagas\\analysisResults\\"

tree = cfs.get_folder_tree(filePath)
fileList= [s for s in tree if "Pre" in s]
fileList.sort()

#flag to create the noise matrix
noiseF=1
noiseCount=0
#cNoiseF=1  

#turn analysis of responses to certain stimulus on/off.
noiseFlag=True
cnoiseFlag=True
bgFlag=True
chirpFlag=True
dsFlag=True
darkdsFlag=True
spotFlag=True
flickerFlag=True

coorFlag=1

for filePrefix in fileList:
    print(filePrefix)
    splitFile = filePrefix.split("\\")

    fileName = splitFile[-1]
    #get stimulus used from filename:
    #reverse file name and find the first "_"
    index = fileName[::-1].index("_")
    #now subtract the full lenght of the string to get the index
    # in the non-reversed string
    index = len(fileName)-index
    sufix = fileName[index:fileName.index(".")]

#    if "noise" in filePrefix or "cnoise" in filePrefix:
#        storeName2 = filePrefix[:]
#        storeName2.replace("Data","user")
#        storeName2.replace(folderName,"analysisResults\\"+folderName)
#        storeName2.replace("\\Pre","")
#        storeName2.replace(".h5","_resp.h5")
#        hdStore2 = pd.HDFStore(storeName2,"w") 
        
    
    #substitute the "Pre" folder with "Raw" folder
    rawFile = filePrefix.replace("\\Pre","\\Raw")
    #replace "h5" with "smh"
    rawFile = rawFile.replace(".h5",".smh")
    #grab the string up to the raw folder
    rawFolder = rawFile[0:rawFile.index("Raw")+4]
    

    #get ini file, in order to get which eye was used
    experimentIni = filePrefix[0:filePrefix.index("\\Pre")+1]
    experimentIni = cfs.get_folder_tree(experimentIni)
    experimentIni = experimentIni[".ini" in experimentIni]

    #get all raw recordings filenames in experiment folder
    recList = cfs.get_folder_tree(rawFolder)
    #keep only the header files
    recList = [s for s in recList if ".smh" in s]
    #keep only the files related to one field
    recList = [s for s in recList if fileName[0:2] in s]
    if coorFlag==1:
        #grab only the ones related to coordinate recording (edges + optic disk)
        coorList = [s for s in recList if fileName[0:index] not in s]
        coorFlag = 0
        

    
    if len(coorList) == 5:
        fieldOut = cfs.process_field_location(iniFile = experimentIni,
                                          edgesFolder = rawFolder,
                                          fieldPath = rawFile ,pattern=None,
                                          fileList = coorList[:]) 
        fieldx = fieldOut["x"].dropna().values[0]
        fieldy = fieldOut["y"].dropna().values[0]
        
        fieldx=pd.Series(fieldx,index=["0"],name=("x "+fieldOut.index[0],sufix))
        fieldy=pd.Series(fieldy,index=["1"],name=("y "+fieldOut.index[1],sufix))
        
    else:
        print("something wrong with coordinate file headers")
            
    imported, data = ipd.importPreprocessedData(filePath,fileName)
    #close file
    imported.close()
                                        
    #get all traces
    allTraces = data["Traces0_raw"]
    allTraces = allTraces.transpose()
    dim=1
        
    #get OS_parameter dict
    osPar = data["OS_Parameters"].to_dict()
    #get sampling frequency
    sampRate = osPar["Value"]["'samp_rate_Hz'"]
            
    #get stimulator delay
    stimDel = osPar["Value"]["'StimulatorDelay'"]
    #get all rois
    rois = data["ROIs"]
    dicHead = rs.read_in_header(filePath = rawFile)

    xcoor = cfs.coor2series(dicHead,"XCoord_um",sufix=sufix)
    ycoor = cfs.coor2series(dicHead,"YCoord_um",sufix=sufix)
    zcoor = cfs.coor2series(dicHead,"ZCoord_um",sufix=sufix)
    
    storeName = filePrefix[:]
    storeName = storeName.replace("Data","user")
    storeName = storeName.replace(folderName,"analysisResults\\"+folderName)
    storeName = storeName.replace("\\Pre","")
    storeName = storeName.replace(".h5","_panda.h5")
    
    hdStore = pd.HDFStore(storeName,"w") 

#    hdStore = pd.HDFStore(storePath+folderName+subFolder+fileName[0:-3]+"_panda.h5","w")
        
    for i in range(0,np.size(allTraces,dim)):
#    for i in [33,40,55,57,30,65,25,46,62,31,52,2,27,9,24,60,59,22,64,14,44,20,17,43,53,19,21]:
        temp = "cell"+str(i+1)
        print(temp) 
        pixels,roix,roiy = cfs.calc_area(rois,i)
        pixels=pixels*((110/64)*(110/64)) 
        pixels = pd.Series(pixels,name=("area_um",sufix))
        samp = pd.Series(sampRate,name=("sampRate",sufix))
        stimDel = pd.Series(stimDel,name=("stimulator_delay",sufix))      
                
        #set the index to call from the roi dictionary
        if i+1 < 10:
            indx="00"+str(i+1)
        elif i+1 < 100:
            indx="0"+str(i+1)
        else:
            indx=str(i+1) 


        if "noise" not in  filePrefix and "cnoise" not in filePrefix:
            if "chirp" in filePrefix or "bg" in filePrefix or "spot" in filePrefix:
                trig = 2
                flag = 1
            else:
                trig = 1
                flag = 1
        else:  
            trig=1
            flag=0

        rawTrace = allTraces["ROI"+indx]
#        plt.figure()
#        plt.subplot(2,1,1)
#        plt.plot(rawTrace)
#        plt.subplot(2,1,2)
#        plt.plot(np.diff(rawTrace))
#        plt.legend(["cell"+str(i+1)])
        
        traceTime = data["Tracetimes0"][i,:]
        triggerTime = data["Triggers"]["Trigger Time"]
#        allData = 0
#        del allData
        allData = cfs.raw2panda(rawTrace=rawTrace,traceTime=traceTime, 
                                    triggerTime=triggerTime,
                                    trigMode=trig,sampRate=sampRate,
                                    stimName=sufix,trialFlag=flag)                                 

        allData = allData.append(pixels)
        allData = allData.append(samp)
        allData = allData.append(stimDel)
        allData = allData.append(xcoor)
        allData = allData.append(ycoor)
        allData = allData.append(zcoor)
        allData = allData.append(fieldx)
        allData = allData.append(fieldy)
        
        if "bg" in filePrefix and bgFlag is True:
            allData = cfs.process_bg(allData)
             
           
        if "ds" in filePrefix and dsFlag is True:
            allData = cfs.process_ds(allData,sufix)
        
        if "darkds" in filePrefix and darkdsFlag is True:
            allData = cfs.process_ds(allData,sufix)
            
        if "chirp" in filePrefix and chirpFlag is True:
            allData = cfs.process_chirp(allData,natureData)
                    
        if "flicker" in filePrefix and flickerFlag is True:

            stim,tStim = cfs.create_step(sampFreq=sampRate,sizes=1,onTime=2.0,offTime=2.0)
            stim = pd.Series(stim.flatten(),name=("stimTrace",sufix))
            tStim = pd.Series(tStim.flatten(),name = ("stimVector",sufix))
            allData = allData.append(stim)
            allData = allData.append(tStim)      
                    
        if "spot" in filePrefix and spotFlag is True:
                    
            stim,tStim = cfs.create_step(sampFreq=sampRate,sizes=2,onTime=2.0,offTime=2.0)
            stim = pd.Series(stim.flatten(),name=("stimTrace",sufix))
            tStim = pd.Series(tStim.flatten(),name = ("stimVector",sufix))
                
            allData = allData.append(stim)
            allData = allData.append(tStim)               
                
        if ("noise" in filePrefix and noiseFlag is True and "cnoise" not in filePrefix) or \
            ("cnoise" in filePrefix and cnoiseFlag is True):
            
            stimuliPath="Z:\\User\\Chagas\\RGCs_stimuli\\"
            storeName2 = filePrefix[:]
            storeName2 = storeName2.replace("Data","user")
            storeName2 = storeName2.replace(folderName,"analysisResults\\"+folderName)
            storeName2 = storeName2.replace("\\Pre","")
            storeName2 = storeName2.replace(".h5","_resp.h5")
            hdStore2 = pd.HDFStore(storeName2,"w") 
        
            if "noise" in filePrefix and "cnoise" not in filePrefix:
                noiseFile = "BWNoise_official.txt"
                    
                    #cNoise=False
            else:
                noiseFile = "colorNoise.txt"
#                        cNoise=True
#                        hdStore1 = pd.HDFStore(storePath+folderName+subFolder+filePrefix+sufix+"_cnoise.h5","w")
                
            noiseList = cfs.read_stimulus(stimuliPath+noiseFile)
                
            noise = cfs.reconstruct_noise(noiseList)
                
            

            #if noiseF==1 and noiseCount<2:
            if ~os.path.isfile(storePath+folderName+subFolder+filePrefix+sufix+"_stim.h5"):
                noiseStimPath = storeName[0:storeName.index(".")-6]
                hdStore1 = pd.HDFStore(noiseStimPath+"_stim.h5","w")
                for j in range(len(noise)):
                    rFrame = pd.DataFrame(noise[j,:,:,0]) 
                    gFrame = pd.DataFrame(noise[j,:,:,1]) 
                    bFrame = pd.DataFrame(noise[j,:,:,2]) 
                    exec("hdStore1['red_frame_"+str(j)+"'] = rFrame")
                    exec("hdStore1['green_frame_"+str(j)+"'] = gFrame")
                    exec("hdStore1['red_frame_"+str(j)+"'] = bFrame")
                        
#                        noiseCount=noiseCount+1
                hdStore1.close()
#                  noiseF=0
                    
                
                #pyplot.plot(trace)
            if len(triggerTime)>=1000:
                    
                trace,triggerInd,trigger = cfs.get_traces_n_triggers(allData)
                trace = trace[~np.isnan(trace)]
                              
                velTrace,normTrace,sd = cfs.get_vel_trace(trace)
                
                indexes = cfs.get_peaks_above_sd(trace = normTrace,sd = sd,onlypos=1)
                    
                #create 5x5 gaussian window with 1 as standard deviation
                gauss = cfs.create_2d_gaussian_window(5,5,1)
                



                                
                allTimesG = list()
                allTimesB = list()                                     
                for j in range(0,-10,-1):
                    if j<0:
                        sign="neg"
                    else:
                        sign="pos"
                        
                    rawG=cfs.STA(spkInd=indexes,triggerInd=triggerInd.dropna(),
                                       stimMatrix=noise[:,:,:,1],responseTrace=velTrace,
                                       timeDelay=j)#,gaussianFilter=gauss)
                        
#                    rawG = rawG/np.std(rawG)
                        
                    rawB=cfs.STA(spkInd=indexes,triggerInd=triggerInd.dropna(),
                                       stimMatrix=noise[:,:,:,2],responseTrace=velTrace,
                                       timeDelay=j)#,gaussianFilter=gauss)
                        
                    avgG = np.mean(rawG,axis=0)

#                    avgG = signal.convolve2d(in1=np.mean(rawG,axis=0),in2=gauss,
#                                         mode="same",boundary='symm')
                    
                    
                    avgB = np.mean(rawB,axis=0)
#                    avgB = signal.convolve2d(in1=np.mean(rawB,axis=0),in2=gauss,
#                                         mode="same",boundary='symm')
                    
                    tempG ="avg_green_RF_"+sign+str(abs(j))
                    tempB ="avg_blue_RF_"+sign+str(abs(j))

                        
                    allTimesG.append(avgG)
#                    allAvgG.append(avgG)
                    allTimesB.append(avgB)
#                    avgG=pd.DataFrame(avgG)
#                    avgB=pd.DataFrame(avgB)

                            
                    exec("hdStore2['"+temp+"_"+tempG+"'] = avgG")
                    exec("hdStore2['"+temp+"_"+tempB+"'] = avgB")
            

                                                     
                                                     
                allTimesG = (allTimesG/np.max(allTimesG))*255
                allTimesB = (allTimesB/np.max(allTimesB))*255
                for nd in range(0,9):
                
                    plt.figure()
                    temp =np.dstack((np.zeros(np.shape(allTimesG[0])),allTimesG[nd],allTimesG[nd])) 
                    plt.imshow(allTimesG[nd],
                           interpolation="None",#cmap="gray",
                           vmax=np.max(temp),vmin=np.min(temp),origin="upper")
                
                    gCon = list()
                    bCon = list()
                    
                    for i in range(len(allTimesG)):
                        gCon.append(signal.convolve2d(in1=allTimesG[i],
                                                      in2=gauss,
                                                      mode="same",boundary='symm'))
                        bCon.append(signal.convolve2d(in1=allTimesB[i],
                                                      in2=gauss,
                                                      mode="same",boundary='symm'))
                gCon = gCon-np.min(gCon)
                gCon = gCon/np.max(gCon)
                
                bCon = bCon-np.min(bCon)
                bCon = bCon/np.max(bCon)
                
#                k=1
#                plt.matshow(gCon[k],cmap="gray",vmax=1,vmin=0)
#                plt.matshow(bCon[k],cmap="gray",vmax=1,vmin=0)
                
                
                maxFrameG,maxRowG,maxColG = np.where(allTimesG==np.amax(allTimesG))
                
                idxMaxG = pd.DataFrame([maxFrameG,maxRowG,maxColG],
                                       index=[["frame","row","column"],[sufix,sufix,sufix]],
                                        columns=["max_green"])    
                
                minFrameG,minRowG,minColG = np.where(allTimesG==np.amin(allTimesG))
                
                idxMinG = pd.DataFrame([minFrameG,minRowG,minColG],
                                       index=[["frame","row","column"],[sufix,sufix,sufix]],
                                        columns=["min_green"])
                
                maxFrameB,maxRowB,maxColB = np.where(allTimesB==np.amax(allTimesB))
#                
                idxMaxB = pd.DataFrame([maxFrameB,maxRowB,maxColB],
                                       index=[["frame","row","column"],[sufix,sufix,sufix]],
                                       columns=["max_blue"])
#                
                minFrameB,minRowB,minColB = np.where(allTimesB==np.amin(allTimesB))
#                
                idxMinB = pd.DataFrame([minFrameB,minRowB,minColB],
                                       index=[["frame","row","column"],[sufix,sufix,sufix]],
                                        columns=["min_blue"])
                
#                for k in range(0,5):#len(allTimesG)):
#                    plt.matshow(test[k],
#                                vmax=np.max(test),
#                                vmin=np.min(test),cmap="gray_r")
#                    plt.matshow(testA[k],
#                                vmax=np.max(testA),
#                                vmin=np.min(testA),cmap="gray_r")
#                    if k==7:
#                        
#                        plt.colorbar()
                        
                allData = allData.append(idxMaxG)
                allData = allData.append(idxMaxB)
            #end if suifx is noise
            hdStore2.close()
            
        hdStore[temp] = allData
        del allData
#                allData.to_hdf(storePath+folderName+subFolder+filePrefix+sufix+"_panda.h5",
#                           "/"+temp+"/",append=True)

    hdStore.close()
            #del hdStore
#    if "hdStore2" in locals():
#        hdStore2.close()
##            del allData

#k=1
#plt.matshow(allTimesG[k],vmin =np.amin(allTimesG), vmax=np.amax(allTimesG))
#plt.colorbar()