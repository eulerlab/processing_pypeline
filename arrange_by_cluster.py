# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 14:27:43 2017

@author: andre

"""
#next two lines are for ipython automatically reload changed classes
%load_ext autoreload
%autoreload 2 


import os
import pandas as pd
import matplotlib.pyplot as plt
import classFuncs as cfs
import numpy as np
import seaborn as sns

clusterClasses = pd.DataFrame({"c1" :"OFF local, OS",        "c2":"OFF DS",  
                               "c3" :"OFF step",             "c4":"OFF slow",             
                               "c5" :"OFF slow",             "c6":"OFF alpha sust",
                               "c7" :"OFF alpha sust",       "c8":"OFF alpha sust",
                               "c9" :"(ON-)OFF 'JAM-B' mix","c10":"OFF sust",
                               "c11":"OFF alpha trans",     "c12":"OFF alpha trans",     
                               "c13":"OFF 'mini' alpha trans",
                               "c14":"ON-OFF local-edge w3","c15":"ON-OFF local",
                               "c16":"ON-OFF local",        "c17":"ON-OFF DS 1",
                               "c18":"ON-OFF DS 1",         "c19":"ON-OFF DS 2", 
                               "c20":"(ON-)OFF local, OS",  "c21":"ON step","c22":"ON DS trans", 
                               "c23":"ON local trans, OS",  "c24":"ON local trans, OS",
                               "c25":"ON local trans, OS",  "26":"ON trans",
                               "c27":"ON trans",            "c28":"ON trans large",
                               "c29":"ON high freq",        "c30":"ON low freq",
                               "c31":"ON sust",             "c32":"ON sust",
                               "c33":"ON mini-alpha",       "c34":"ON alpha",
                               "c35":"ON DS sust 1",        "c36":"ON DS sust 2",
                               "c37":"ON slow",             "c38":"ON constrast suppr",
                               "c39":"ON constrast suppr",  "c40":"ON DS sust 3",
                               "c41":"ON local sust OS",    "c42":"OFF suppr 1",
                               "c43":"OFF suppr 1",         "c44":"OFF suppr 1",
                               "c45":"OFF suppr 1",         "c46":"OFF suppr 1",
                               "c47":"OFF suppr 2",         "c48":"OFF suppr 2",
                               "c49":"OFF suppr 2"},index=["definition"])

rootFolder ="Z:\\User\\Chagas\\analysisResults\\"

tree = cfs.get_folder_tree(rootFolder)

#get all folders that have processed data --> ".h5" tag
fieldFolders=list()
for file in tree:
    split = file.split("\\")
    if ".h5" in split[-1]:
        folder = "\\".join(split[0:-1])
        fieldFolders.append(folder)

#remove redundant info and make it a list
fieldFolders = list(set(fieldFolders))
fh=None
ax=None
#run through all folders containing recordings
for folder in fieldFolders:
    experimentTree = cfs.get_folder_tree(folder)
    experimentTree.sort()
    #copy the contents of the folder    
    treeCopy = experimentTree[:]
    #run through all files in there
    for file in experimentTree:
        fieldList=list()
        #split the complete path with "_" tag
        #this way one can separate files that came from the same field
        split = file.split("_")
        #run through the copy of the folder
        for sameField in treeCopy:
            #if the split path is the same beginning as the current file
            if split[0] in sameField:
                #append this filepath to list
                fieldList.append(sameField)
                #and remove it from original tree to avoid redundance
                experimentTree.pop(experimentTree.index(sameField))
        #from here, analyse/plot fields

    #get the filenames containing the responses to chirp and moving bar
    chirpField = [s for s in fieldList if "chirp" in s]
    dsField = [s for s in fieldList if "ds" in s and "dark" not in s]
    bgField = [s for s in fieldList if "bg" in s]

    responseFields = [chirpField[0], dsField[0]]

    #get a list of cells that had a minimal quality. In this case a response index above 
    #0.45 for chirp or 0.6 for the moving bars
    cleanListChirp,qualChirp = cfs.clean_field(filePath=responseFields[0],minQual=0.45)
    cleanListDs,qualDS  = cfs.clean_field(filePath=responseFields[1],minQual=0.6)
    cleanList = list(set(cleanListChirp + cleanListDs))
    print("useful data: " +str(len(cleanList)/len(qualChirp)*100))
    
    majorTable = pd.DataFrame()             
    for cell in cleanList:
        index=list()
        columns=list()
        toframe = list()
        
        data = pd.read_hdf(chirpField[0],cell)
        data = data.transpose()
        
        index.append(cell[1:])
        
        columns.append("area")
        toframe.append(data["area_um"].dropna().values[0][0])
        
        columns.append("chirp_qual")
        toframe.append(data["qualIndex"].values[0][0])
        
        ffi = data["qualIndex"].values[0][0]

        columns.append("clus_corr1")
        toframe.append(data["clusCorrs"].values[0][0])
            
        columns.append("cluster1")
        toframe.append(data["clusIndx"].values[0][0])
            
        columns.append("clus_corr2")
        toframe.append(data["clusCorrs"].values[1][0])
            
        columns.append("cluster2")
        toframe.append(data["clusIndx"].values[1][0])
            
        data = pd.read_hdf(dsField[0],cell)
        data = data.transpose()
        
        columns.append("ds_qual")
        toframe.append(data["qualIndex"].values[0][0])
        
        ffi = (data["qualIndex"].values[0][0]-ffi)/(data["qualIndex"].values[0][0]+ffi)    
        columns.append("full_field_index")
        toframe.append(ffi)
        
        columns.append("dir_selec")
        toframe.append(data["dirSelec"].values[0][0])
            
        columns.append("ds_p_value")
        toframe.append(data["ds_stat_signif"].values[0][0])
            
        if data["ds_stat_signif"].values[0][0] <= 0.05:
            ds=1
        else:
            ds=0
                
        columns.append("ds")
        toframe.append(ds)
        
        columns.append("on_off")
        toframe.append(data["ooi"].values[0][0])
        
        
        
        
        data = pd.read_hdf(bgField[0],cell)
        data = data.transpose()
        
        columns.append("color_on_index")
        toframe.append(data["colorOnInd"].values[0][0])
        
        columns.append("color_off_index")
        toframe.append(data["colorOffInd"].values[0][0])
        
        majorTable =majorTable.append(pd.DataFrame(data=[toframe],
                                               columns=columns,
                                               index=index))


    
    
    test = cfs.cluster_distribution(hdFilePath = responseFields[0], keys=cleanList)
    labels,numbers,cells = cfs.get_cluster_histogram(test,plot=0)    
    intLabels = [s.strip('c') for s in labels]
    intLabels = list(map(int,intLabels))
    
    data = pd.read_hdf(chirpField[0],cleanListChirp[0])
    data = data.transpose()
    loc=list()
    if "ventral" in data.columns:
        loc.append("ventral")
    else:
        loc.append("dorsal")
    if "nasal" in data.columns:
        loc.append("nasal")
    else:
        loc.append("temp")
        
    fieldx=data.loc["x"].dropna()
    fieldy=data.loc["y"].dropna()
    field = pd.DataFrame([fieldx,fieldy])
#       data[loc[0]][loc[1]].dropna()
        
    color = [np.random.random(1)[0],np.random.random(1)[0],np.random.random(1)[0]]
    print(field)
    print("\n")
        
    fh,ax= cfs.plot_field_location(field=field,fh=fh,ax=ax,colour=color)

    
#            hdfStore = pd.HDFStore(dataPath+filePrefix+dataFile)
#            keys = hdfStore.keys()
#
#            data = pd.read_hdf(dataPath+filePrefix+dataFile,keys[0])
#            data = data.transpose()
#
#            hdfStore.close()




           
        


#tree.pop(tree.index(file))

