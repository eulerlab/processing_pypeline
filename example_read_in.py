# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 09:08:02 2017

@author: andre
"""
#the next two lines are used so that Ipython automatically reloads updated libraries
%load_ext autoreload
%autoreload 2  


import os

os.chdir("E:\\github\\processing_pypeline\\")
import readScanM as rsm

filePath = "E:\\github\\processing_pypeline\\example_data\\"
headerName =  "rgc_ogb1_regular_scan.smh"
binary =  "rgc_ogb1_regular_scan.smp"

dicHeader = rsm.read_in_header(filePath = filePath+headerName)

