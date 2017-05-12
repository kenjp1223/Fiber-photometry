# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:24:41 2017

@author: owner
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cv2
import json
import peakutils
import pylab as pyl

dir1 = (xx)
os.chdir(dir1)
name = 'xx'
x = Fiberphotometry(name)

plt.ylabel('mV')
plt.xlabel('Time (s)')
plt.plot(x.time,x.dat1,'green')
plt.plot(x.time,x.dat3,'m')
plt.plot(stim_timing,[max(x.dat3)+0.5]*len(stim_timing),'.')
stim_timing = x.time[peakutils.indexes(x.dat3,thres = 1/max(x.dat3), min_dist = 1000)]
