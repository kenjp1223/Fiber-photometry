# -*- coding: utf-8 -*-
"""
Created on Wed Aug 09 14:30:47 2017

@author: owner
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cv2
import peakutils
from scipy import interpolate
import scipy.ndimage as ndi
from mpl_toolkits.axes_grid1 import make_axes_locatable

def grid_density_gaussian_filter(w, h, r,data):
    w = int(w) +1
    h = int(h) +1
    img = np.zeros((h,w))
    cal = np.zeros((h,w))
    for i in range(len(data)):
        a = data[i][0]
        b = data[i][1]
        c = data[i][2]
        img[b][a] += 1
        cal[b][a] += c
    cal = np.divide(cal,img) *100
    return img,cal, ndi.gaussian_filter(img, (r,r)), ndi.gaussian_filter(cal, (r,r)) ## gaussian convolution

dir1 = 'C:\Users\owner\Desktop\Fiber photometry data\\170816'
os.chdir(dir1)
name = 'GCamp 18_7'
trial = ['_habit','_stim']
for i in trial:
    df = pd.read_csv(name + i + '_track result.csv',header = 0,index_col = 0)
    
    r = np.array(df['ROI'])
    x = np.array(df['X'])
    y = np.array(df['Y'])
    t = np.array(df['Time (s)'])
    t = t-t[0]
    data1 = np.array(df.In1)
    data2 = np.array(df.In2)
    ta = np.zeros(len(t))
    tb = np.zeros(len(t))
    ta[np.where((x<125)&(y<125))] = 1
    tb[np.where((x>np.max(x)-125)&(y<125))] = 1
    ta[ta==0] = np.nan
    tb[tb==0] = np.nan
    a_timing = []
    b_timing = []


    for ii in range(len(t)):
        if ((x[ii]<125) & (y[ii]<125)):
            if ii ==0:
                a_timing.append(t[ii])
            elif ((x[ii-1]>=125) | (y[ii-1]>=125)):
                a_timing.append(t[ii])


    for ii in range(len(t)):
        if ((x[ii]>np.max(x)-125)&(y[ii]<125)):
            if ii ==0:
                b_timing.append(t[ii]-t[0])
            elif ((x[ii-1]<=np.max(x)-125)|(y[ii-1]>=125)):
                b_timing.append(t[ii]-t[0])

    for ii in a_timing:
        if (ta[np.where(t==ii)[0] -1] == 0)&(ta[np.where(t==ii)[0] -3] != 0):
            a_timing.remove(ii)

    for ii in b_timing:
        if (tb[np.where(t==ii)[0] -1] == 0)&(tb[np.where(t==ii)[0] -3] != 0):
            b_timing.remove(ii)
    frame = (np.max(t) - np.min(t))/(len(t)-1)
    w = r[1] - r[0] #width of ROI
    h = r[3] - r[2] #height of ROI
    data1_n = (data1 - np.median(data1))/np.median(data1)
    data2_n = (data2 - np.median(data2))/np.median(data2)
    data = zip(x,y,data2_n)
    time_A = np.count_nonzero(x <= (w/3))*frame
    time_B = np.count_nonzero((x < (w*2/3)) & (x > (int(w)/3)))*frame
    time_C = np.count_nonzero(x >= (w*2/3))*frame
    print name + i + ': Time in zone A, ' + str(np.round(time_A,1)) +'(s)' 
    print name + i + ': Time in zone B, ' + str(np.round(time_B,1)) +'(s)' 
    print name + i + ': Time in zone C, ' + str(np.round(time_C,-1)) +'(s)' 
    #Draw figure
    #H = plt.hist2d(x,y,bins = 60,range =np.array([(0,int(w)+1),(0,int(h)+1)]))
    #plt.clf()
    #heatmap = ndi.gaussian_filter(H[0],(1,1))
    img,cal,g_img,g_cal = grid_density_gaussian_filter(w,h,5,data)
    fig = plt.figure()
    
    ax1 = fig.add_subplot(111)
    
    
    #heat = (H[0].T)[-1::-1]
    ax1.set_ylim(0,int(h)+1)
    ax1.set_xlim(0,int(w)+1) 
    #ax.imshow(heatmap.T,origin = 'upper',interpolation = 'bilinear', extent = [0,int(w)+1,0,int(h)+1],vmin =0, vmax = 120)
    image1 = ax1.imshow(g_img[-1::-1],origin ='lower',vmax = 30)
    ax1.axvline(x = (w/3), ymin = 0, ymax = h-1,ls =  '--', color = 'white')
    ax1.axvline(x = (w/3)*2,ymin = 0, ymax = h-1,ls = '--', color = 'white')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    plt.title(name + i + ' density heatmap')
    #H[3].set_clim(0,150)
    image1.set_clim(0,5)
    divider = make_axes_locatable(ax1)
    ax_cb1 = divider.new_horizontal(size="4%", pad=0.05)
    fig.add_axes(ax_cb1)
    cbar1 = plt.colorbar(image1,cax = ax_cb1)
    #cbar1.ax.set_yticklabels(['0','1','2','3','<4'])  # vertically oriented colorbar
    
    #cb = fig.colorbar(H[3],ticks = [0,1/frame,2/frame,3/frame,4/frame,5/frame])
    cbar1.set_label('Time (s)')
    plt.savefig(name + i + ' density heatmap.png',format = 'png', dpi=300)
    plt.close()
    
    fig = plt.figure()
    ax2 = fig.add_subplot(111)    
    #heat = (H[0].T)[-1::-1]
    ax2.set_ylim(0,int(h)+1)
    ax2.set_xlim(0,int(w)+1) 
    #ax.imshow(heatmap.T,origin = 'upper',interpolation = 'bilinear', extent = [0,int(w)+1,0,int(h)+1],vmin =0, vmax = 120)
    image2 = ax2.imshow(cal[-1::-1],origin ='lower',cmap = 'coolwarm',vmax = 30)
    ax2.axvline(x = (w/3), ymin = 0, ymax = h-1,ls =  '--', color = 'white')
    ax2.axvline(x = (w/3)*2,ymin = 0, ymax = h-1,ls = '--', color = 'white')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    plt.title(name + i + ' GCaMP heatmap')
    #H[3].set_clim(0,150)
    image2.set_clim(-10,15)
    divider = make_axes_locatable(ax2)
    ax_cb2 = divider.new_horizontal(size="4%", pad=0.05)
    fig.add_axes(ax_cb2)
    cbar2 = plt.colorbar(image2,cax = ax_cb2, cmap = 'coolwarm',ticks=[-10,-5,0,5,10,15])
    cbar2.ax.set_yticklabels(['-10','-5','0','5','10','15'])  # vertically oriented colorbar
    
    #cb = fig.colorbar(H[3],ticks = [0,1/frame,2/frame,3/frame,4/frame,5/frame])
    cbar2.set_label('dF/F (%)')
    plt.savefig(name + i + ' GCaMP heatmap.png',format = 'png', dpi=300)
    plt.close()
    
    fig = plt.figure(figsize = (10,5))
    ax = fig.add_subplot(111)
    ax.set_ylabel('dF/F (%)')
    ax.set_xlabel('Time (s)')
    ax.set_ylim([-10,30])
    ax.plot(t,np.array(data)[:,2]*100,'g')
    ax.plot(t,(ta*max(np.array(data)[:,2])+0.01)*100,'b',label = 'Time near left cage')
    ax.plot(t,(tb*max(np.array(data)[:,2])+0.01)*100,'m',label = 'Time near right cage')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, borderaxespad=0.)
    
    #plt.subplots_adjust(left = adjl, right = adjr)
    ax.set_xlim(0,600)
    plt.savefig(name + i + ' Trace.png',format = 'png', dpi=300)
    plt.close()

'''
     #resample data for AnalogIn 1 and 2 to for behavior data
name = 'GCamp 18_7' 
data = pd.read_csv(name + '.csv',header = 0,index_col = 4)
time = np.array(data['Time(s)'])
dat1 = np.array(data['AnalogIn-1'])
dat2 = np.array(data['AnalogIn-2'])
dat3 = np.array(data['AnalogIn-3'])
start_timing = time[peakutils.indexes(dat3,thres = 2/max(dat3), min_dist = 20000)]
first = [[]] * len(start_timing)
last = [[]] * len(start_timing)
habit_period = [[]] * len(start_timing)
df = pd.read_csv(name + trial[0] + '_track result.csv',header = 0,index_col = 0)
#r = np.array(df['ROI'])
x = np.array(df['X'])
y = np.array(df['Y'])
t = np.array(df['Time (s)'])
t = t-t[0]
#gets habit data for [0] and stim data for [1]
#error = t[0]-start_timing[0]
for i in range(len(start_timing)):
    first[i] = time > start_timing[i]
    last[i] = time < (start_timing[i] + 600)
    habit_period[i] = np.where(first[i]&last[i])
times = [time[habit_period[0]],time[habit_period[1]]]
dats1 = [dat1[habit_period[0]],dat1[habit_period[1]]]
dats2 = [dat2[habit_period[0]],dat2[habit_period[1]]]

if max(t) > (times[0][-1]-times[0][0]):
    print 'time too large'
    
xnew = np.linspace(0, max(t), len(t)) # Equi-distant time vector (for FFT)
f1 = [interpolate.interp1d(times[i] - times[i][0], dats1[i]) for i in range(2)]
f2 = [interpolate.interp1d(times[i] - times[i][0], dats2[i]) for i in range(2)]
dats1_i = [f1[i](xnew) for i in range(2)]   # use interpolation function returned by `interp1d`
dats2_i = [f2[i](xnew) for i in range(2)]   # use interpolation function returned by `interp1d`

for i in range(2):
    df = pd.read_csv(name + trial[i] + '_track result.csv',header = 0,index_col = 0)
    df['In1'] = dats1_i[i]
    df['In2'] = dats2_i[i]
    df.to_csv(name + trial[i] + "_track result.csv")
'''