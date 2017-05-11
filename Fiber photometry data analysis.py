# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 10:02:11 2017

@author: Ken
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cv2
import json
import peakutils
import pylab as pyl

#Select file path THIS SHOULD BE CHANGED FOR EACH ANALYSIS
dir1 = 'C:\Users\owner\Desktop\Fiber photometry data\MeA 170424'
os.chdir(dir1)
#select filelist and stim list
stims = ['control','ESP1','urine']
names = ['GCamp10_' + s for s in stims]
times = [[]]*len(stims)
plots = [[]]*len(stims)
peaks = [[]]*len(stims)
colorlist = ["m", "g", "b", "c", "m", "y", "k", "w"]

def save_figs(TIME,DATA,NAME,STIM_TIMING,NORM):
    if (min(DATA) < 0) &(max(DATA)<0):
        plt.ylim([1.2*min(DATA),0.1])
    elif min(DATA)<0:
        plt.ylim([1.2*min(DATA),1.2*max(DATA)])
    else:
        plt.ylim([0.8*min(DATA),1.2*max(DATA)])
    plt.title(NAME)
    plt.xlabel("Time(s)")
    if 'control_light' in NAME or 'red' in NAME:
        color = 'm'
    elif 'green' in NAME:
        color = 'g'
    else:
        color = 'g'
    if STIM_TIMING is not None:
        plt.plot(STIM_TIMING,[max(DATA)*1.05]*len(STIM_TIMING),'ko',ms = 4)
    if NORM == 1:
        plt.ylabel('dF/F')
        plt.plot(TIME,DATA,color)
        plt.savefig(NAME+"_norm.png",format = 'png', dpi=300)
        plt.close()
    elif NORM == 0:
        plt.ylabel('mV')
        plt.plot(TIME,DATA,color)
        plt.savefig(NAME+".png",format = 'png', dpi=300)
        plt.close()
    else:
        print 'define norm, if norm enter 1, else enter 0'
        plt.close()
        
#extract data from interleaved, for dat1, num = 1, for dat2, num = 0
def interleaved_extract(TIME,DATA,TTL1,TIPS,num):
    global OUT_DATA,OUT_TIME
    OUT_DATA = []
    OUT_TIME = []
    for i in range(len(TIPS)):
        if i == len(TIPS)-1:
            break
        #cut periods out 25ms
        elif (TIME[TIPS[i+1]] - TIME[TIPS[i]] >(0.024)) & (TIME[TIPS[i+1]] - TIME[TIPS[i]] <(0.027)) & (TTL1[TIPS[i]] == num):
            OUT_DATA.append(DATA[(TIPS[i]+TIPS[i+1])/2])
            OUT_TIME.append(TIME[(TIPS[i]+TIPS[i+1])/2])
        elif OUT_DATA == []:
            OUT_DATA = []
        else:
            OUT_DATA.append(OUT_DATA[-1])
            OUT_TIME.append(TIME[(TIPS[i]+TIPS[i+1])/2])
            
#gaussian filter
def smooth_data(DATA,window_len=11,window='hanning'):
    if DATA.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if DATA.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return DATA
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[2*DATA[0]-DATA[window_len-1::-1],DATA,2*DATA[-1]-DATA[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:  
        w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

class Fiberphotometry(object):
    def __init__(self,name):
        self.data = pd.read_csv(name + '.csv',index_col=12)
        if 'control_light' in name:
            #read 473nm detector
            self.dat2= np.array(self.data['AnalogIn-1'])
            #read 560nm detector
            self.dat1= np.array(self.data['AnalogIn-2'])
        else:
            self.dat1= np.array(self.data['AnalogIn-1'])
            self.dat2= np.array(self.data['AnalogIn-2'])
        self.time = np.array(self.data['Time(s)'])
        self.dat3 = np.array(self.data["AnalogIn-3"])
        self.ttl1 = np.array(self.data["TTL-1"])
        self.ttl2 = np.array(self.data["TTL-2"])
        #get doric imaging studio time stamps
        if os.path.exists(dir1 + '\\' + name + '_notes.json'):
            f = open(name + '_notes.json', 'r')
            stamp_data = json.load(f)
            keyList = stamp_data['Timestamps']
            len_data =  len(keyList)
            k = 0
            self.stamps = []
            while k < len_data:
                self.stamps.append(keyList[k]['Time'])
                k += 1

    def extract_data_continous(self,type,SMOOTH):
        #extracts data for 5min habit + 5 min free moving
        self.habit_last = np.array(self.stamps).max()
        up = self.time < self.habit_last
        last = self.time > (self.habit_last - 300)
        habit_period = np.where(up & last)
        if type == '5min_free':
            up1 = self.time > self.habit_last
            last1 = self.time < (self.habit_last + 300)
        elif type == '5min_fixed':
            up1 = self.time > self.habit_last
            last1 = self.time < max(self.time)
        else:
            return 'something wrong'
        if SMOOTH == 'smooth':
            self.dat1 = smooth_data(self.dat1,window_len=11,window='hanning')
        assay_period = np.where(up1 & last1)
        self.dat1_habit = self.dat1[habit_period]
        self.time_habit = self.time[habit_period]
        self.time_habit = self.time_habit - self.habit_last
        #get 300 sec data after time stamp
        self.dat1 = self.dat1[assay_period]
        self.time = self.time[assay_period] - self.habit_last
        self.dat3 = self.dat3[assay_period]
        self.stim_timing = self.time[peakutils.indexes(self.dat3,thres = 2/max(self.dat3), min_dist = 1000)]
        #if len(self.stim_timing) != 3:
            #print 'more stim timing than 3'
        #for debug
        '''
        plt.plot(self.time,self.dat3)
        plt.plot(self.stim_timing,self.dat3[peakutils.indexes(self.dat3,thres = 2/max(self.dat3), min_dist = 1000)],'mo')
        plt.show()
        plt.close()
        '''
        med_dat1 = np.median(self.dat1_habit)
        self.dat1_norm = ((self.dat1 - med_dat1)/med_dat1)

        
    def extract_data_interleaved(self,type,SMOOTH):
        #extracts data for 5min habit + 5 min free moving or head fixed stim
        #get the frequency
        self.tips = []
        for i in range(len(self.ttl1)):
            if i == 0:
                self.tips = self.tips
            elif self.ttl1[i] != self.ttl1[i-1]:
                self.tips.append(i)
            else:
                self.tips = self.tips
        
        interleaved_extract(self.time,self.dat1,self.ttl1,self.tips,1)
        self.dat1 = np.array(OUT_DATA)
        self.time1 = np.array(OUT_TIME)
        interleaved_extract(self.time,self.dat2,self.ttl1,self.tips,0)
        self.dat2 = np.array(OUT_DATA)
        self.time2 = np.array(OUT_TIME)
        if type == '5min_fixed':
            self.habit_last = np.array(self.stamps).max()
            up3 = self.time > self.habit_last
            up1 = self.time1 > self.habit_last
            up2 = self.time2 > self.habit_last
            last1 = self.time1 < max(self.time1)
            last2 = self.time2 < max(self.time2)
            last3 = self.time < max(self.time)
            assay_period3 = np.where(up3 & last3)
            self.dat3 = self.dat3[assay_period3]
            self.stim_timing = self.time[peakutils.indexes(self.dat3,thres = 2/max(self.dat3), min_dist = 1000)]
            if len(self.stim_timing) != 3:
                print 'more stim timing than 3'
        elif type == '5min_free':
            self.habit_last = max(self.time[peakutils.indexes(self.dat3,thres = 2/max(self.dat3), min_dist = 1000)])
            up1 = self.time1 > self.habit_last
            up2 = self.time2 > self.habit_last
            last1 = self.time1 < (self.habit_last + 300)
            last2 = self.time2 < (self.habit_last + 300)
        else:
            print 'something wrong'

        
        #smooth data by 1D gaussian filter
        if SMOOTH == 'smooth':
            self.dat1 = smooth_data(self.dat1,window_len=11,window='hanning')
            self.dat2 = smooth_data(self.dat2,window_len=11,window='hanning')        
        #calculate dF/F (GREEN and RED) as dat1/dat2_norm
       
        if len(self.dat1) < len(self.dat2):
            self.dat2 = np.delete(self.dat2,(len(self.dat1) - len(self.dat2)),0)
            self.time2 = np.delete(self.time2,(len(self.time1) - len(self.time2)),0)
        elif len(self.dat1) > len(self.dat2):
            self.dat1 = np.delete(self.dat1,(len(self.dat2) - len(self.dat1)),0)
            self.time1 = np.delete(self.time1,(len(self.time2) - len(self.time1)),0)
        
        reg = np.polyfit(self.dat2,self.dat1,1)
        self.controlfit = reg[0]*self.dat2+ reg[1]
        self.normDat = (self.dat1 - self.controlfit)/self.controlfit
        assay_period1 = np.where(up1 & last1)
        assay_period2 = np.where(up2 & last2)
        
        self.dat1_habit = self.dat1[np.where((self.time1 < self.habit_last)& (self.time1>(self.habit_last -300)))]
        self.time1_habit = self.time1[np.where((self.time1 < self.habit_last)& (self.time1>(self.habit_last -300)))]
        self.time1_habit = self.time1_habit - self.habit_last
        self.dat2_habit = self.dat2[np.where((self.time2 < self.habit_last)& (self.time2>(self.habit_last -300)))]
        self.time2_habit = self.time2[np.where((self.time2 < self.habit_last)& (self.time2>(self.habit_last -300)))]
        self.time2_habit = self.time2_habit - self.habit_last
        
        med_dat1 = np.median(self.dat1_habit)
        med_dat2 = np.median(self.dat2_habit)
        self.dat1_norm = ((self.dat1 - med_dat1)/med_dat1)
        self.dat2_norm = ((self.dat2 - med_dat2)/med_dat2)
        #get 300 sec data after time stamp
        self.dat1_assay = self.dat1[assay_period1]
        self.dat1_norm_assay = self.dat1_norm[assay_period1]        
        self.time1_assay = self.time1[assay_period1] - self.habit_last
        self.dat2_assay = self.dat2[assay_period2]
        self.dat2_norm_assay = self.dat2_norm[assay_period2]                
        self.time2_assay = self.time2[assay_period2] - self.habit_last
        self.normDat_assay = self.normDat[assay_period1]

        #for debug
       #gets data chunks 10 sec around head fixed stimulation
    def plot_stim_related_signals(self):
        self.fig = plt.figure(figsize=(5,5))
        self.peaks = []
        for i in range(len(self.stim_timing)):
            up = (self.stim_timing[i]+10)>self.time
            last = (self.stim_timing[i]-10)<self.time
            period = np.where(up & last)
            self.time_chunk = self.time[period] - self.stim_timing[i]
            dat1_chunk = self.dat1[period]            
            #get std for -10 to -5 sec of stim
            last1 = (self.stim_timing[i]-5)>self.time
            period1 = np.where(last & last1)
            dat1_med = np.median(self.dat1[period1])
            #calculate dF/F)
            dat1_chunk_norm = (dat1_chunk - dat1_med)/dat1_med
            self.peaks = np.append(self.peaks,max(dat1_chunk_norm))
            plt.plot(self.time_chunk,dat1_chunk_norm*100,'lightgray')
            if i == 0:
                self.ave = dat1_chunk_norm
            else:
                self.ave = self.ave + dat1_chunk_norm
        self.ave = self.ave / len(self.stim_timing)
        plt.ylim([-5,10])
        plt.plot(self.time_chunk,self.ave*100,'g')
        plt.title(name)
        plt.ylabel('dF/F (%)')
        plt.xlabel('Time (s)')
        plt.savefig(name+"_chunks.png",format = 'png', dpi=300)
        plt.close()
    def plot_stim_related_signals_int(self):
        self.fig = plt.figure(figsize=(5,5))
        for i in range(len(self.stim_timing)):
            up1 = (self.stim_timing[i]+10)>self.time1_assay
            last1 = (self.stim_timing[i]-10)<self.time1_assay
            up2 = (self.stim_timing[i]+10)>self.time2_assay
            last2 = (self.stim_timing[i]-10)<self.time2_assay
            period1 = np.where(up1 & last1)
            period2 = np.where(up2 & last2)
            self.time1_chunk = self.time1_assay[period1] - self.stim_timing[i]
            time2_chunk = self.time2_assay[period2] - self.stim_timing[i]
            dat1_chunk = self.dat1_assay[period1]            
            dat2_chunk = self.dat2_assay[period2]
            if len(dat1_chunk) < len(dat2_chunk):
                dat2_chunk = np.delete(dat2_chunk,(len(dat1_chunk) - len(dat2_chunk)),0)
                time2_chunk = np.delete(time2_chunk,(len(self.time1_chunk) - len(time2_chunk)),0)
            elif len(dat1_chunk) > len(dat2_chunk):
                dat1_chunk = np.delete(dat1_chunk,(len(dat2_chunk) - len(dat1_chunk)),0)
                self.time1_chunk = np.delete(self.time1_chunk,(len(time2_chunk) - len(self.time1_chunk)),0)
            #get std for -10 to -5 sec of stim
            up1_5sec = min(self.time1_chunk) + 5 > self.time1_chunk
            last1_5sec = min(self.time1_chunk) < self.time1_chunk
            period1_5sec = np.where(up1_5sec & last1_5sec)
            #period2_5sec = np.where(last2 & last2_5sec)
            #normalize by calculating dF/F
            #dat1_med = np.median(self.dat1_assay[period1_5sec])
            #dat1_chunk_norm = (dat1_chunk - dat1_med)
            #dat2_med = np.median(self.dat2_assay[period2_5sec])
            #dat2_chunk_norm = (dat2_chunk - dat2_med)
            reg = np.polyfit(dat2_chunk,dat1_chunk,1)
            controlfit = reg[0]*dat2_chunk+ reg[1]
            normDat_chunk = (dat1_chunk - controlfit)/controlfit
            normDat_med = np.median(normDat_chunk[period1_5sec])
            normDat_chunk = (normDat_chunk - normDat_med)
            plt.plot(self.time1_chunk,normDat_chunk*100,'lightgray')
            if i == 0:
                self.ave = normDat_chunk
            else:
                self.ave = self.ave + normDat_chunk
        self.ave = self.ave / len(self.stim_timing)
        plt.ylim([-5,10])
        plt.plot(self.time1_chunk,self.ave*100,'g')
        plt.title(name)
        plt.ylabel('dF/F (%)')
        plt.xlabel('Time (s)')
        plt.savefig(name+"_chunks.png",format = 'png', dpi=300)
        plt.close()
    def plot_signals_int_free(self):
        save_figs(self.time1_assay,self.dat1_assay,name+'_green',None,0)
        save_figs(self.time2_assay,self.dat2_assay,name+'_red',None,0)
        save_figs(self.time1_assay,self.dat1_norm_assay,name+'_green',None,1)
        save_figs(self.time2_assay,self.dat2_norm_assay,name+'_red',None,1)
        save_figs(self.time1,self.normDat,name+'_norm',None,1)    
        save_figs(self.time1_habit,self.dat1_habit,name+'_green_habit',None,0) 
        save_figs(self.time2_habit,self.dat2_habit,name+'_red_habit',None,0)    
    def plot_signals(self):
        save_figs(self.time,self.dat1,name,self.stim_timing,0)
        save_figs(self.time,self.dat1_norm,name,self.stim_timing,1)
        save_figs(self.time_habit,self.dat1_habit,name+'_habit',None,0)
    def plot_signals_int(self,STIM):
        if STIM == 'fixed':
            stim = self.stim_timing
        elif STIM == 'free':
            stim = None
        save_figs(self.time1_assay,self.dat1_assay,name+'_green',stim,0)
        save_figs(self.time2_assay,self.dat2_assay,name+'_red',stim,0)
        save_figs(self.time1_assay,self.dat1_norm_assay,name+'_green_norm',stim,1)
        save_figs(self.time2_assay,self.dat2_norm_assay,name+'_red_norm',stim,1)
        save_figs(self.time1,self.normDat,name+'_norm',stim,1)
        save_figs(self.time1_assay,self.normDat_assay,name+'_norm_assay',stim,1)
        save_figs(self.time1_habit,self.dat1_habit,name+'_green_habit',None,0) 
        save_figs(self.time2_habit,self.dat2_habit,name+'_red_habit',None,0)    
    def extract_peaks(self):
        #(A) get peaks as 3 times of standard deciation of dat1_habit
        dev_dat1 = np.std(self.dat1_habit)
        med_dat1 = np.median(self.dat1_habit)
        self.time_peaks = self.time[np.where(((self.dat1 - med_dat1)> 3 * dev_dat1))]
        self.dat1_peaks = self.dat1[(self.dat1 - med_dat1)> 3* dev_dat1]
        self.normdat1_peaks = (self.dat1_peaks - med_dat1) / med_dat1
        '''#(B) get large peaks as more than 80% of (max-min)dat1 
        time_peaks = timec[np.where(((max(dat1) * 0.80 + min(dat1)*0.20)< dat1))]
        dat1_peaks = dat1[((max(dat1) * 0.8 + min(dat1)*0.20)< dat1)]
        normdat1_peaks = (dat1_peaks - med_dat1) / med_dat1'''
        '''#(C) use peakutils method
        peaks = peakutils.indexes(self.dat1,thres = 2/max(self.dat1), min_dist = 50)
        self.time_peaks = self.time[peaks]
        self.dat1_peaks = self.dat1[peaks]'''


for name,color,i in zip(names,colorlist,range(len(plots))):    
    x = Fiberphotometry(name)    
    x.extract_data_continous('5min_fixed','smooth',)
    #x.extract_data_interleaved('5min_fixed,'smooth')
    #x.plot_signals_int('free')
    x.plot_stim_related_signals()
    #x.extract_peaks()
    #peaks[i] = x.peaks
    #plt.plot(x.time_chunk,x.ave*100,color)
    #plt.close()
    #plots[i] = x.ave
    #times[i] = x.time_chunk

#collect random

'''
plt.figure(figsize=(5,5))
plt.ylim([-5,30])
plt.ylabel('dF/F (%)')
plt.xlabel('Time (s)')
for i in range(len(plots)):
    plt.plot(times[i],plots[i]*100,colorlist[i])
plt.legend(stims)
plt.savefig("170503 GCaMP11_ave.png",format = 'png', dpi=300)
'''
'''
x = Fiberphotometry('GCamp10_int_control')
x.extract_data_interleaved('5min_fixed','smooth')
'''
