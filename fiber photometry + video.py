# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 16:42:10 2017
This is a program for fiber photometry data from free moving mouse.
The experiment paradigm should be as following
    1. 5 min habituation
    2. Stim entry with TIME STAMP (press 'AnalogIn-3' button)
    3. 5 min (or more) video data

The program extracts data from photo detector (473nm) for habituation period and 5 min with stimuli.
The habituation data is used to calculate mean and standard deviation of basal Ca activity.
dF/F is calculated using mean.
    dF/F = (data(t) - mean_data) / mean_data
The peak is calculated as following
    A) (Should used) peak = (data - mean_data) over 3 times std_data
    B) peak = data over 80% of max_data

The required data are
    1. Ca data ('NAME'.csv), IF the trials used control light (560nm), name the file as 'NAME_control_light'.csv
    2. croped video (to 232 x 232 pix), ('name'.avi)

The inputs to the code is 'NAME' and the directory of file (line42 and 45)
   
Outputs the following figs, movies
    1. 5 min y = raw Ca data, x = time plot
    2. 5 min y = dF/F, x = time plot
    3. Tracks moving mouse and plots the track of the center of gravity over the input video, overlays a red circle when there is a Ca peak 

The code works as following
    1. The first half of the code extracts necessary data from the raw input data.
        returns 'dat1' as 473 nm data, 'timec' as time, 'dat1_peaks''time_peaks' as the time and 473 nm data for peaks.
    2. The second half produces a video showing the track of mouse and Ca peaks.
        This is defined as defenition "open_field_extract(video,time_peaks,name)'.

@author: owner
"""



import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cv2

#Select filie name without extension THIS SHOULD BE CHANGED FOR EACH ANALYSIS
name = 'GCamp7_snake_1'
out_name = name + '_out.avi'
#Select file path THIS SHOULD BE CHANGED FOR EACH ANALYSIS
dir1 = 'C:\Users\owner\Desktop\SF1 170411'
os.chdir(dir1)
#import data file. (index_col can be any number except 0-4)
data = pd.read_csv(name + '.csv',index_col=5)
video = cv2.VideoCapture(name+'.avi')


INTERVAL = 33
ESC_KEY = 0x1b

def save_figs(TIME,DATA,NAME,NORM):
    plt.ylim([1.5*min(DATA)-0.5*max(DATA),1.5*max(DATA)-0.5*min(DATA)])
    plt.title(NAME)
    plt.xlabel("Time(s)")
    if NORM == 1:
        plt.ylabel('dF/F')
        plt.plot(TIME,DATA)
        plt.savefig(name+"_norm.png",format = 'png', dpi=300)
    elif NORM == 0:
        plt.ylabel('mA')
        plt.plot(TIME,DATA)
        plt.savefig(name+".png",format = 'png', dpi=300)
    else:
        print 'define norm, if norm enter 1, else enter 0'
    plt.close()

#extract data
time = np.array(data['Time(s)'])
#get data for time stamp
dat3 = np.array(data["AnalogIn-3"])
dat3_c = np.array(np.where((dat3<16) & (dat3>2)))
entry_time = time[dat3_c.max()]
#get data, for control light exp (560nm) get detecor2 data as dat1
if 'control_light' in name:
    #read 473nm detector
    dat2= np.array(data['AnalogIn-1'])
    #read 560nm detector
    dat1= np.array(data['AnalogIn-2'])
else:
    dat1= np.array(data['AnalogIn-1'])
    dat2= np.array(data['AnalogIn-2'])
#get dat1 during habituation, 300 sec before time stamp
dat1_habit = dat1[np.where((time<(entry_time)) & (time > (entry_time - 300)))]
#get 300 sec data after time stamp
dat1 = dat1[np.where((time<(entry_time + 300)) & (time>  entry_time))]
dat2 = dat2[np.where((time<(entry_time + 300)) & (time>  entry_time))]
timec = time[(time < (entry_time + 300)) & (time>entry_time)] - entry_time
med_dat1 = np.median(dat1_habit)
normdat1 = ((dat1 - med_dat1)/med_dat1)
'''#(A) get peaks as 3 times of standard deciation of dat1_habit
dev_dat1 = np.std(dat1_habit)
time_peaks = timec[np.where(((dat1 - med_dat1)> 3 * dev_dat1))]
dat1_peaks = dat1[(dat1 - med_dat1)> 3* dev_dat1]
normdat1_peaks = (dat1_peaks - med_dat1) / med_dat1'''
#(B) get large peaks as more than 80% of (max-min)dat1 
time_peaks = timec[np.where(((max(dat1) * 0.80 + min(dat1)*0.20)< dat1))]
dat1_peaks = dat1[((max(dat1) * 0.8 + min(dat1)*0.20)< dat1)]
normdat1_peaks = (dat1_peaks - med_dat1) / med_dat1

#save two types of plot
save_figs(timec,dat1,name,0)
save_figs(timec,normdat1,name,1)

'''
#save plots showing peaks
plt.plot(timec,dat1)
plt.ylabel('mA')
plt.xlabel('Time(s)')
plt.plot(time_peaks,dat1_peaks,'ro')
plt.savefig(name+"_peaks.png",format = 'png', dpi=300)
plt.close()
plt.plot(timec,normdat1)
plt.ylabel('dF/F')
plt.xlabel('Time(s)')
plt.plot(time_peaks,normdat1_peaks,'ro')
plt.savefig(name+'_norm_peaks.png',format = 'png', dpi = 300)
'''

def open_field_extract(VIDEO,TIME_PEAKS,NAME):
    x = [0]
    y = [0]
    count = 0
    fps = VIDEO.get(cv2.CAP_PROP_FPS)
    end_flag,c_frame = VIDEO.read()
    h, w, channels = c_frame.shape
    Track = np.zeros_like(c_frame)
    rec = cv2.VideoWriter(NAME+'_out.avi',cv2.VideoWriter_fourcc(*'XVID'),fps, (w, h))
    shamTrack = np.zeros_like(c_frame)
    rlist = np.array([round(i,3) for i in TIME_PEAKS])
    while end_flag == True:
        # グレースケール変換
        g_frame = cv2.cvtColor(c_frame, cv2.COLOR_BGR2GRAY)
        # 二値化
        ret,th1 = cv2.threshold(g_frame,65,255,cv2.THRESH_BINARY)
        th1 = cv2.medianBlur(th1,21)
        #便宜的な背景
        ret,bg = cv2.threshold(g_frame,0,255,cv2.THRESH_BINARY)
        th1 = (bg - th1)
        image, contours, hierarchy = cv2.findContours(th1,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
        max_area = 0
        j = 0
        for j in range(len(contours)):
            area = cv2.contourArea(contours[j])
            if max_area < area:
                max_area=area;
                max_area_contours =j;
        cnt = contours[max_area_contours]
        M = cv2.moments(cnt)
        x_frame = int(M['m10']/M['m00'])
        y_frame = int(M['m01']/M['m00'])
        frame = cv2.drawContours(c_frame, contours, -1, (0,255,0), 3)
        frame = cv2.circle(frame,(x_frame,y_frame), 5, (255,255,255), -1)
        #for debug
        #cv2.imshow(NAME,frame)
        if len(x) == 1:
            Track = cv2.line(Track, (x_frame,y_frame), (x_frame,y_frame),(255, 255, 255), 1,cv2.LINE_AA)
            shamTrack = cv2.line(shamTrack, (x_frame,y_frame), (x_frame,y_frame),(255, 255, 255), 1,cv2.LINE_AA)
        else:
            Track = cv2.line(Track, (x_frame,y_frame), (x[-1],y[-1]),(255, 255, 255), 1,cv2.LINE_AA)
            shamTrack = cv2.line(shamTrack, (x_frame,y_frame), (x[-1],y[-1]),(255, 255, 255), 1,cv2.LINE_AA)
        x = np.append(x, x_frame)
        y = np.append(y, y_frame)
        #cv2.imshow("track of " + ORG_FILE_NAME, imgTrack)
        time_frame = count/fps
        if any(round(time_frame,3) == rlist):
            Track = cv2.circle(Track,(x_frame,y_frame),3,(0,0,255),-1)
            shamTrack = cv2.circle(shamTrack,(x_frame,y_frame),3,(255,255,0),-1)
        imgTrack = cv2.addWeighted(c_frame,1,Track, 1,0.)
        #for debug
        #cv2.imshow(NAME + '_track',imgTrack)
        count = count + 1
        # フレーム書き込み
        rec.write(imgTrack)
        # Escキーで終了
        key = cv2.waitKey(INTERVAL)
        if key == ESC_KEY:
            break
        # 次のフレーム読み込み
        end_flag, c_frame = video.read()
    # 終了処理
    h,w = shamTrack.shape[:2]
    crop_Track = shamTrack[0:h, ((w-h)/2):((w+h)/2)]
    #crop_Track1 = cv2.rectangle(crop_Track, (w-1, 0), ((w-w/10), h/10), (120, 0, 120), 10)
    cv2.imwrite(NAME + ' _track.png',255 - crop_Track)
    #cv2.imwrite(NAME + ' _track1.png',255 - crop_Track1)
    cv2.destroyAllWindows()
    video.release()
    rec.release()

open_field_extract(video,time_peaks,name)
