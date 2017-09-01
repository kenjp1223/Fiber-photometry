# -*- coding: utf-8 -*-
"""
Created on Tue Aug 08 10:25:29 2017

@author: owner
"""


import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import cv2
import peakutils

INTERVAL = 33
ESC_KEY = 0x1b

dir1 = 'C:\Users\owner\Desktop\Fiber photometry data\\170831'
os.chdir(dir1)
name = 'SC 8_1'
video = cv2.VideoCapture(name+'.mp4')
error = 20.533

'''#MATLAB
data = pd.read_csv(name + '.csv',header = None)
time = np.array(data[0])
dat3 = np.array(data[5])
start_timing = time[peakutils.indexes(dat3,thres = 2/max(dat3), min_dist = 20000)]'''
#Doric
data = pd.read_csv(name + '.csv',header = 0,index_col = 4)
time = np.array(data['Time(s)'])
dat1 = np.array(data['AnalogIn-1'])
dat2 = np.array(data['AnalogIn-2'])
dat3 = np.array(data['AnalogIn-3'])
start_timing = time[peakutils.indexes(dat3,thres = 2/max(dat3), min_dist = 20000)]


first = [[]] * len(start_timing)
last = [[]] * len(start_timing)
habit_period = [[]] * len(start_timing)
for i in range(len(start_timing)):
    first[i] = time > start_timing[i] 
    last[i] = time < (start_timing[i] + 600 )
    habit_period[i] = np.where(first[i]&last[i])
dats1 = [dat1[habit_period[0]],dat1[habit_period[1]]]
dats2 = [dat2[habit_period[0]],dat2[habit_period[1]]]


#Process for gettin ROI

from pylab import *
from skimage import data
from skimage.viewer.canvastools import RectangleTool
from skimage.viewer import ImageViewer


def get_rect_coord(extents):
    global viewer,coord_list
    coord_list.append(extents)

def get_ROI(im):
    global viewer,coord_list

    selecting=True
    while selecting:
        viewer = ImageViewer(im)
        coord_list = []
        rect_tool = RectangleTool(viewer, on_enter=get_rect_coord) 
        print "Draw your selections, press ENTER to validate one and close the window when you are finished"
        viewer.show()
        finished=raw_input('Is the selection correct? [y]/n: ')
        if finished!='n':
            selecting=False
    return coord_list

def crop_video(VIDEO,START,NAME,ERROR):
    count = 0
    #video1 = cv2.VideoCapture('GCamp 18_1.mp4')
    #end_flag1,c_frame1 = video1.read()
    fps = VIDEO.get(cv2.CAP_PROP_FPS)
    end_flag,c_frame = VIDEO.read()
    r = get_ROI(c_frame)[0]
    c_frame = c_frame[int(r[2]):int(r[3]), int(r[0]):int(r[1])]  
    #c_frame1 = c_frame1[int(r[2]):int(r[3]), int(r[0]):int(r[1])]  
    #g_frame_zero = cv2.cvtColor(c_frame1, cv2.COLOR_BGR2GRAY)
    #backtorgb_zero = cv2.cvtColor(g_frame_zero,cv2.COLOR_GRAY2RGB)     
    h, w, channels = c_frame.shape
    rec_h = cv2.VideoWriter(NAME+'_out_habit.mp4',cv2.VideoWriter_fourcc(*'MP4V'),fps, (w, h))
    rec_s = cv2.VideoWriter(NAME+'_out_stim.mp4',cv2.VideoWriter_fourcc(*'MP4V'),fps, (w, h))
    rlist = np.array([round(i,3) for i in START])
    g_frame_zero = cv2.cvtColor(c_frame, cv2.COLOR_BGR2GRAY)
    backtorgb_zero = cv2.cvtColor(g_frame_zero,cv2.COLOR_GRAY2RGB)
    #g_frame_zero = cv2.medianBlur(g_frame_zero,21)
    time_frame = count/fps
    for i in rlist:
        x = [0]
        y = [0]
        Time = [0]
        Track = np.zeros_like(c_frame)
        shamTrack = np.zeros_like(c_frame)
        print i - (START[0]- ERROR)
        while time_frame < i - (START[0]- ERROR):
            count = count + 1
            time_frame = count/fps
            end_flag, c_frame = VIDEO.read()
            c_frame = c_frame[int(r[2]):int(r[3]), int(r[0]):int(r[1])]            
        while time_frame < 600 + i - (START[0]- ERROR):
            # グレースケール変換
            g_frame = cv2.cvtColor(c_frame, cv2.COLOR_BGR2GRAY)
            # 二値化
            #ret,th1 = cv2.threshold(g_frame,20,255,cv2.THRESH_BINARY)
            #g_frame = cv2.medianBlur(g_frame,21)
            img_diff = 255 - cv2.absdiff(g_frame,g_frame_zero)
            ret,th1 = cv2.threshold(img_diff,185,255,cv2.THRESH_BINARY)
            th1 = cv2.medianBlur(th1,21)            
            backtorgb = cv2.cvtColor(th1,cv2.COLOR_GRAY2RGB)
            backtorgb_g = cv2.cvtColor(g_frame,cv2.COLOR_GRAY2RGB)
            image, contours, hierarchy = cv2.findContours(255-th1,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
            max_area = 0
            j = 0
            if len(contours) != 0:           
                for j in range(len(contours)):
                    area = cv2.contourArea(contours[j])
                    if max_area < area:
                        max_area = area
                        max_area_contours = j
                cnt = contours[max_area_contours]
                M = cv2.moments(cnt)
                if M['m00'] == 0:
                    x_frame = x[-1]
                    y_frame = y[-1]
                else:
                    x_frame = int(M['m10']/M['m00'])
                    y_frame = int(M['m01']/M['m00'])
            else:
                x_frame = x[-1]
                y_frame = y[-1]
            
            #if video taken in red light, use backtorgb instead of c_frame
            frame = cv2.drawContours(backtorgb_g, contours, -1, (0,255,0), 3)
            frame = cv2.circle(frame,(x_frame,y_frame), 5, (255,255,255), -1)
            #for debug
            cv2.imshow(NAME + ' contour',frame)
            if len(x) == 1:
                Track = cv2.line(Track, (x_frame,y_frame), (x_frame,y_frame),(0,153,244), 1,cv2.LINE_AA)
                shamTrack = cv2.line(shamTrack, (x_frame,y_frame), (x_frame,y_frame),(0,153,244), 1,cv2.LINE_AA)
            else:
                Track = cv2.line(Track, (x_frame,y_frame), (x[-1],y[-1]),(0,153,244), 1,cv2.LINE_AA)
                shamTrack = cv2.line(shamTrack, (x_frame,y_frame), (x[-1],y[-1]),(0,153,244), 1,cv2.LINE_AA)  
            x = np.append(x, x_frame)
            y = np.append(y, y_frame)
            Time = np.append(Time,time_frame)
            #cv2.imshow("track of " + ORG_FILE_NAME, imgTrack)
            imgTrack = cv2.addWeighted(backtorgb_g,0.5,Track,1,0.)
            #for debug
            cv2.imshow(NAME,img_diff)
            cv2.imshow(NAME + 'track', imgTrack)
            count = count + 1
            time_frame = count/fps

            # フレーム書き込み
            if i == rlist[0]:
                rec_h.write(imgTrack)
            elif i == rlist[1]:
                rec_s.write(imgTrack)
            # Escキーで終了
            key = cv2.waitKey(INTERVAL)
            if key == ESC_KEY:
                break
            # 次のフレーム読み込み
            end_flag, c_frame = VIDEO.read()
            c_frame = c_frame[int(r[2]):int(r[3]), int(r[0]):int(r[1])]    
        while len(r) < len (x[1:]):
            r = np.append(r,None)
        df = pd.DataFrame({'X':x[1:],'Y':y[1:],'Time (s)':Time[1:],'ROI':np.array(r)},columns = ['Time (s)','X','Y','ROI'])
        if i == rlist[0]:
            trial = '_habit'
        elif i == rlist[1]:
            trial = '_stim'
        df.to_csv(NAME + trial + "_track result.csv")
        cv2.destroyAllWindows()
        Trace = cv2.addWeighted(backtorgb_zero,0.5,shamTrack,0.8,0.)
        cv2.imwrite(NAME + trial +' _track.png',Trace)

        # 終了処理
        #cv2.imwrite(NAME + ' _track.png',255 - crop_Track)
        #cv2.imwrite(NAME + ' _track1.png',255 - crop_Track1)
    cv2.destroyAllWindows()
    VIDEO.release()
    rec_h.release()
    rec_s.release()
    
crop_video(video,start_timing,name,error)
