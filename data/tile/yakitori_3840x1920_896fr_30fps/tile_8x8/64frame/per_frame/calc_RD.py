#!/usr/local/bin/python

import string
import re
import numpy as np
import math

# variables
# constant
FR_NUM = 900
FPS = 30
VER_NUM = 7
TILE_NUM = 64
TILE_W = 480
TILE_H = 240
W=3840
H=1920
# arrays
tile_arr = np.arange(TILE_NUM)
tsize = np.zeros((TILE_NUM, FR_NUM, VER_NUM))
tPSNR = np.zeros((TILE_NUM, FR_NUM, VER_NUM))
fPSNR = np.zeros((FR_NUM, VER_NUM))
fsize = np.zeros((FR_NUM, VER_NUM))
vPSNR = np.zeros((1, VER_NUM))
vsize = np.zeros((1, VER_NUM))
# read tile info
for tid in range(0, TILE_NUM):
    log = "tile_" + str(tid) + ".txt"
    f = open(log, "r")
    line = f.readline()
    p = re.compile('\S+')
    fid = 0 # frame id
    while line:
        arr = p.findall(line)
        #print arr
        for vid in range(0,VER_NUM): 
            tsize[tid][fid][vid] = arr[2*vid]
            tPSNR[tid][fid][vid] = arr[2*vid + 1]
            #print tsize[tid][fid][vid], tPSNR[tid][fid][vid]
        fid = fid + 1
        line = f.readline()
    #print log, fid
#print tsize[0,:,0]
# calculate BR and PSNR of each frame
for fid in range(0, FR_NUM):
    for vid in range(0,VER_NUM):
        fPSNR[fid][vid] = 0
        fsize[fid][vid] = 0
        for tid in range(0, TILE_NUM):
            fsize[fid][vid] = fsize[fid][vid] + tsize[tid][fid][vid]
            fPSNR[fid][vid] = fPSNR[fid][vid] + (255 * 255)/(np.power(10,tPSNR[tid][fid][vid]/10.0))
        fPSNR[fid][vid] = 10 * math.log10(255*255/(TILE_W*TILE_H*1.0/(W*H*1.0)*fPSNR[fid][vid]))
# calculate R-D of each version
print fsize[:,0]
print fPSNR[:,0]
for vid in range(0, VER_NUM):
    vsize[0][vid] = np.sum(fsize[:,vid])/1000.0/(FR_NUM * 1.0/FPS)
    vPSNR[0][vid] = np.average(fPSNR[:,vid])
    print vsize[0][vid], vPSNR[0][vid]
