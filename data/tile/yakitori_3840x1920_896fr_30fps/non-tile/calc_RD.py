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
TILE_W = 480
TILE_H = 240
W=3840
H=1920
# arrays
fPSNR = np.zeros((FR_NUM, VER_NUM))
fsize = np.zeros((FR_NUM, VER_NUM))
vPSNR = np.zeros((1, VER_NUM))
vsize = np.zeros((1, VER_NUM))
# read tile info
log = "DASH_frame.txt"
f = open(log, "r")
line = f.readline()
p = re.compile('\S+')
fid = 0 # frame id
while line:
    arr = p.findall(line)
    #print arr
    for vid in range(0,VER_NUM): 
        fsize[fid][vid] = arr[2*vid]
        fPSNR[fid][vid] = arr[2*vid + 1]
    fid = fid + 1
    line = f.readline()
for vid in range(0, VER_NUM):
    vsize[0][vid] = np.sum(fsize[:,vid])/1000.0/(FR_NUM * 1.0/FPS)
    vPSNR[0][vid] = np.average(fPSNR[:,vid])
    print vsize[0][vid], vPSNR[0][vid]
