#!/usr/local/bin/python

#libs
import string
import re
import numpy as np

#read and print log files
# trace_arr = np.arange(1,52,1)
trace_arr = np.array([94,263])
#bw_arr = np.array([6000, 9000, 15000, 20000])
#bw_arr = np.array([3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 85000, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000])
bw_arr = np.array([10000])
method_arr = np.array([1,3, 14])
mname_arr = np.array(['DASH', 'ROI', 'OPTIMAL']) 
methodNum = len(method_arr)
inter_arr = np.array([24])
buff_arr = np.array([24])
video = "timelapseNY"
err = "no"
dir_loc = "video_timelapseNY/3840x1920/fixedBW/realtrace/"
#
START = 48 
END = 95
FRAME_NO = END - START
# metrics
q_deg_seg = np.arange(1000, dtype='f'); # store q_deg of every interval
q_std_seg = np.arange(1000, dtype='f'); # store q_deg of every interval
# full
f_avg_full = np.arange(5, dtype='f'); # average V-PSNR
f_std_full = np.arange(5, dtype='f'); # std.dev V-PSNR
q_deg_full = np.arange(5, dtype='f'); # average of quality degration in each interval
q_deg_2_full = np.arange(5, dtype='f'); # average of quality degration in each interval
# moving intervals only
f_avg_mov = np.arange(5, dtype='f'); # average V-PSNR
f_std_mov = np.arange(5, dtype='f'); # std.dev V-PSNR
q_deg_mov = np.arange(5, dtype='f'); # average of quality degration in each interval
q_deg_2_mov = np.arange(5, dtype='f'); # average of quality degration in each interval

for inter in np.nditer(inter_arr):
    for buff in np.nditer(buff_arr):
        for bw in np.nditer(bw_arr):
            log_metric = "metric_video_" + video + "_INTER_" + str(inter) + "_BUFF_" + str(buff) + "_S_" + str(START) + "_E_" + str(END) + ".txt"
            fmetric = open(log_metric,"w")
            fmetric.write('Trace\t')
            for method_id in range(0, methodNum):
                fmetric.write('%s\t' %(mname_arr[method_id]))
            for method_id in range(0, methodNum):
                fmetric.write('%s\t' %(mname_arr[method_id]))
            fmetric.write('\n')
            for trace in np.nditer(trace_arr):
                cnt = 0
                for method in np.nditer(method_arr):
                    #print trace, bw, method, inter, buff
                    VPSNR = np.zeros(shape=(2000, 30))
                    log = dir_loc +  "log_frame_TRACE_" + str(trace) + "_BW_" + str(bw)+ "_METHOD_" + str(method) + "_INTER_" + str(inter) + "_BUFF_" + str(buff) + "_EST_0.txt"
                    print log
                    f = open(log, "r")
                    # skip the first line
                    line = f.readline()
                    line = f.readline();
                    p = re.compile('\S+')
                    while line:
                        arr=p.findall(line)
                        #print arr[0], method, arr[2]
                        VPSNR[int(arr[0])][int(method)] = arr[2]
                        line = f.readline()
                    f.close()
                    # calcualte metrics
                    f_avg_full[cnt] = np.average(VPSNR[START:(END+1),[method]])
                    f_std_full[cnt] = np.std(VPSNR[START:(END+1),[method]])
                    cnt = cnt + 1
                fmetric.write('%d\t' % trace)
                for i in range (0, methodNum):
                    fmetric.write('%.1f\t' % (f_avg_full[i]))
                for i in range (0, methodNum):
                    fmetric.write('%.1f\t' % (f_std_full[i]))
                fmetric.write('\n')
            fmetric.close()
