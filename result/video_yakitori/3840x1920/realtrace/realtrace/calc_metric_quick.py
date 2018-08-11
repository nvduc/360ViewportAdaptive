#!/usr/local/bin/python

#libs
import string
import re
import numpy as np

#read and print log files
#trace_arr = np.arange(1,62,1)
trace_arr = np.array([27, 22, 28])
#bw_arr = np.array([6000, 9000, 15000, 20000])
#bw_arr = np.array([3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 85000, 9000, 9500, 10000, 10500, 11000, 11500, 12000, 12500, 13000, 13500, 14000, 14500, 15000])
bwtrace_arr = np.array([0])
method_arr = np.array([1,3,4,5,6,12,14])
#method_arr = np.array([20, 21])
mname_arr = np.array(['DASH','ROI', 'OPT-1', 'OPT-2','Ghent', 'Prob360', 'OPTIMAL']) 
#mname_arr = np.array(['OPT-1', 'OPT-2'])
methodNum = len(method_arr)
inter_arr = np.array([32])
buff_arr = np.array([32, 96])
video = "yaki_est"
err = "no"
dir_loc = "10M"
VP_EST_METHOD = 1
# moving parts marking
p1 = np.arange(224,352,1)
p2 = np.arange(576,1120,1)
p3 = np.arange(1344,1664,1)
move_part = np.concatenate((p1,p2,p3))
# segments
p1_seg = np.arange(7,11,1)
p2_seg = np.arange(18,35,1)
p3_seg = np.arange(42,52,1)
mov_part_seg = np.concatenate((p1_seg,p2_seg,p3_seg))
#
START = 64
END = 1791 
FRAME_NO = END - START
# metrics
q_deg_seg = np.arange(1000, dtype='f'); # store q_deg of every interval
q_std_seg = np.arange(1000, dtype='f'); # store q_deg of every interval
# full
f_avg_full = np.arange(100, dtype='f'); # average V-PSNR
f_std_full = np.arange(100, dtype='f'); # std.dev V-PSNR
q_deg_full = np.arange(100, dtype='f'); # average of quality degration in each interval
q_deg_2_full = np.arange(100, dtype='f'); # average of quality degration in each interval
# moving intervals only
f_avg_mov = np.arange(100, dtype='f'); # average V-PSNR
f_std_mov = np.arange(100, dtype='f'); # std.dev V-PSNR
q_deg_mov = np.arange(100, dtype='f'); # average of quality degration in each interval
q_deg_2_mov = np.arange(100, dtype='f'); # average of quality degration in each interval

for inter in np.nditer(inter_arr):
    for buff in np.nditer(buff_arr):
        for bwtrace in np.nditer(bwtrace_arr):
            log_metric = "metric_video_" + video + "_BWTRACE_" + str(bwtrace) +  "_INTER_" + str(inter) + "_BUFF_" + str(buff) + "_VP_EST_" + str(VP_EST_METHOD)+ "_S_" + str(START) + "_E_" + str(END) + ".txt"
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
                    VPSNR = np.zeros(shape=(2000, 100))
                    log =  "log_frame_TRACE_" + str(trace) + "_BWTRACE_" + str(bwtrace)+ "_METHOD_" + str(method) + "_INTER_" + str(inter) + "_BUFF_" + str(buff) + "_EST_" + str(VP_EST_METHOD) + ".txt"
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
                    f_avg_mov[cnt] = np.average(VPSNR[move_part, [method]])
                    f_std_mov[cnt] = np.std(VPSNR[move_part, [method]])
                    # quality degradataion
                    segNo = FRAME_NO / inter
                    q_deg_full[cnt] = 0
                    q_deg_2_full[cnt] = 0
                    for sid in range(0, segNo):
                        q_deg_seg[sid] = 0
                        for i in range(1, inter):
                            q_deg_seg[sid] = q_deg_seg[sid] + VPSNR[(sid+2) * inter, [method]] - VPSNR[(sid+2) * inter + i, [method]]
                        q_deg_seg[sid] = q_deg_seg[sid]/(inter-1)
                        q_std_seg[sid] = np.std(VPSNR[sid*inter:((sid+1) * inter),[method]])
                        q_deg_full[cnt] = q_deg_full[cnt] + q_deg_seg[sid]
                        q_deg_2_full[cnt] = q_deg_2_full[cnt] + q_std_seg[sid]
                    q_deg_full[cnt] = q_deg_full[cnt]*1.0/segNo
                    q_deg_2_full[cnt] = q_deg_2_full[cnt]*1.0/segNo
                    q_deg_mov[cnt] = np.average(q_deg_seg[mov_part_seg])
                    q_deg_2_mov[cnt] = np.average(q_std_seg[mov_part_seg])
                    cnt = cnt + 1
                fmetric.write('%d\t' % trace)
                for i in range (0, methodNum):
                    fmetric.write('%.1f\t' % (f_avg_full[i]))
                for i in range (0, methodNum):
                    fmetric.write('%.1f\t' % (f_std_full[i]))
                fmetric.write('\n')
            fmetric.close()
