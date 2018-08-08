#include "AdaptLogic.h"
#include "Metadata.h"
#include "Projection.h"
#include "common.h"
#include <math.h>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include "common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

AdaptLogic::AdaptLogic(Metadata meta){
  metadata = meta;
  // throughput
  thrp.est_thrp_act = new double[metadata.video.NO_SEG_FULL];
  thrp.seg_thrp = new double[metadata.video.NO_SEG_FULL];
  thrp.est_seg_thrp = new double[metadata.video.NO_SEG_FULL];
  // viewport2
  vp2.cur_vp = init2dArrayInt(metadata.video.NO_SEG_FULL, 2);
  vp2.est_vp = init2dArrayInt(metadata.video.NO_SEG_FULL, 2);
  vp2.est_err = init2dArrayInt(metadata.video.NO_SEG_FULL, 2);
  vp2.est_frame_vp = init2dArrayInt(metadata.video.NO_FRAME, 2);
  // tile version
  tile_ver = init2dArrayInt(metadata.video.NO_SEG_FULL, metadata.vp.No_tile);
  tile_size = init2dArrayInt(metadata.video.NO_SEG_FULL, metadata.vp.No_tile);
  decide_width = new int[metadata.video.NO_SEG_FULL];
  //
  proc_time = new int[metadata.video.NO_SEG_FULL];
  //
  vpsnr = new double[metadata.video.NO_FRAME];
  seg_br = new double[metadata.video.NO_SEG_FULL];
  // 
  VPSNR_thres = new double[metadata.video.NO_SEG_FULL];
  avgVPSNR = new double[metadata.video.NO_SEG_FULL];
}
AdaptLogic::~AdaptLogic(){
}
int* AdaptLogic::get_next_segment(int index){
  /* calculate processing time */
  struct timeval t_start, t_end;
  int i;
  gettimeofday(&t_start, NULL);
  //
  thrp_estimator(index);
  vp_estimator(index);
  //
  switch(metadata.adapt.TILE_SELECT_METHOD){
    case 0:
      tile_ver[index] = DASH_CMP(index);
      break;
    case 1:
      tile_ver[index] = DASH_ERP(index);
      break;
    case 2:
      tile_ver[index] = EXT_ALL(index, 0); // ROI
      break;
    case 3:
      tile_ver[index] = ISM(index);
      break;
    case 4:
      tile_ver[index] = Ghent(index);
      break;
    case 5:
      tile_ver[index] = Ireland(index);
      break;
  }
  /* calculate processing time */
  gettimeofday(&t_end, NULL);
  proc_time[index] = tvdiff_us(&t_end, &t_start);
  return tile_ver[index];
}
void AdaptLogic::calc_result(){
  int index, i;
  for(index = 0; index < metadata.video.NO_SEG_FULL; index ++){
    /* calculate vpsnr and segment bitrate */
    if(metadata.adapt.TILE_SELECT_METHOD == 0){
      seg_br[index] = metadata.video.CMP_DASH_BR[index][tile_ver[index][0]];
      for(i=0; i < metadata.video.INTERVAL; i++)
        vpsnr[index * metadata.video.INTERVAL + i] = metadata.video.CMP_DASH_PSNR[index][tile_ver[index][0]];
    }else{
      if(metadata.adapt.TILE_SELECT_METHOD == 1){ // DASH
        seg_br[index] = metadata.video.DASH_BR[index][tile_ver[index][0]];
        for(i=0; i < metadata.video.INTERVAL; i++)
          vpsnr[index * metadata.video.INTERVAL + i] = metadata.video.DASH_PSNR[index][tile_ver[index][0]];
      }else{  // tiling-based
        // bitrate
        seg_br[index] = 0;
        for(i=0; i < metadata.vp.No_tile; i++)
          seg_br[index] += metadata.video.TILE_BR[index][i][tile_ver[index][i]];
        // vpsnr
        for(i=0; i < metadata.video.INTERVAL; i++)
          vpsnr[index * metadata.video.INTERVAL + i] = est_vp_psnr(metadata.video.TILE_MSE[index], metadata.vp.No_tile, tile_ver[index], metadata.trace.frame_vp[0][index * metadata.video.INTERVAL + i]);
      }
    }
  }
}
char* AdaptLogic::get_method_name(int method_id){
  char* full_name = new char[1024];
  char* name[100] = {"DASH-CMP","DASH-ERP","ROI", "ISM", "Ghent", "Ireland","test", "test2", "test3"};
  sprintf(full_name, "%df_%dx%d_%s", metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, name[method_id-1]);
  return full_name;
}
void AdaptLogic::thrp_estimator(int index){
  if(index == 0){
    thrp.est_thrp_act[index] = thrp.seg_thrp[index];
  }else{
    switch(metadata.adapt.THRP_EST_METHOD){
      case 1: // last throughput-based
        thrp.est_thrp_act[index] = (1- metadata.adapt.alpha) * thrp.est_thrp_act[index -1] + metadata.adapt.alpha * thrp.seg_thrp[index];
        break;
    }
  }
  thrp.est_seg_thrp[index] = (1-metadata.adapt.thrp_est_margin) * thrp.est_thrp_act[index];
}
void AdaptLogic::vp_estimator(int index){
  int i;
  double v_phi;
  double v_theta;
  switch(metadata.adapt.VP_EST_METHOD){
    case 1:// linear regression
      if(index < metadata.video.BUFF/metadata.video.INTERVAL){
        for(i=0; i < metadata.video.INTERVAL; i++){
          vp2.est_frame_vp[index * metadata.video.INTERVAL + i][0] = 0;
          vp2.est_frame_vp[index * metadata.video.INTERVAL + i][1] = 0;
        }
        vp2.est_err[index][0] = 0;
        vp2.est_err[index][1] = 0;
        vp2.cur_vp[index][0] = 0;
        vp2.cur_vp[index][1] = 0;	    
      }else{
        if(index == metadata.video.BUFF/metadata.video.INTERVAL){
          for(i=0; i < metadata.video.INTERVAL; i++){
            vp2.est_frame_vp[index * metadata.video.INTERVAL + i][0] = metadata.trace.frame_vp[0][0][0];
            vp2.est_frame_vp[index * metadata.video.INTERVAL + i][1] = metadata.trace.frame_vp[0][0][1];
          }
        }else{
          /* calculate velocity in phi directions (degree/frame) */
          v_phi = (metadata.trace.frame_vp[0][(index - metadata.video.BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL - 1][0] - metadata.trace.frame_vp[0][(index- metadata.video.BUFF/metadata.video.INTERVAL -1) * metadata.video.INTERVAL][0]) * 1.0 / metadata.video.INTERVAL;
          v_theta = (metadata.trace.frame_vp[0][(index - metadata.video.BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL - 1][1] - metadata.trace.frame_vp[0][(index- metadata.video.BUFF/metadata.video.INTERVAL -1) * metadata.video.INTERVAL][1]) * 1.0 / metadata.video.INTERVAL;
          /* estimate viewport position at each frame */
          for(i=0; i < metadata.video.INTERVAL; i++){
            // est_metadata.trace.frame_vp[0][index * metadata.video.INTERVAL + i][0] = metadata.trace.frame_vp[0][(index - BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL][0]  + (i + BUFF) * v_phi;
            // est_metadata.trace.frame_vp[0][index * metadata.video.INTERVAL + i][1] = metadata.trace.frame_vp[0][(index - BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL][1]  + (i + BUFF) * v_theta;
            vp2.est_frame_vp[index * metadata.video.INTERVAL + i][0] = metadata.trace.frame_vp[0][(index - metadata.video.BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL][0]  + i * v_phi;
            vp2.est_frame_vp[index * metadata.video.INTERVAL + i][1] = metadata.trace.frame_vp[0][(index - metadata.video.BUFF/metadata.video.INTERVAL)* metadata.video.INTERVAL][1]  + i * v_theta;
          }
        }
        vp2.est_err[index][0] = metadata.trace.frame_vp[0][(index-1) * metadata.video.INTERVAL][0] - vp2.est_frame_vp[(index-1) * metadata.video.INTERVAL][0];
        if(vp2.est_err[index][0] < -180)
          vp2.est_err[index][0] += 360;
        else if(vp2.est_err[index][0] > 180)
          vp2.est_err[index][0] -= 360;
        // theta
        vp2.est_err[index][1] = metadata.trace.frame_vp[0][(index-1) * metadata.video.INTERVAL][1] - vp2.est_frame_vp[(index-1) * metadata.video.INTERVAL][1];
        if(vp2.est_err[index][1] < -90)
          vp2.est_err[index][1] += 180;
        else if(vp2.est_err[index][1] > 90)
          vp2.est_err[index][1] -= 180;
        //
        vp2.cur_vp[index] = metadata.trace.frame_vp[0][(index-metadata.video.BUFF/metadata.video.INTERVAL)*metadata.video.INTERVAL];
      }
      break;
    case 2: // know all
      for(i=0; i < metadata.video.INTERVAL; i++){
        vp2.est_frame_vp[index * metadata.video.INTERVAL + i] = metadata.trace.frame_vp[0][index * metadata.video.INTERVAL + i];
      }
      vp2.cur_vp[index] = metadata.trace.frame_vp[0][index * metadata.video.INTERVAL];
      vp2.est_err[index][0] = 0;
      vp2.est_err[index][1] = 0;	
      break;
  }
  //
  vp2.est_vp[index] = vp2.est_frame_vp[index * metadata.video.INTERVAL];
}
int* AdaptLogic::test(int index){
  int* vmask = NULL;
  int* vmask_tmp;
  int frameID;
  int tid;
  if(index >=  metadata.video.BUFF / metadata.video.INTERVAL){
    vmask = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL]);
    for(frameID = 1; frameID < metadata.video.INTERVAL; frameID ++){
      /* get the visible mask */
      vmask_tmp = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL + frameID]);
      /* update the aggr. visible mask */
      for(tid=0; tid < metadata.vp.No_tile; tid ++){
        if(vmask_tmp[tid] == 1)
          vmask[tid] = 1;
      }
    }
    // printf("#[EXT_ALL]: index=%d vmask:\n", index);
    // showVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v);
  }
  /* Use ROI to find the version of each tile */
  return ROI(metadata.video.TILE_BR[index], metadata.vp.No_tile, metadata.video.NO_VER, vmask, thrp.est_seg_thrp[index]);
}
int* AdaptLogic::test2(int index){
  int* vmask = NULL;
  int* vmask_tmp;
  int* tileVer;
  int frameID;
  int tid;
  double avgVPSNR;
  double VPSNR_thres = 40.0;
  double alpha = 0.8;
  if(index >=  metadata.video.BUFF / metadata.video.INTERVAL){
    vmask = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL]);
    for(frameID = 1; frameID < metadata.video.INTERVAL; frameID ++){
      /* get the visible mask */
      vmask_tmp = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL + frameID]);
      /* update the aggr. visible mask */
      for(tid=0; tid < metadata.vp.No_tile; tid ++){
        if(vmask_tmp[tid] == 1)
          vmask[tid] = 1;
      }
    }
  }
  /* Use ROI to find the version of each tile */
  tileVer =  ROI(metadata.video.TILE_BR[index], metadata.vp.No_tile, metadata.video.NO_VER, vmask, thrp.est_seg_thrp[index]);
  /* Estimate the average V-PSNR of the interval */
  avgVPSNR = 0.0;
  for(frameID = 0; frameID < metadata.video.INTERVAL; frameID ++){
    avgVPSNR += est_vp_psnr(metadata.video.TILE_MSE[index], metadata.vp.No_tile, tileVer, vp2.est_frame_vp[index * metadata.video.INTERVAL + frameID]);
  }
  avgVPSNR /= metadata.video.INTERVAL;
  /* Adjust the tile version if the est.avg.VPSNR exceeds the threshold */
  if(avgVPSNR > VPSNR_thres){
    printf("#[test2] index=%d avgVPSNR=%.2f(dB) VPSNR_thres=%.2f(dB)\n", index, avgVPSNR, VPSNR_thres);
    return ROI_adapt(metadata.video.TILE_BR[index], metadata.vp.No_tile, metadata.video.NO_VER, vmask, thrp.est_seg_thrp[index], alpha);
  }else{
    return tileVer;
  }
}
//int* AdaptLogic::test3(int index){
//int* vmask = NULL;
//int* vmask_tmp;
//int* tileVer = NULL;
//int frameID;
//int tid;
//double alpha = 0.8;
//if(index >=  metadata.video.BUFF / metadata.video.INTERVAL){
//vmask = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL]);
//for(frameID = 1; frameID < metadata.video.INTERVAL; frameID ++){
///* get the visible mask */
//vmask_tmp = get_visible_tile(vp2.est_frame_vp[index * metadata.video.INTERVAL + frameID]);
///* update the aggr. visible mask */
//for(tid=0; tid < metadata.vp.No_tile; tid ++){
//if(vmask_tmp[tid] == 1)
//vmask[tid] = 1;
//}
//}
//}
//// calculate the thresholds
//if(index == 0){
//VPSNR_thres[index] = 0;
//}else{
///* Calculate the instant VPSNR value of last interval if ROI were selected */
//avgVPSNR[index-1] = 0;
//for(frameID = 0; frameID < metadata.video.INTERVAL; frameID ++){
//avgVPSNR[index-1] += est_vp_psnr(metadata.video.TILE_MSE[index-1], metadata.vp.No_tile, tile_ver[index-1], metadata.trace.frame_vp[(index - 1) * metadata.video.INTERVAL + frameID]);
//}
//avgVPSNR[index-1] /= metadata.video.INTERVAL;
//if(index == 1){
//VPSNR_thres[index] = avgVPSNR[index-1];
//}
//}
///* Use ROI to find the version of each tile sothat the avgVPSNR <= VPSNR_thres */
//tileVer =  ROI_adapt_2(metadata.video.TILE_BR[index], metadata.video.TILE_MSE[index], metadata.vp.No_tile, metadata.video.NO_VER, vmask, thrp.est_seg_thrp[index], &vp2.est_frame_vp[index * metadata.video.INTERVAL], VPSNR_thres[index]);
//// 
//return tileVer;
//}
int* AdaptLogic::DASH_ERP(int index){
  int* tile = new int[metadata.vp.No_tile];
  int selectVer;
  selectVer = 0;
  // printf("est_thrp=%.2f\n", thrp.est_seg_thrp[index]);
  // for(int i=0; i < metadata.video.NO_VER; i++)
  // printf("DASH_BR[%d][%d] = %.2f\n", index, i, metadata->DASH_BR[index][i]);
  while(selectVer < metadata.video.NO_VER && metadata.video.DASH_BR[index][selectVer] < thrp.est_seg_thrp[index]) selectVer ++;
  //
  selectVer = (selectVer > 0)?(selectVer-1):selectVer;
  printf("#[DASH_ERP] ver: %d\n", selectVer);
  //
  for(int i=0; i < metadata.vp.No_tile; i++){
    tile[i] = selectVer;
  }
  return tile;
}
int* AdaptLogic::DASH_CMP(int index){
  int* tile = new int[metadata.vp.No_tile];
  int selectVer;
  selectVer = 0;
  printf("est_thrp=%.2f\n", thrp.est_seg_thrp[index]);
  for(int i=0; i < metadata.video.NO_VER; i++){
    printf("DASH_BR[%d][%d] = %.2f PSNR=%.2f\n", index, i, metadata.video.CMP_DASH_BR[index][i],metadata.video.CMP_DASH_PSNR[index][i]);
  }
  while(selectVer < metadata.video.NO_VER && metadata.video.CMP_DASH_BR[index][selectVer] < thrp.est_seg_thrp[index]) selectVer ++;
  //
  selectVer = (selectVer > 0)?(selectVer-1):selectVer;
  printf("#[DASH_CMP] ver: %d\n", selectVer);
  //
  for(int i=0; i < metadata.vp.No_tile; i++){
    tile[i] = selectVer;
  }
  return tile;
}
int* AdaptLogic::EXT_ALL(int index, int ext_width){
  int* vmask = NULL;
  if(index >=  metadata.video.BUFF / metadata.video.INTERVAL){
    vmask = get_visible_tile(vp2.est_vp[index]);
    // printf("#[EXT_ALL]: index=%d vmask:\n", index);
    // showVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v);
    vmask = extVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, ext_width);
    // printf("#[EXT_ALL]: index=%d vmask_est\n", index);
    // showVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v);
  }
  return ROI(metadata.video.TILE_BR[index], metadata.vp.No_tile, metadata.video.NO_VER, vmask, thrp.est_seg_thrp[index]);
}
int* AdaptLogic::ROI(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp){
  int i,j;
  int select_ver;
  int* tile = new int[NO_TILE];
  double sum;
  double delta_R;
  if(vmask == NULL){ // should be changed to EQUAL later
    for(i=0; i < NO_TILE; i++){
      tile[i] = 0;
    }
  }else{
    // allocate the minimum version for tiles not in T_h
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] == 0){
        tile[i] = 0;
        sum += TILE_BR[i][tile[i]];
      }
    }
    delta_R = est_thrp - sum;
    // Apply equal to the remaining tiles
    for(j=0; j < NO_VER; j++){
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          sum += TILE_BR[i][j];
        }
      }
      if(sum > delta_R)
        break;
    }
    //
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] > 0){
        tile[i] = j - 1;
        sum += TILE_BR[i][tile[i]];
      }
    }
    delta_R -= sum;
    //
    bool FLAG = true;
    while(delta_R > 0){
      FLAG = true;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 1){
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R){
            delta_R -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
            tile[i] += 1;
            FLAG = false;
          }
        }
      }
      if(FLAG) break; // break if no improvement is done in the last round.
    }
    while(delta_R > 0){
      FLAG = true;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 0){
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R){
            delta_R -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
            tile[i] += 1;
            FLAG = false;
          }
        }
      }
      if(FLAG) break; // break if no improvement is done in the last round.
    }
  }
  return tile;
}
int* AdaptLogic::ROI_adapt(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp, double alpha){
  int i,j;
  int select_ver;
  int* tile = new int[NO_TILE];
  double sum;
  double delta_R, delta_R_2;
  if(vmask == NULL){ // should be changed to EQUAL later
    for(i=0; i < NO_TILE; i++){
      tile[i] = 0;
    }
  }else{
    // allocate the minimum version for tiles not in T_h
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] == 0){
        tile[i] = 0;
        sum += TILE_BR[i][tile[i]];
      }
    }
    // calculate the remaining bandwidth
    delta_R = est_thrp - sum;
    delta_R_2 = alpha * delta_R;
    // Apply equal to the remaining tiles
    for(j=0; j < NO_VER; j++){
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          sum += TILE_BR[i][j];
        }
      }
      if(sum > delta_R_2)
        break;
    }
    //
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] > 0){
        tile[i] = j - 1;
        sum += TILE_BR[i][tile[i]];
      }
    }
    delta_R_2 -= sum;
    // increase versions of visible tiles (if possible)
    bool FLAG = true;
    while(delta_R > 0){
      FLAG = true;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 1){
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R_2){
            delta_R_2 -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
            tile[i] += 1;
            FLAG = false;
          }
        }
      }
      if(FLAG) break; // break if no improvement is done in the last round.
    }
    // increase versions of invisible tiles
    delta_R = delta_R * (1-alpha);
    while(delta_R > 0){
      FLAG = true;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 0){
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R){
            delta_R -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
            tile[i] += 1;
            FLAG = false;
          }
        }
      }
      if(FLAG) break; // break if no improvement is done in the last round.
    }
  }
  return tile;
}
int* AdaptLogic::ROI_adapt_2(double** TILE_BR, double** TILE_MSE, int NO_TILE, int NO_VER, int* vmask, double est_thrp, int** est_vp, double VPSNR_thres){
  int i,j;
  int select_ver;
  int* tile = new int[NO_TILE];
  double sum;
  double delta_R, delta_R_tmp;
  int max_ver;
  double avgVPSNR;
  int frameID;
  if(vmask == NULL){ // should be changed to EQUAL later
    for(i=0; i < NO_TILE; i++){
      tile[i] = 0;
    }
  }else{
    // allocate the minimum version for invisible tiles
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] == 0){
        tile[i] = 0;
        sum += TILE_BR[i][tile[i]];
      }
    }
    delta_R = est_thrp - sum;
    // Apply equal to the remaining tiles
    for(j=0; j < NO_VER; j++){
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          sum += TILE_BR[i][j];
        }
      }
      if(sum > delta_R)
        break;
    }
    //pp
    max_ver = j;
    // search for the maximum version that satisfy the the VPSNR_thres constraint
    for(j=max_ver; j > 0; j--){
      // allocate the minimum version for invisible tiles
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 0){
          tile[i] = 0;
        }
      }
      delta_R_tmp = delta_R;
      // allocate version 'j' for all visible tiles
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          tile[i] = j - 1;
          sum += TILE_BR[i][tile[i]];
        }
      }
      delta_R_tmp -= sum;
      // improve the version of visible tiles if possible
      bool FLAG = true;
      if(j==max_ver){
        while(delta_R_tmp > 0){
          FLAG = true;
          for(i=0; i < NO_TILE; i++){
            if(vmask[i] == 1){
              if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R_tmp){
                delta_R_tmp -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
                tile[i] += 1;
                FLAG = false;
              }
            }
          }
          if(FLAG) break; // break if no improvement is done in the last round.
        }
      }
      // improve the version of invisible tiles if possible
      while(delta_R_tmp > 0){
        FLAG = true;
        for(i=0; i < NO_TILE; i++){
          if(vmask[i] == 0){
            if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R_tmp){
              delta_R_tmp -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
              tile[i] += 1;
              FLAG = false;
            }
          }
        }
        if(FLAG) break; // break if no improvement is done in the last round.  
      }
      /* estimate the resulting VPSNR values */
      avgVPSNR = 0.0;
      for(frameID = 0; frameID < metadata.video.INTERVAL; frameID ++){
        avgVPSNR += est_vp_psnr(TILE_MSE, metadata.vp.No_tile, tile, est_vp[frameID]);
      }
      avgVPSNR /= metadata.video.INTERVAL;
      if(avgVPSNR < VPSNR_thres) break;
    }
  }
  return tile;
}
int* AdaptLogic::ISM(int index){
  int width;
  int**tileVer = new int*[3];
  int* selectTileVer = NULL;
  double max_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  int INTERVAL = metadata.video.INTERVAL;
  int NO_TILE = metadata.vp.No_tile;
  for(width=0; width <= 2; width ++){
    if(width == 0){
      tileVer[width] = EXT_ALL(index, 0); // ROI
    }else{
      tileVer[width] = ISM(index, width);
    }
    printf("est_vp: (%d, %d)\n", vp2.est_frame_vp[index][0], vp2.est_frame_vp[index][1]);
    tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
    tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr = est_vp_psnr(metadata.video.TILE_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
    tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr_last = est_vp_psnr(metadata.video.TILE_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    //
    printf("index=%d width=%d: %.2f\n",index, width, vp_psnr+vp_psnr_last);
    if((vp_psnr + vp_psnr_last) > max_psnr){
      selectTileVer = tileVer[width];
      max_psnr = vp_psnr + vp_psnr_last;
      decide_width[index] = width;
    }
  }
  printf("index=%d select: max_psnr=%.2f w=%d\n", index, max_psnr, decide_width[index]);
  return selectTileVer;
}
int* AdaptLogic::ISM(int index, int ext_width){
  int i, j, ii, jj;
  int tiles = metadata.vp.No_tile;
  int *tileVer = NULL;
  int mean_ver = metadata.video.NO_VER;
  int mean_ver_2 = metadata.video.NO_VER;
  int vp_ver;
  int ext_1_ver;
  double BW = thrp.est_seg_thrp[index];
  double remainBW;
  bool FLAG;
  int *selectTileVer = new int[tiles];
  double max_vp_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  //
  int width = 0;
  int INTERVAL = metadata.video.INTERVAL;
  double** TILE_MSE = metadata.video.TILE_MSE[index];
  double** TILE_BR = metadata.video.TILE_BR[index];
  int NO_TILE = metadata.vp.No_tile;
  int NO_VER = metadata.video.NO_VER;
  //
  printf("#[ISM][seg #%d]: ext_width:%d\n", index, ext_width);
  //
  int *vmask = NULL;
  if(index >=  metadata.video.BUFF / metadata.video.INTERVAL){
    vmask = get_visible_tile(vp2.est_vp[index]);
    vmask = extVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, ext_width);
  }
  // get the selected version when all tile are of equal quality
  tileVer = ROI(metadata.video.TILE_BR[index
      ], NO_TILE, NO_VER, vmask, thrp.est_seg_thrp[index]);
  if(vmask == NULL)
    return tileVer;

  // calculate viewport-psnr given 'tileVer'
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  // first frame & last frame
  tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
  tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
  if(tmp_vp[0] >= 180)
    tmp_vp[0] -= 360;
  if(tmp_vp[0] < -180)
    tmp_vp[0] += 360;
  //
  if(tmp_vp[1] >= 90)
    tmp_vp[1] -= 180;
  if(tmp_vp[1] <= -90)//
    tmp_vp[1] += 180;
  vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
  tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
  tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
  if(tmp_vp[0] >= 180)
    tmp_vp[0] -= 360;
  if(tmp_vp[0] < -180)
    tmp_vp[0] += 360;
  //
  if(tmp_vp[1] >= 90)
    tmp_vp[1] -= 180;
  if(tmp_vp[1] <= -90)//
    tmp_vp[1] += 180;
  vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
  // update selection
  if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
    std::copy (tileVer, tileVer + tiles, selectTileVer);
    max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
  }
  //p
  if(ext_width == 1){
    // determine quality of visible part
    for(i=0; i < tiles; i++)
      if(vmask[i] > 0 && tileVer[i] < mean_ver)
        mean_ver = tileVer[i];
    // search for the best (viewport-version, extesion-version) pair
    vp_ver = mean_ver + 1;
    while(vp_ver < NO_VER){
      remainBW = BW;
      // assign quality for visible and non-visible tiles
      for(i=0; i < tiles; i++){
        if(vmask[i] == 1) tileVer[i] = vp_ver; // assign a quality of 'vp_ver' for visible tiles
        else
          tileVer[i] = 0; // minimum quality for non-visible tiles
        remainBW -= TILE_BR[i][tileVer[i]]; // update remaining bandwidth
      }
      //
      if(remainBW <= 0){
        // reduce vp tiles's version to meet the bandwidth constraint
        for(i=0; i < tiles; i++){
          if(vmask[i] == 1){
            tileVer[i] -= 1; // assume that all tiles have same priority
            remainBW += TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
            if(remainBW > 0)
              break;
          }
        }
        // calculate viewport-psnr given 'tileVer'
        int tmp_vp[2];
        // first frame & last frame
        tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
        tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
        if(tmp_vp[0] >= 180)
          tmp_vp[0] -= 360;
        if(tmp_vp[0] < -180)
          tmp_vp[0] += 360;
        //
        if(tmp_vp[1] >= 90)
          tmp_vp[1] -= 180;
        if(tmp_vp[1] <= -90)//
          tmp_vp[1] += 180;
        vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
        tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
        tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
        if(tmp_vp[0] >= 180)
          tmp_vp[0] -= 360;
        if(tmp_vp[0] < -180)
          tmp_vp[0] += 360;
        //
        if(tmp_vp[1] >= 90)
          tmp_vp[1] -= 180;
        if(tmp_vp[1] <= -90)//
          tmp_vp[1] += 180;
        vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
        // update selection
        if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
          std::copy (tileVer, tileVer + tiles, selectTileVer);
          max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
        }
        break; // not enough bandwidth
      }
      // assign remaining bandwidth for ext_1 tiles
      FLAG = false;
      while(remainBW > 0 && !FLAG){
        FLAG = true;
        for(i=0; i < NO_TILE; i++){
          if(vmask[i] == 2 && remainBW > 0){
            if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
              remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
              tileVer[i] ++;
              FLAG = false;
            }
          }
        }
      }
      // assign remaining bandwidth to non-visible tiles
      if(remainBW > 0){
        FLAG = false;
        while(remainBW > 0 && !FLAG){
          FLAG = true;
          for(i=0; i < NO_TILE; i++){
            if(vmask[i] == 0 && remainBW > 0){
              if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
                remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
                tileVer[i] ++;
                FLAG = false;
              }
            }
          }
        }
      }
      // calculate viewport-psnr given 'tileVer'
      // first frame & last frame
      tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
      tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
      if(tmp_vp[0] >= 180)
        tmp_vp[0] -= 360;
      if(tmp_vp[0] < -180)
        tmp_vp[0] += 360;
      //
      if(tmp_vp[1] >= 90)
        tmp_vp[1] -= 180;
      if(tmp_vp[1] <= -90)//
        tmp_vp[1] += 180;
      vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
      tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
      if(tmp_vp[0] >= 180)
        tmp_vp[0] -= 360;
      if(tmp_vp[0] < -180)
        tmp_vp[0] += 360;
      //
      if(tmp_vp[1] >= 90)
        tmp_vp[1] -= 180;
      if(tmp_vp[1] <= -90)//
        tmp_vp[1] += 180;
      vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      // update selection
      if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
        std::copy (tileVer, tileVer + tiles, selectTileVer);
        max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
      }
      //
      vp_ver ++;
    }//endwhile
  }
  else{
    // determine quality of visible part
    for(i=0; i < tiles; i++)
      if(vmask[i] > 0 && tileVer[i] < mean_ver)
        mean_ver = tileVer[i];
    // search for the best (viewport-version, ext1-version, ext2-version) tube
    vp_ver = mean_ver + 1;
    while(vp_ver < NO_VER){
      remainBW = BW;
      // assign quality for visible and non-visible tiles
      for(i=0; i < tiles; i++){
        if(vmask[i] == 1) tileVer[i] = vp_ver; // assign a quality of 'vp_ver' for visible tiles
        else
          tileVer[i] = 0; // minimum quality for non-visible tiles
        remainBW -= TILE_BR[i][tileVer[i]]; // update remaining bandwidth
      }
      //
      if(remainBW <= 0) break;
      // determine quality if 'ext1' and 'ext2' are assigned same version
      // equally assign remaining bandwidth for 'ext_1' and 'ext_2'
      FLAG = false;
      while(remainBW > 0 && !FLAG){
        FLAG = true;
        for(i=0; i < NO_TILE; i++){
          if((vmask[i] == 2 || vmask[i] == 3) && remainBW > 0){
            if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
              remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
              tileVer[i] ++;
              FLAG = false;
            }
          }
        }
      }
      // find version
      for(i=0; i < tiles; i++)
        if((vmask[i] == 2 || vmask[i] == 3) && tileVer[i] < mean_ver_2)
          mean_ver_2 = tileVer[i];
      // assign remaining bandwidth to non-visible tiles
      if(remainBW > 0){
        FLAG = false;
        while(remainBW > 0 && !FLAG){
          FLAG = true;
          for(i=0; i < NO_TILE; i++){
            if(vmask[i] == 0 && remainBW > 0){
              if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
                remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
                tileVer[i] ++;
                FLAG = false;
              }
            }
          }
        }
      }
      //
      // calculate viewport-psnr given 'tileVer'
      // first frame & last frame
      tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
      tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
      if(tmp_vp[0] >= 180)
        tmp_vp[0] -= 360;
      if(tmp_vp[0] < -180)
        tmp_vp[0] += 360;
      //
      if(tmp_vp[1] >= 90)
        tmp_vp[1] -= 180;
      if(tmp_vp[1] <= -90)//
        tmp_vp[1] += 180;
      vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
      tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
      if(tmp_vp[0] >= 180)
        tmp_vp[0] -= 360;
      if(tmp_vp[0] < -180)
        tmp_vp[0] += 360;
      //
      if(tmp_vp[1] >= 90)
        tmp_vp[1] -= 180;
      if(tmp_vp[1] <= -90)//
        tmp_vp[1] += 180;
      vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      // update selection
      if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
        std::copy (tileVer, tileVer + tiles, selectTileVer);
        max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
      }
      //
      ext_1_ver = mean_ver_2 + 1;
      while(ext_1_ver <= vp_ver){
        remainBW = BW;
        // assign quality for visible and non-visible tiles
        for(i=0; i < tiles; i++){
          if(vmask[i] == 1) tileVer[i] = vp_ver; // assign a quality of 'vp_ver' for visible tiles
          else
            if(vmask[i] == 2) tileVer[i] = ext_1_ver; // assign quality for ext_1 tiles
            else
              tileVer[i] = 0;
          remainBW -= TILE_BR[i][tileVer[i]]; // update remaining bandwidth
        }
        //
        if(remainBW <= 0) break;
        // assign remaining bandwidth for 'ext_2' tiles
        FLAG = false;
        while(remainBW > 0 && !FLAG){
          FLAG = true;
          for(i=0; i < NO_TILE; i++){
            if(vmask[i] == 3 && remainBW > 0){
              if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
                remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
                tileVer[i] ++;
                FLAG = false;
              }
            }
          }
        }
        // assign remaining bandwidth to non-visible tiles
        if(remainBW > 0){
          FLAG = false;
          while(remainBW > 0 && !FLAG){
            FLAG = true;
            for(i=0; i < NO_TILE; i++){
              if(vmask[i] == 0 && remainBW > 0){
                if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
                  remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
                  tileVer[i] ++;
                  FLAG = false;
                }
              }
            }
          }
        }
        // calculate viewport-psnr given 'tileVer'
        // first frame & last frame
        tmp_vp[0] = vp2.est_frame_vp[index * INTERVAL][0] + vp2.est_err[index][0];
        tmp_vp[1] = vp2.est_frame_vp[index * INTERVAL][1] + vp2.est_err[index][1];
        if(tmp_vp[0] >= 180)
          tmp_vp[0] -= 360;
        if(tmp_vp[0] < -180)
          tmp_vp[0] += 360;
        //
        if(tmp_vp[1] >= 90)
          tmp_vp[1] -= 180;
        if(tmp_vp[1] <= -90)//
          tmp_vp[1] += 180;
        vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
        tmp_vp[0] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][0] + vp2.est_err[index][0];
        tmp_vp[1] = vp2.est_frame_vp[(index + 1) * INTERVAL - 1][1] + vp2.est_err[index][1];
        if(tmp_vp[0] >= 180)
          tmp_vp[0] -= 360;
        if(tmp_vp[0] < -180)
          tmp_vp[0] += 360;
        //
        if(tmp_vp[1] >= 90)
          tmp_vp[1] -= 180;
        if(tmp_vp[1] <= -90)//
          tmp_vp[1] += 180;
        vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
        // update selection
        if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
          std::copy (tileVer, tileVer + tiles, selectTileVer);
          max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
        }
        //
        ext_1_ver ++;
      }
      vp_ver++;
    }
  }
  return selectTileVer;
}
int* AdaptLogic::Ireland(int index){
  int i,j, ext_w, k;
  int select_ver;
  int* tile_ver = new int[metadata.vp.No_tile];
  double* tile_w = new double[metadata.vp.No_tile];
  double* tile_br = new double[metadata.vp.No_tile];
  double* tile_d = new double[metadata.vp.No_tile];
  double* k_i = new double[metadata.vp.No_tile];
  int* vmask, *vmask_cur, *pixel;
  int* vmask_ext;
  double BR_budget = thrp.est_seg_thrp[index];
  double BR_used = 0;
  int NT = metadata.vp.No_tile;
  int NT_W = metadata.vp.No_tile_v;
  int NT_H = metadata.vp.No_tile_h;
  int W_T = metadata.vp.tile_W;
  int W_H = metadata.vp.tile_H;
  int VP_T = 960;
  int VP_H = 960;
  int W = 3840;
  int H = 1920;
  double lamda = 0.8, u, v, m, n;
  int tid;
  /* marking */
  vmask = get_visible_tile(vp2.est_vp[index]);
  pixel = get_visible_pixel(vp2.est_vp[index]);
  // printf("#[Ireland][index=%d] visible budget=%.2f (%d, %d)\n", index, BR_budget, vp2.est_vp[index][0], vp2.est_vp[index][1]);
  // printf("#[Ireland][index=%d] vmask:\n", index);
  showArrayInt(vmask, NT);
  // printf("#[Ireland][index=%d] pixel:\n", index);
  showArrayInt(pixel, NT);
  /* assign bitrate for visible tiles */
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    if(vmask[i] == 1){
      /* calculate tile weight and allocated bitrate*/
      tile_w[i] = pixel[i] * 1.0 / (VP_T * VP_H);
      tile_br[i] = lamda * BR_budget * tile_w[i];
      /* decide the version of each tile */
      j = 0;
      while(j < metadata.video.NO_VER && metadata.video.TILE_BR[index][i][j] < tile_br[i]) j++;
      if(j==0){ /*allocated bitrate is not enough even for the lowest version */
        // printf("#[Ireland] tile #%d\n", i);
        tile_ver[i] = 0;//
      }else{
        tile_ver[i] = j-1;
      }
      BR_used += metadata.video.TILE_BR[index][i][tile_ver[i]];
    }
  }
  // printf("#[Ireland][index=%d] Visible tiles version:\n", index);
  showArrayInt(tile_ver, NT);
  // printf("#[Ireland][index=%d] visible budget=%.2f used=%.2f\n", index,lamda * BR_budget, BR_used);
  // for(i=0; i < NT_H; i++){
  //   for(j=0; j < NT_W; j++)
  // printf("%d ", tile_ver[i*NT_W + j]);
  // printf("\n");
  // }
  /* calculate distance to the viewport center of each outside tiles */
  double max_d = 0;
  double sum_d = 0;
  /* convert (phi,theta) to ERP's location */
  // u = vp2.est_vp[index][0]/(2*M_PI) + 0.5;
  // v = 0.5 - vp2.est_vp[index][1]/M_PI;
  // m = u * W  - 0.5;
  // n = v * H - 0.5;
  // for(k=0; k < metadata.vp.No_face; k++){
  //  for(i=0; i < NT_H; i++){
  //    for(j=0; j < NT_W; j++){
  //      if(vmask[k* NT_H * NT_W + i * NT_W + j] == 0){
  // 	/* calculate distace to viewport center*/
  // 	tile_d[k* NT_H * NT_W + i * NT_W + j] = calc_distance(vp2.est_vp[index], k, i, j, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, metadata.vp.face_W, metadata.vp.face_H);
  // 	/* find the max. distance */
  // 	if(tile_d[i*NT_W + j] > max_d)
  //  	 max_d = tile_d[i*NT_W + j];//
  //      }
  //    }
  //  }
  // }
  for(tid = 0; tid < NT; tid++){
    if(vmask[tid] == 0){
      tile_d[tid] = calc_distance(vp2.est_vp[index], tid, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, metadata.vp.face_W, metadata.vp.face_H);
      /* find the longest distance */
      if(tile_d[tid] > max_d)
        max_d = tile_d[tid];
    }else{
      tile_d[tid] = 0; // distance is zero if it is a visible tile
    }
  }
  // printf("#[Ireland][index=%d] Invisible tiles distances:\n", index);
  showArrayDouble(tile_d, NT);
  /* calculate  k_i */
  for(i=0; i < NT; i++){
    if(vmask[i] == 0){
      k_i[i] = max_d * 1.0 / tile_d[i];
      sum_d += k_i[i]; 
    }
  }
  BR_used = 0;
  /* assign the bitrate for outside tiles */
  for(i=0; i < NT; i++){
    if(vmask[i] == 0){
      tile_br[i] = (1-lamda) * BR_budget * k_i[i]/sum_d;
      j = 0;
      while(j < metadata.video.NO_VER && metadata.video.TILE_BR[index][i][j] < tile_br[i]) j++;
      if(j==0){ /*allocated bitrate is not enough even for the lowest version */
        tile_ver[i] = 0;//
      }else{
        tile_ver[i] = j-1;
      }
      BR_used += metadata.video.TILE_BR[index][i][tile_ver[i]];
    }
  }
  // printf("#[Ireland][index=%d] outside budget=%.2f used=%.2f\n", index,(1-lamda) * BR_budget, BR_used);
  // printf("#[Ireland][index=%d] Invisible tiles version:\n", index);
  // showArrayDouble(tile_br, NT);
  showArrayInt(tile_ver, NT);
  // for(i=0; i < NT; i++)
  // printf("%.2f\t", metadata.video.TILE_BR[index][i][tile_ver[i]]);
  // printf("\n");
  // printf("#[Ireland][index=%d] Invisible tiles distances:\n", index);
  // for(i=0; i < NT_H; i++){
  //   for(j=0; j < NT_W; j++)
  //     printf("%d ", tile_ver[i*NT_W + j]);
  //   printf("\t");
  //   for(j=0; j < NT_W; j++)
  //     printf("%.2f ", metadata.video.TILE_BR[index][i*NT_W + j][tile_ver[i*NT_W + j]]);
  //   printf("\t");
  //   for(j=0; j < NT_W; j++)
  //     printf("%.2f ", tile_br[i*NT_W + j]);
  //   printf("\t");
  //   for(j=0; j < NT_W; j++)
  //     printf("%.2f ", tile_d[i*NT_W + j]);
  //   printf("\n");
  // }
  return tile_ver;
}
int* AdaptLogic::Ghent(int index){
  int i,j, ext_w;
  int select_ver;
  int* vmask, *vmask_cur;
  int* vmask_ext;
  double BR_budget = thrp.est_seg_thrp[index];
  int NT = metadata.vp.No_tile;
  int NT_W = metadata.vp.No_tile_v;
  int NT_H = metadata.vp.No_tile_h;
  int* tile_ver = new int[NT];

  printf("#[Ghent] index=%d\n", index);
  // marking tiles
  ext_w = 1;
  vmask = get_visible_tile(vp2.est_vp[index]);
  /* add visible tiles of the current viewport to 'vmask'*/
  vmask_cur = get_visible_tile(vp2.cur_vp[index]);
  for(i=0; i < NT; i++)
    if(vmask_cur[i] == 1)
      vmask[i] = 1;
  /* extend vmask by 1 tile in all directions --> adjcent tiles*/
  vmask_ext = extVmask(vmask, metadata.vp.No_face, metadata.vp.No_tile_h, metadata.vp.No_tile_v, 1);
  //check point #1:
  // printf("#[Ghent][index=%d] BR_budget=%.2f, est_vp=(%d,%d), est_vp=(%d, %d)\n", index, BR_budget, vp2.est_vp[index][0], vp2.est_vp[index][1], vp2.cur_vp[index][0], vp2.cur_vp[index][1]);
  // for(i=0; i < NT_H; i++){
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", vmask[i*NT_W + j]);
  // 	printf("\t");
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", vmask_cur[i*NT_W + j]);
  // 	printf("\t");
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", vmask_ext[i*NT_W + j]);
  // 	printf("\n");
  // }
  // assigning tile's bitrates
  /* assign the lowest version to all tiles*/
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    BR_budget -= metadata.video.TILE_BR[index][i][tile_ver[i]]; 
  }
  if(BR_budget < 0) /* return if avail. bw is not enough*/
    return NULL;
  /* assign highest possible version for viewport tiles */
  double sum_br;
  for(j=1; j < metadata.video.NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1)
        sum_br += metadata.video.TILE_BR[index][i][j] - metadata.video.TILE_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1){
        BR_budget -= (metadata.video.TILE_BR[index][i][j-1] - metadata.video.TILE_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  //check point #2:
  // printf("#[Ghent][index=%d] viewport tiles, BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", tile_ver[i*NT_W + j]);
  // 	printf("\n");
  // }
  /* assign highest possible version for adjcent tiles */
  for(j=1; j < metadata.video.NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2)
        sum_br += metadata.video.TILE_BR[index][i][j] - metadata.video.TILE_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2){
        BR_budget -= (metadata.video.TILE_BR[index][i][j-1] - metadata.video.TILE_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  // printf("#[Ghent][index=%d] adjecent BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", tile_ver[i*NT_W + j]);
  // 	printf("\n");
  // }
  /* assign highest possible version for outside tiles */
  for(j=1; j < metadata.video.NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0)
        sum_br += metadata.video.TILE_BR[index][i][j] - metadata.video.TILE_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0){
        BR_budget -= (metadata.video.TILE_BR[index][i][j-1] - metadata.video.TILE_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  // printf("#[Ghent][index=%d] outside BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  // 	for(j=0; j < NT_W; j++)
  // 		printf("%d ", tile_ver[i*NT_W + j]);
  // 	printf("\n");
  // }
  return tile_ver;
}
int* AdaptLogic::get_visible_tile(int* vp){
  int phi,theta;
  if(vp[0] < 0)
    phi = vp[0] + 360;
  else
    phi = vp[0];
  // convert from [-90; 90] to [0; 180]
  theta = vp[1] + 90;
  // printf("#[get_visible_tile]: (%d, %d) -> (%d, %d)\n", vp[0], vp[1], phi, theta);
  return metadata.vp.vmask[phi][theta]; 
}
int* AdaptLogic::get_visible_pixel(int* vp){
  int phi,theta;
  if(vp[0] < 0)
    phi = vp[0] + 360;
  else
    phi = vp[0];
  // convert from [-90; 90] to [0; 180]
  theta = vp[1] + 90;
  // printf("#[getVisibleTile]: (%d, %d) -> (%d, %d)\n", vp[0], vp[1], phi, theta);
  return metadata.vp.pixel[phi][theta]; 
}
double AdaptLogic::est_vp_psnr(double** TILE_MSE, int NO_TILE, int* tileVer, int* est_vp){
  double v_psnr = 0;
  int* vmask = get_visible_tile(est_vp);
  if(vmask == NULL){
    printf("vmask null\n");
    printf("vp: (%d, %d)\n", est_vp[0], est_vp[1]);
  }
  int* pixel = get_visible_pixel(est_vp);
  if(pixel == NULL){
    printf("(phi, theta) = (%d, %d)\n", est_vp[0], est_vp[1]);
    printf("pixel null\n");
  }
  // estimate viewport MSE as a weighted sum of tiles' MSE. Tiles' weight = % of tile in viewport
  double est_vp_mse = 0;
  for(int i=0; i < metadata.vp.No_tile; i++){
    // printf("%d %d %d %d\n", index, i, tileVer[i], pixel[i]);
    if(vmask[i] > 0){
      // printf("#[est_vp_psnr] tile:%d\n", i);
      // printf("#[est_vp_psnr] pixel: %d\n", pixel[i]);
      // printf("#[est_vp_psnr] tileVer: %d\n", tileVer[i]);
      // printf("#[est_vp_psnr] mse: %.2f\n", TILE_MSE[i][tileVer[i]]);
      // printf("#[est_vp_psnr] tile:%d mse: %.2f pixel: %d\n", i, TILE_MSE[i][tileVer[i]], pixel[i]);
      est_vp_mse += TILE_MSE[i][tileVer[i]] * pixel[i] * 1.0 / (metadata.vp.vp_W * metadata.vp.vp_H);
    }
  }
  v_psnr = 10 * log10((255*255)/est_vp_mse);
  // printf("%.2f\t%.2f\n", est_vp_mse, est_vp_psnr_first);
  return v_psnr;
}
double AdaptLogic::calc_distance(int* vp, int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H){
  double dis, u, v, m, n;
  int W_T = face_W / No_tile_h;
  int W_H = face_H / No_tile_v;
  int W  = face_W * No_face;
  int H  = face_H * No_face;
  double face_phi, face_theta;
  /* calculate the spherical position (phi, theta) of the tile's centre*/
  get_tile_center_pos(tid, No_face, No_tile_h, No_tile_v, face_W, face_H, &face_phi, &face_theta);
  /* calculate the distance to the viewport's centre */
  dis = get_Euclide_distance(vp[0], vp[1], face_phi, face_theta);
  return dis;
}
double AdaptLogic::get_Euclide_distance(double phi, double theta, double phi2, double theta2){
  double dist;
  double x,y,z, x2, y2, z2;
  //
  x = cos(theta) * cos(phi);
  y = sin(theta);
  z = -cos(theta) * sin(phi);
  //
  x2 = cos(theta2) * cos(phi2);
  y2 = sin(theta2);
  z2 = -cos(theta2) * sin(phi2);
  //
  dist = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + (z-z2)*(z-z2));
  printf("#[get_Euclide_distance] dis: %.2f (x,y,z)=(%.2f,%.2f,%.2f) (x2,y2,z2)=(%.2f,%.2f,%.2f)\n",dist,x,y,z,x2,y2,z2);
  return dist;
}
void AdaptLogic::get_tile_center_pos(int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H, double* phi, double *theta){
  double m,n,u,v;
  double x,y,z;
  int tile_W = face_W / No_tile_h;
  int tile_H = face_H / No_tile_v;
  int face_id;
  int tile_id;
  int tile_id_h;
  int tile_id_v;
  // convert from tile id to (face_id, tid_h, tid_v)
  get_face_tid(No_face, No_tile_h, No_tile_v, tid, &face_id, &tile_id);
  tile_id_h = tile_id % No_tile_h;
  tile_id_v = tile_id / No_tile_h;
  //
  if(No_face == 2){ // ERP-2+MxN
    if(tid == 0 || tid == 1){
      m = face_W / 2;
      n = (tid == 0)?(metadata.vp.erp_H - face_H)/4:((metadata.vp.erp_H - face_H)*3/4 + face_H);
    }else{
      tile_id_v = (tid - 2) / No_tile_h;
      tile_id_h = (tid - 2) % No_tile_h;
      m = tile_id_h * tile_W + tile_W/2;
      n = (metadata.vp.erp_H - face_H)/2 + tile_id_v * tile_H + tile_H/2;
    }
    u = (m+0.5)/metadata.vp.erp_W;
    v = (n+0.5)/metadata.vp.erp_H;
    *phi = (u-0.5) * 2 * M_PI;
    *theta = (0.5 - v) * M_PI;
    printf("#[get_tile_center_pos]: tid:%d m=%.0f n=%.0f phi=%.2f theta=%.2f\n", tid, m,n,*phi, *theta);
  }
  if(No_face == 1){ // ERP-MxN
    m = tile_id_h * tile_W + tile_W/2;
    n = tile_id_v * tile_H + tile_H/2;
    u = (m+0.5)/face_W;
    v = (n+0.5)/face_H;
    *phi = (u-0.5) * 2 * M_PI;
    *theta = (0.5 - v) * M_PI;
  }else{
    if(No_face == 6){ // CMP
      m = tile_id_h * tile_W + tile_W/2;
      n = tile_id_v * tile_H + tile_H/2;
      u = (m+0.5)*2/face_W -1;
      v = (n+0.5)*2/face_H -1;
      switch(face_id){
        case 0:
          x = 1.0;
          y = -v;
          z = -u;
          break;
        case 1:
          x = -1.0;
          y = -v;
          z = u;
          break;
        case 2:
          x = u;
          y = 1.0;
          z = v;
          break;
        case 3:
          x = u;
          y = -1.0;
          z = -v;
          break;
        case 4:
          x = u;
          y = -v;
          z = 1.0;
          break;
        case 5:
          x = -u;
          y = -v;
          z = -1.0;
          break;
      }
      *phi = atan2(-z, x);
      *theta = asin(y/sqrt(x*x + y*y + z*z));
      printf("#[get_tile_center_pos]: tid:%d phi=%.2f theta=%.2f tid_h=%d tid_v=%d\n", tid, *phi * 180.0 / M_PI, *theta * 180.0 / M_PI, tile_id_h, tile_id_v);
    }
  }
}
