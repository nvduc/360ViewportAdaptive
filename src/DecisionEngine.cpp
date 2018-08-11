#include "DecisionEngine.h"
#include "Metadata.h"
#include "Projection.h"
#include "common.h"
#include "GaussianFilter.h"
#include "gurobi_c++.h"
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
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <malloc.h>
// File in 'exper' branch
DecisionEngine::DecisionEngine(Metadata* meta, AdaptInfo adaptInfo){
  int i,j,k, tid;
  metadata = meta;
  INTER = adaptInfo.INTER;
  // yakitori_erp_896_frame
  switch(INTER){
    case 4:
      INTER_ID = 0;
      break;
    case 16:
      INTER_ID = 1;
      break;
    case 32:
      INTER_ID = 2;
      break;
    case 64:
      INTER_ID = 3;
      break;
  }
  // Timelapse_NY_24_frame
  /*
  switch(INTER){
    case 24:
      INTER_ID = 0;
  }
  */
  BW = adaptInfo.BW;
  BUFF = adaptInfo.BUFF;
  TILING = adaptInfo.TILING;
  NO_SEG = metadata->video_info.NO_FRAME / INTER;
  NO_VER = metadata->video_info.NO_VER;
  HTRACE_ID = adaptInfo.HTRACE_ID;
  No_tile = metadata->video_info.tile[TILING].No_tile;
  No_tile_h = metadata->video_info.tile[TILING].No_tile_h; 
  No_tile_v = metadata->video_info.tile[TILING].No_tile_v; 
  NO_FRAME = metadata->video_info.NO_FRAME;
  NO_FRAME_ORIGIN = metadata->video_info.NO_FRAME_ORIGIN;
  FPS = metadata->video_info.FPS;
  No_face = metadata->video_info.tile[TILING].No_face;
  tile_W = metadata->video_info.tile[TILING].tile_W;
  tile_H = metadata->video_info.tile[TILING].tile_H;
  FoV = metadata->video_info.tile[TILING].FoV;
  vmask = metadata->video_info.tile[TILING].vmask;
  pixel = metadata->video_info.tile[TILING].pixel;
  TILE_SEG_BR = metadata->video_info.tile[TILING].TILE_SEG_BR[INTER_ID];
  TILE_SEG_SIZE = metadata->video_info.tile[TILING].TILE_SEG_SIZE[INTER_ID];
  TILE_SEG_MSE = metadata->video_info.tile[TILING].TILE_SEG_MSE[INTER_ID];
  TILE_SEG_PSNR = metadata->video_info.tile[TILING].TILE_SEG_PSNR[INTER_ID];
  No_user = metadata->video_info.htrace.No_user;
  htrace = metadata->video_info.htrace.trace[HTRACE_ID];
  TILE_SELECT_METHOD = adaptInfo.TILE_SELECT_METHOD;
  VP_EST_METHOD = adaptInfo.VP_EST_METHOD;
  face_W = metadata->video_info.tile[TILING].face_W;
  face_H = metadata->video_info.tile[TILING].face_H;
  vp_H = metadata->video_info.tile[TILING].vp_H;
  vp_W = metadata->video_info.tile[TILING].vp_W;
  W = metadata->video_info.W;
  H = metadata->video_info.H;
  B_max = BUFF*1.0/FPS;
  B_target = B_max;
  B_cri = B_max;
  // DASH-info
  DASH_SEG_BR = init2dArrayDouble(NO_SEG, NO_VER);
  DASH_SEG_PSNR = init2dArrayDouble(NO_SEG, NO_VER);
  for(i=0; i < NO_SEG; i++){
    for(j=0; j < NO_VER; j++){
      DASH_SEG_BR[i][j] = metadata->video_info.tile[0].TILE_SEG_BR[INTER_ID][i][0][j];
      DASH_SEG_PSNR[i][j] = metadata->video_info.tile[0].TILE_SEG_PSNR[INTER_ID][i][0][j];
    }
  }
  // load DASH's tile MSE values 
  DASH_TILE_PSNR = init3dArrayDouble(NO_SEG, No_tile, NO_VER);
  DASH_TILE_MSE = init3dArrayDouble(NO_SEG, No_tile, NO_VER);
  // BellLab
  ang_speed = new double[NO_SEG];
  net_delay = new double[NO_SEG];
  Uti = init3dArrayDouble(NO_SEG, No_tile, NO_VER);       // utility
  Cost = init3dArrayDouble(NO_SEG, No_tile, NO_VER);      // Cost
  for(tid=0; tid < No_tile; tid++){
    for(i=0; i < NO_SEG; i++){
      for(j=0; j < NO_VER; j++){
        if(j==0){
          Uti[i][tid][j] = TILE_SEG_MSE[i][tid][j] - 65025;
          Cost[i][tid][j] = TILE_SEG_BR[i][tid][j];
        }else{
          Uti[i][tid][j] = TILE_SEG_MSE[i][tid][j-1] - TILE_SEG_MSE[i][tid][j];
          Cost[i][tid][j] = TILE_SEG_BR[i][tid][j] - TILE_SEG_BR[i][tid][j-1];
        }
      }
    }
  }
  // only for version 2 (QP=40)
  vector <double> v;
  int rows;
  int cols;
  int ret;
  char buff[1024];
  for(i=0; i < No_tile; i++){
    sprintf(buff, "data/tile/%s_%dx%d_%dfr_%dfps/tile_1x1/%dframe/tile_seg_psnr/tile_%d.txt", metadata->video_info.name.c_str(), W, H, NO_FRAME_ORIGIN,FPS,INTER,i);
    ret = metadata->import_matrix_from_txt_file(buff, v, rows, cols);
    if(ret != 0){
      printf("#[load_visible_mask] Cannot open file !\n");
      exit(1);
    }
    // for(j=0; j < NO_SEG; j++){
    //  for(k=0; k < NO_VER; k++){
    //    switch(k){
    //      case 0:
    //        DASH_TILE_PSNR[j][i][k] = 0;
    //        break;
    //      case 1:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols + 4];
    //        break;
    //      case 2:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols + 5];
    //        break;
    //      case 3:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols + 3];
    //        break;
    //      case 4:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols + 2];
    //        break;
    //      case 5:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols + 1];
    //        break;
    //      case 6:
    //        DASH_TILE_PSNR[j][i][k] = v[j*cols];
    //        break;
    //    }
    //    DASH_TILE_MSE[j][i][k] = 255.0 * 255.0 / (pow(10, DASH_TILE_PSNR[j][i][k]/10.0));
    //  }
    // }
    for(j=0; j < NO_FRAME_ORIGIN/INTER; j++){
      for(k=0; k < NO_VER; k++){
        DASH_TILE_PSNR[j][i][k] = v[j*cols + NO_VER - 1 - k]; // should be modifed depenending on the log
        DASH_TILE_MSE[j][i][k] = 255.0 * 255.0 / (pow(10, DASH_TILE_PSNR[j][i][k]/10.0));
      }
    }
    // repeat during the session
    for(j=NO_FRAME_ORIGIN/INTER; j < NO_SEG; j++){
      for(k=0; k < NO_VER; k++){
        DASH_TILE_PSNR[j][i][k] = DASH_TILE_PSNR[j%(NO_FRAME_ORIGIN/INTER)][i][k];
        DASH_TILE_MSE[j][i][k] = DASH_TILE_MSE[j%(NO_FRAME_ORIGIN/INTER)][i][k];
      }
    }
  }
 // compare tile-wise psnr of different versions
  // printf("=========\n");
  // for(k=0; k < NO_VER; k++){
  //   for(i=0; i < No_tile_v; i++){
  //     // DASH
  //     for(j=0; j < No_tile_h; j++){
  //       printf("%.2f ", DASH_TILE_PSNR[0][i * No_tile_v + j][k]);
  //     }
  //     printf("\t");
  //     // tiling
  //     for(j=0; j < No_tile_h; j++){
  //       printf("%.2f ", TILE_SEG_PSNR[0][i * No_tile_v + j][k]);
  //     }
  //     printf("\n");
  //   }
  //   printf("=========\n");
  // }
  // calculate average versions' bitrates
  VER_AVG_BR = new double[NO_VER];
  for(k=0; k < NO_VER; k++){
    VER_AVG_BR[k] = 0;
    for(i=0; i < NO_SEG; i++){
      for(j=0; j < No_tile; j++){
        VER_AVG_BR[k] += TILE_SEG_SIZE[i][j][k];
      }
    }
    VER_AVG_BR[k] /= (No_tile * NO_FRAME * 1.0/FPS);
    printf("#[DecisionEngine] init() ver #%d=%.2f\n", k, VER_AVG_BR[k]);
  }
  // throughput
  est_thrp_act = new double[NO_SEG];
  seg_thrp = new double[NO_SEG];
  est_seg_thrp = new double[NO_SEG];
  // buffer
  buff_level = new double[NO_SEG];
  stall_time = new double[NO_SEG];
  seg_down_time = new double[NO_SEG];
  seg_down_start_time = new double[NO_SEG];
  seg_down_finis_time = new double[NO_SEG];
  time_to_next_request = new double[NO_SEG];
  B_0 = 1; // initial buffering
  // viewport2
  cur_vp = init2dArrayInt(NO_SEG, 2);
  est_vp = init2dArrayInt(NO_SEG, 2);
  est_err = init2dArrayInt(NO_SEG, 2);
  est_err_frame = init2dArrayInt(NO_FRAME, 2);
  est_frame_vp = init2dArrayInt(NO_FRAME, 2);
  est_frame_vp_2 = init2dArrayInt(NO_FRAME, 2);
  speed = init2dArrayDouble(NO_SEG, 2);
  for(i=0; i < NO_SEG; i++){
    est_err[i][0] = 0;
    est_err[i][1] = 0;
  }
  // tile version
  tile_ver = init2dArrayInt(NO_SEG, metadata->video_info.tile[TILING].No_tile);
  tile_size = init2dArrayInt(NO_SEG, metadata->video_info.tile[TILING].No_tile);
  decide_width = new int[NO_SEG];
  //
  proc_time = new int[NO_SEG];
  //
  vpsnr = new double[NO_FRAME];
  seg_br = new double[NO_SEG];
  all_seg_br = new double[NO_SEG];
  vp_br = new double[NO_SEG];
  seg_vp_psnr = new double[NO_SEG];
  // 
  VPSNR_thres = new double[NO_SEG];
  //
  visiTileOut = new int[NO_FRAME];
  lowQLPixelPercent = init2dArrayDouble(NO_FRAME, NO_VER);
  //
  ext_tile = new int[NO_SEG];
  useful_ext_tile = new int[NO_FRAME];
  useful_ext_tile_percent = new double[NO_FRAME];
  useful_ext_tile_br = new double[NO_FRAME];
  avg_visi_tile_ver = new double[NO_FRAME];
  avg_ext_tile_ver = new double[NO_FRAME];
  ext_tile_br = new double[NO_FRAME];
  ext_tile_br_useful = new double[NO_FRAME];
  s_ver = new double[NO_FRAME];
  est_err_ang = new double[NO_FRAME];
  est_err_ang_2 = new double[NO_FRAME];
  visiTileOut_percent = new double[NO_FRAME];
  visiTileOut_br = new double[NO_FRAME];
  last_frame_id = new int[NO_FRAME];
  // load real bandwidth traces
  sprintf(buff, "data/bw_trace/trace_%d.txt", adaptInfo.BWTRACE_ID);
  ret = metadata->import_matrix_from_txt_file(buff, v, rows, cols);
  bw_trace = init2dArrayDouble(rows, 2);
  for(i=0; i < rows; i++){
    bw_trace[i][0] = v[i*2];
    bw_trace[i][1] = v[i*2+1];
  }
  printf("#[DecisionEngine] init() finished #\n");
}
DecisionEngine::~DecisionEngine(){
}
int* DecisionEngine::get_next_segment(int index){
  /* calculate processing time */
  struct timeval t_start, t_end;
  int i;
  gettimeofday(&t_start, NULL);
  //
  thrp_estimator(index);
  vp_estimator(index);
  //
  switch(TILE_SELECT_METHOD){
    case 1:
      tile_ver[index] = DASH(index);
      break;
    case 2:
      tile_ver[index] = EQUAL(index); // EQUAL
      break;
    case 3:
      tile_ver[index] = EXT_ALL(index, 0); // all tiles in a area have the same version.
      break;
    case 4:
      tile_ver[index] = ISM(index);
      break;
    case 5:
      tile_ver[index] = ISM_ext(index);
      break;
    case 6:
      tile_ver[index] = Ghent(index);
      break;
    case 7:
      tile_ver[index] = Ireland(index);
      break;
    case 8:
      tile_ver[index] = EXT_ALL_same_ver(index, 0); // ROI
      break;
    case 9:
      tile_ver[index] = ISM_same_ver(index);
      break;
    case 10:
      tile_ver[index] = Ghent_opt(index);
      break;
    case 11:
      tile_ver[index] = BellLab(index);
      break;
    case 12:
      tile_ver[index] = ProbDASH(index);
      break;
    case 13:
      tile_ver[index] = Ireland_v2(index);
      break;
    case 14:
      tile_ver[index] = OPTIMAL(index);
      break;
  }
  /* calculate processing time */
  gettimeofday(&t_end, NULL);
  proc_time[index] = tvdiff_us(&t_end, &t_start);
  seg_br[index] = 0;
  // showTileVersion(tile_ver[index], No_tile_h, No_tile_v);
  for(int i=0; i < No_tile; i++){
    seg_br[index] += TILE_SEG_BR[index][i][tile_ver[index][i]];
    if(TILE_SELECT_METHOD == 1) break;
  }
  return tile_ver[index];
}

double DecisionEngine::calc_seg_down_time(double t_start, double BR, double SD, double** bw_trace){
  int N_0;
  double seg_size = BR * SD; // kbits
  double down_time = 0;
  // calculate N_0
  N_0 = 0;
  while(t_start > bw_trace[N_0][0]) N_0++;
  down_time = seg_size/bw_trace[N_0][1];
  if(down_time <= (bw_trace[N_0][0] - t_start)){
    return down_time;
  }
  seg_size -= (bw_trace[N_0][0] - t_start)*bw_trace[N_0][1];
  down_time =  (bw_trace[N_0][0] - t_start);
  N_0 ++;
  while((seg_size / bw_trace[N_0][1]) > (bw_trace[N_0][0] - bw_trace[N_0-1][0])){
    seg_size -= bw_trace[N_0][1]*(bw_trace[N_0][0] - bw_trace[N_0-1][0]);
    down_time +=  (bw_trace[N_0][0] - bw_trace[N_0-1][0]);    
    N_0 ++;
  }
  down_time += seg_size / bw_trace[N_0][1];
  return down_time;
}
// for real bandwidth trace
void DecisionEngine::down_next_segment_2(int index){
  // request and receive tiles' versions
  // calculate download metrics and update buffer
  // throughput 
  if(index==0)
    seg_down_start_time[index] = 0;
  else
    seg_down_start_time[index] = seg_down_start_time[index-1] + seg_down_time[index-1] + time_to_next_request[index-1]; 
  seg_down_time[index] = calc_seg_down_time(seg_down_start_time[index], seg_br[index],INTER*1.0/FPS, bw_trace);
  seg_thrp[index] = seg_br[index]*1.0*(INTER*1.0/FPS) / seg_down_time[index];
  if(index < BUFF/INTER){
    if(index == 0)
      buff_level[index] =  INTER*1.0/FPS;
    else
      buff_level[index] = buff_level[index-1] + INTER*1.0/FPS;
    stall_time[index] = 0;
    last_frame_id[index] = 0;
  }
  else{
    if(BUFFERING_STATE == false){
      if(buff_level[index-1] < seg_down_time[index])
        stall_time[index] = seg_down_time[index] - buff_level[index-1];
      else
        stall_time[index] = 0;
      buff_level[index] = buff_level[index-1] + INTER*1.0/FPS - seg_down_time[index] + stall_time[index];
      if(stall_time[index] > 0)
        BUFFERING_STATE = true;
      // update last played frame
      last_frame_id[index] = last_frame_id[index-1]+ (int)(min(seg_down_time[index], buff_level[index-1]) * FPS);
    }else{
      stall_time[index] = seg_down_time[index];
      buff_level[index] += INTER*1.0/FPS;
    }
  } 
  // switch to normal state after a rebuffering
  if(BUFFERING_STATE == true && buff_level[index] == B_max)
    BUFFERING_STATE = false; 
  time_to_next_request[index] = 0;
  // Ghent method
  if(TILE_SELECT_METHOD == 6 || TILE_SELECT_METHOD == 10){
    if(BUFFERING_STATE == false){
      if(buff_level[index] >= B_cri){
        last_frame_id[index] += (int)((buff_level[index] - B_cri) * FPS); 
        buff_level[index] = B_cri;
        time_to_next_request[index] += (buff_level[index] - B_cri);
      }
    }
  }
}
void DecisionEngine::down_next_segment(int index){
  // request and receive tiles' versions
  // calculate download metrics and update buffer
  // throughput 
  seg_down_time[index] = seg_br[index]*1.0*(INTER*1.0/FPS) / seg_thrp[index];
  if(index < BUFF/INTER){
    buff_level[index] += INTER*1.0/FPS;
    stall_time[index] = 0;
    last_frame_id[index] = 0;
  }
  else{
    if(BUFFERING_STATE == false){
      if(buff_level[index-1] < seg_down_time[index])
        stall_time[index] = seg_down_time[index] - buff_level[index-1];
      else
        stall_time[index] = 0;
      buff_level[index] = buff_level[index-1] + INTER*1.0/FPS - seg_down_time[index] + stall_time[index];
      if(stall_time[index] > 0)
        BUFFERING_STATE = true;
      // update last played frame
      last_frame_id[index] = last_frame_id[index-1] + (int)(min(seg_down_time[index], buff_level[index-1]) * FPS);
    }else{
      stall_time[index] = seg_down_time[index];
      buff_level[index] += INTER*1.0/FPS;
    }
  } 
  // switch to normal state after a rebuffering
  if(BUFFERING_STATE == true && buff_level[index] == B_max)
   BUFFERING_STATE = false; 
  // Ghent method
  if(TILE_SELECT_METHOD == 6 || TILE_SELECT_METHOD == 10){
    if(BUFFERING_STATE == false){
      if(buff_level[index] >= B_cri){
        last_frame_id[index] += (int)((buff_level[index] - B_cri) * FPS); 
        buff_level[index] = B_cri;
      }
    }
  }
}
void DecisionEngine::calc_result(int NO_SEG){
  int index, i, j, num1, num2, num;
  int* vpixel;
  double* vp_size = new double[INTER];
  // int* required_ext_tile, *ext_tile, *wasted_tile;
  for(index = 0; index < NO_SEG; index ++){
    vp_br[index] = 0;
    seg_vp_psnr[index] = 0;
    /* calculate vpsnr and segment bitrate */
    if(TILE_SELECT_METHOD == 1){ // DASH
      seg_br[index] = DASH_SEG_BR[index][tile_ver[index][0]];
      for(i=0; i < INTER; i++){
        if(index >= 2)
          vpsnr[index * INTER + i] = est_vp_psnr(DASH_TILE_MSE[index], No_tile, tile_ver[index], htrace[index * INTER + i]);
        else
          vpsnr[index * INTER + i] = DASH_SEG_PSNR[index][tile_ver[index][0]];
      }
      seg_vp_psnr[index] = avg(vpsnr+index*INTER, INTER);
    }else{  // tiling-based
      // total bitrate
      seg_br[index] = 0;
      for(i=0; i < No_tile; i++)
        seg_br[index] += TILE_SEG_BR[index][i][tile_ver[index][i]];
      // viewport PSNR
      for(i=0; i < INTER; i++){
        vpsnr[index * INTER + i] = est_vp_psnr(TILE_SEG_MSE[index], No_tile, tile_ver[index], htrace[index * INTER + i]);
        seg_vp_psnr[index] += vpsnr[index * INTER + i];
      }
      seg_vp_psnr[index] /= (1.0 * INTER);
      // avgVPSNR = avg(vpsnr+2*INTER, NO_FRAME - 2*INTER);
        
      // Viewport-bitrate
      // show tiles' bitrate
      // printf("#[calc_result] index=%d\n", index);
      // for(i=0; i < No_tile; i++){
      //  printf("%.2f ", TILE_SEG_BR[index][i][tile_ver[index][i]]);
      //  if((i+1) % No_tile_h == 0)
      //    printf("\n");
      // }
      // printf("#[calc_result]\n");
      for(i=0; i < INTER; i++){
        vpixel = get_visible_pixel(htrace[index * INTER + i]);
        vp_size[i] = 0;
        for(j=0; j < No_tile; j++){
          // printf("%.2f ", vpixel[j] * 1.0 / (tile_W * tile_H));
          vp_size[i] += TILE_SEG_SIZE[index][j][tile_ver[index][j]]*1.0/INTER * (vpixel[j] * 1.0 / (tile_W * tile_H));
          // if((j+1)%No_tile_h == 0)
            // printf("\n");
        }
        vp_br[index] += vp_size[i]; 
      }
      vp_br[index] /= (1.0 * INTER / FPS);
    }
    // calculate number of visible tiles 
    for(i=0; i < INTER; i++){
      // calculate number of required extension tiles
      visiTileOut[index*INTER + i] = calc_visiTileOut(get_visible_tile(est_vp[index]), get_visible_tile(htrace[index*INTER+i]));
      visiTileOut_percent[index*INTER + i] = calc_visiTileOut_percent(index, i);
      visiTileOut_br[index*INTER + i] = calc_visiTileOut_br(index, i);
      avg_visi_tile_ver[index*INTER + i] = calc_avg_visi_tile_ver(index, i);
      avg_ext_tile_ver[index*INTER + i] = calc_ext_tile_avg_ver(index);
      ext_tile_br[index*INTER + i] = calc_ext_tile_br(index);
      ext_tile_br_useful[index*INTER + i] = calc_ext_tile_br_useful(index, i);
      s_ver[index*INTER + i] = calc_sversion(index, i);
      useful_ext_tile[index*INTER + i] = calc_useful_ext_tile(index, i);
      useful_ext_tile_percent[index*INTER + i] = calc_useful_ext_tile_percent(index, i);
      useful_ext_tile_br[index*INTER + i] = calc_useful_ext_tile_br(index, i);
      // 
      for(j=0; j < NO_VER; j++)
        lowQLPixelPercent[index*INTER + i][j] = calc_lowQLPixelPercent(get_visible_pixel(est_vp[index]), get_visible_pixel(htrace[index*INTER+i]), tile_ver[index], j);

    }
    ext_tile[index] = calc_ext_tile(index);
    if(index == 5){
        showTileVersion(tile_ver[index], No_tile_h, No_tile_v);
        printf("\n");
        for(i=0; i < INTER; i++){
          printf("\n");
          showTileVersion(get_visible_pixel(htrace[index*INTER + i]), No_tile_h, No_tile_v);
        }
    }
    // if(index >= 3){
    //   required_ext_tile = calc_required_ext_tile(index, &num1);
    //   ext_tile = calc_ext_tile(index, &num2);
    //   wasted_tile = calc_wasted_ext_tile(required_ext_tile, num1, ext_tile, num2, &num);
    // }
  }
  avgVPSNR = avg(vpsnr+2*INTER, NO_SEG*INTER - 2*INTER);
  //}
}
int DecisionEngine::calc_visiTileOut(int* vmask, int* vmask2){
  int ret = 0, i;
  for(i=0; i < No_tile; i++){
    if(vmask2[i] == 1 && vmask[i] == 0)
      ret++;
  }
  return ret;
}
double DecisionEngine::calc_visiTileOut_percent(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  int* vpixel_frame;
  double pecent = 0;
  vmask_frame = get_visible_tile(htrace[index*INTER + frame_id]);
  vpixel_frame = get_visible_pixel(htrace[index*INTER + frame_id]);
  for(i=0; i < No_tile; i++){
    if(vmask_frame[i] == 1 && vmask_est[i] == 0)
      pecent += vpixel_frame[i];
  }
  return (pecent * 100.0 / (vp_W * vp_H));
}
double DecisionEngine::calc_visiTileOut_br(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  int* vpixel_frame;
  double pecent = 0;
  vmask_frame = get_visible_tile(htrace[index*INTER + frame_id]);
  vpixel_frame = get_visible_pixel(htrace[index*INTER + frame_id]);
  for(i=0; i < No_tile; i++){
    if(vmask_frame[i] == 1 && vmask_est[i] == 0)
      pecent += TILE_SEG_BR[index][i][tile_ver[index][i]] - TILE_SEG_BR[index][i][0];
  }
  if(calc_ext_tile(index) == 0)
    return 0;
  return pecent;
}
double DecisionEngine::calc_lowQLPixelPercent(int* vpixel, int* vpixel2, int* tileVer, int VER){
  int i, pixSum = 0;
  for(i=0; i < No_tile; i++){
    // if(vpixel2[i] > 0 && vpixel[i] == 0 && tileVer[i] == VER)
    if(vpixel2[i] > 0 && tileVer[i] == VER)
      pixSum += vpixel2[i];
  }
  return 100 * pixSum * 1.0 / (vp_W * vp_H);
}
int* DecisionEngine::calc_required_ext_tile(int index, int* num){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  bool FLAG;
  std::vector<int> required_ext_tile;
  for(i=0; i < INTER; i++){
    // printf("#i=%d:\n", (i+1));
    vmask_frame = get_visible_tile(htrace[index*INTER +i]);
    // showTileVersion(vmask_frame, No_tile_h, No_tile_v);
    for(j=0; j < No_tile; j++){
      if(vmask_frame[j] == 1 && vmask_est[j] == 0){
        FLAG = true;
        for(jj = 0; jj < required_ext_tile.size(); jj++){
          if(required_ext_tile[jj] == j){
            FLAG = false;
            break;
          }
        }
        if(FLAG){
          required_ext_tile.push_back(j);
        }
      }
    }
  }
  *num = required_ext_tile.size();
  return required_ext_tile.data();
}
int DecisionEngine::calc_ext_tile(int index){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  bool FLAG;
  int num = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0)
      num++;
  if(TILE_SELECT_METHOD == 1 || TILE_SELECT_METHOD == 2 || TILE_SELECT_METHOD == 3)
    num = 0;
  if(TILE_SELECT_METHOD == 4 && decide_width[index] == 0){
    num = 0;
  }
  if(index==5 && TILE_SELECT_METHOD == 11){
    vmask_est = get_visible_tile(est_vp[index]);
    showTileVersion(vmask_est, No_tile_h, No_tile_v);
  }
  return num;
}
double DecisionEngine::calc_ext_tile_avg_ver(int index){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  bool FLAG;
  int num = 0;
  double avg_ver = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0){
      avg_ver += tile_ver[index][i];
      num++;
    }
  if(TILE_SELECT_METHOD == 1 || TILE_SELECT_METHOD == 2 || TILE_SELECT_METHOD == 3)
    num = 0;
  if(TILE_SELECT_METHOD == 4 && decide_width[index] == 0){
    num = 0;
  }
  if(num == 0)
    return 0;
  else
    return (avg_ver/num);
}
double DecisionEngine::calc_ext_tile_br(int index){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame;
  bool FLAG;
  double num = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0)
      num += (TILE_SEG_BR[index][i][tile_ver[index][i]] - TILE_SEG_BR[index][i][0]);
  if(TILE_SELECT_METHOD == 1 || TILE_SELECT_METHOD == 2 || TILE_SELECT_METHOD == 3)
    num = 0;
  if(TILE_SELECT_METHOD == 4 && decide_width[index] == 0){
    num = 0;
  }
  return num;
}
double DecisionEngine::calc_ext_tile_br_useful(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame = get_visible_tile(htrace[index*INTER + frame_id]);
  bool FLAG;
  double num = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0 && vmask_frame[i] == 1)
      num += (TILE_SEG_BR[index][i][tile_ver[index][i]] - TILE_SEG_BR[index][i][0]);
  if(TILE_SELECT_METHOD == 1 || TILE_SELECT_METHOD == 2 || TILE_SELECT_METHOD == 3)
    num = 0;
  if(TILE_SELECT_METHOD == 4 && decide_width[index] == 0){
    num = 0;
  }
  return num;
}
int DecisionEngine::calc_useful_ext_tile(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame = get_visible_tile(htrace[index*INTER+frame_id]);
  bool FLAG;
  int num = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0 && vmask_frame[i] == 1)
      num++;
  return num;
}
double DecisionEngine::calc_useful_ext_tile_percent(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame = get_visible_tile(htrace[index*INTER+frame_id]);
  int* vpixel = get_visible_pixel(htrace[index*INTER + frame_id]);
  bool FLAG;
  double num = 0;
  
}

double DecisionEngine::calc_useful_ext_tile_br(int index, int frame_id){
  int i,j, jj;
  int* vmask_est = get_visible_tile(est_vp[index]);
  int* vmask_frame = get_visible_tile(htrace[index*INTER+frame_id]);
  int* vpixel = get_visible_pixel(htrace[index*INTER + frame_id]);
  bool FLAG;
  double num = 0;
  for(i=0; i < No_tile; i++)
    if(tile_ver[index][i] > 0 && vmask_est[i] == 0 && vmask_frame[i] == 1)
      num += (vpixel[i]*1.0/(vp_W * vp_H)) * (TILE_SEG_BR[index][i][tile_ver[index][i]] - TILE_SEG_BR[index][i][0]);
  return num;
}
double DecisionEngine::calc_avg_visi_tile_ver(int index, int frame_id){
  int i,j,cnt=0;
  int* vmask = get_visible_tile(htrace[index*INTER + frame_id]);
  double avg_ver = 0;
  for(i=0; i < No_tile; i++){
    if(vmask[i] == 1){
      avg_ver += tile_ver[index][i];
      cnt++;
    }
  }
  return avg_ver/cnt;
}

double DecisionEngine::calc_sversion(int index, int frame_id){
  int i,j,cnt=0;
  int* vmask = get_visible_tile(htrace[index*INTER + frame_id]);
  int* vpixel = get_visible_pixel(htrace[index*INTER+frame_id]);
  double avg_ver = 0;
  for(i=0; i < No_tile; i++){
    if(vmask[i] == 1){
      avg_ver += tile_ver[index][i] * vpixel[i] * 1.0 / (vp_W * vp_H);
      cnt++;
    }
  }
  return avg_ver;
}
char* DecisionEngine::get_method_name(int method_id){
  char* full_name = new char[1024];
  char* name[100] = {"DASH","EQUAL","ROI", "ISM", "Ghent", "Ireland"};
  sprintf(full_name, "%df_%dx%d_%s", No_face, No_tile_h, No_tile_v, name[method_id-1]);
  return full_name;
}
void DecisionEngine::thrp_estimator(int index){
  double alpha = 0.8;
  double margin = 0.0;
  if(index == 0){
    est_thrp_act[index] = 0;
    est_seg_thrp[index] = 0;
  }else{
    est_thrp_act[index] = (1- alpha) * est_thrp_act[index -1] + alpha * seg_thrp[index-1];
    est_seg_thrp[index] = (1-margin) * est_thrp_act[index];
    est_seg_thrp[index] = seg_thrp[index-1]; // last throughput-based
  }
  // decide the total bitrate allocated for tiles
  switch(TILE_SELECT_METHOD){
    case 12:
      if(index==0)
        all_seg_br[index] = 0;
      else
        if(BUFFERING_STATE == false)
          all_seg_br[index] = (buff_level[index-1] - B_target + INTER * 1.0/FPS)*est_seg_thrp[index]/(INTER * 1.0 / FPS); 
        else
          all_seg_br[index] = est_seg_thrp[index]/2; 
      break; 
  }
}
void DecisionEngine::vp_estimator(int index){
  int i, VP_EST_WIN;
  double v_phi;
  double v_theta;
  int delta_phi, delta_theta;;
  VP_EST_WIN = INTER;
  net_delay[index] = 50; // ms --> parameterize needed!
  // printf("#[vp_estimator] index=%d\n", index);
  switch(VP_EST_METHOD){
    case 0: // know all
      // est_err[index][0] = 0;
      // est_err[index][1] = 0;  
      est_vp[index][0] = htrace[index * INTER][0];
      est_vp[index][1] = htrace[index * INTER][1];
      for(i=0; i < INTER; i++){
        est_frame_vp[index * INTER + i][0] = htrace[index * INTER + i][0];
        est_frame_vp[index * INTER + i][1] = htrace[index * INTER + i][1];
      }
      cur_vp[index] = est_vp[index];
       // calculate move speed
        // phi
        delta_phi = htrace[index * INTER + INTER - 1][0] - htrace[index * INTER][0];
        if(delta_phi < -180)
          delta_phi += 360;
        else
          if(delta_phi > 180)
            delta_phi -= 360;
        speed[index][0] = delta_phi/(1.0 * INTER);
        // theta
        delta_theta = htrace[index * INTER + INTER - 1][1] - htrace[index * INTER][1];
        // printf("delta_theta=%d\n", delta_theta);
        if(delta_theta < -90)
          delta_theta += 180;
        else
          if(delta_theta > 90)
            delta_theta -= 180;
        speed[index][1] = delta_theta/(1.0 * INTER);
        ang_speed[index] = sqrt(speed[index][0] * speed[index][0] + speed[index][1] * speed[index][1]);
      break;
    case 1:// linear regression
      if(index < BUFF/INTER){
        for(i=0; i < INTER; i++){
          est_frame_vp[index * INTER + i][0] = 0;
          est_frame_vp[index * INTER + i][1] = 0;
        }
        est_err[index][0] = 0;
        est_err[index][1] = 0;
        //
        cur_vp[index][0] = 0;
        cur_vp[index][1] = 0;     
        //
        est_vp[index][0] = 0;
        est_vp[index][1] = 0;
        // angular speed
        ang_speed[index] = 0;
      }else{
        //last_frame_id = index * INTER - BUFF;
        cur_vp[index][0] = htrace[last_frame_id[index-1]][0];
        cur_vp[index][1] = htrace[last_frame_id[index-1]][1];
        // calculate estimation errors 
        // phi
        if(last_frame_id[index-1] < BUFF){
          est_err[index][0] = 0;
          est_err[index][1] = 0;
        }else{
          est_err[index][0] = -est_frame_vp[last_frame_id[index-1]][0] + htrace[last_frame_id[index-1]][0];
          if(est_err[index][0] < -180)
            est_err[index][0] += 360;
          else if(est_err[index][0] > 180)
            est_err[index][0] -= 360;
          // theta
          est_err[index][1] = - est_frame_vp[last_frame_id[index-1]][1] + htrace[last_frame_id[index-1]][1];
          if(est_err[index][1] < -90)
            est_err[index][1] += 180;
          else if(est_err[index][1] > 90)
            est_err[index][1] -= 180;
        } 
        //
        // calculate move speed
        // phi
        if(last_frame_id[index-1] >= VP_EST_WIN){
          delta_phi = htrace[last_frame_id[index-1]][0] - htrace[last_frame_id[index-1]-VP_EST_WIN][0];
        }
        else{
          delta_phi = htrace[last_frame_id[index-1]][0] - htrace[0][0];
        }
        if(delta_phi < -180)
          delta_phi += 360;
        else
          if(delta_phi > 180)
            delta_phi -= 360;
        if(last_frame_id[index-1] == 0)
          speed[index][0] = 0;
        else
          speed[index][0] = delta_phi/(1.0 * ((last_frame_id[index-1] >= VP_EST_WIN)?VP_EST_WIN:last_frame_id[index-1]));
        // theta
        if(last_frame_id[index-1] >= VP_EST_WIN)
          delta_theta = htrace[last_frame_id[index-1]][1] - htrace[last_frame_id[index-1]-VP_EST_WIN][1];
        else
          delta_theta = htrace[last_frame_id[index-1]][1] - htrace[0][1];
        // printf("delta_theta=%d\n", delta_theta);
        if(delta_theta < -90)
          delta_theta += 180;
        else
          if(delta_theta > 90)
            delta_theta -= 180;
        if(last_frame_id[index-1] == 0)
          speed[index][1] = 0;
        else
          speed[index][1] = delta_theta/(1.0 * ((last_frame_id[index-1] >= VP_EST_WIN)?VP_EST_WIN:last_frame_id[index-1]));

        // estimate viewport 
        for(i=0; i < INTER; i++){
          est_frame_vp[index * INTER + i][0] = (int) (cur_vp[index][0] + speed[index][0] * (buff_level[index-1]*FPS + i));
          est_frame_vp[index * INTER + i][1] = (int) (cur_vp[index][1] + speed[index][1] * (buff_level[index-1]*FPS + i));
          // phi
          while(est_frame_vp[index * INTER + i][0] >= 180)
            est_frame_vp[index * INTER + i][0] -= 360;
          while(est_frame_vp[index * INTER + i][0] < -180)
            est_frame_vp[index * INTER + i][0] += 360;
          // theta
          while(est_frame_vp[index * INTER + i][1] >= 90)
            est_frame_vp[index * INTER + i][1] -= 90;
          while(est_frame_vp[index * INTER + i][1] <= -90)//
            est_frame_vp[index * INTER + i][1] += 90;
          // printf("(%d,%d) - (%d, %d)\n", est_frame_vp[index * INTER + i][0], est_frame_vp[index * INTER + i][1], htrace[index * INTER + i][0],htrace[index * INTER + i][1]);
          est_err_ang[index*INTER+i] = acos(sin(htrace[index*INTER+i][1]*M_PI/180) * sin(est_frame_vp[index*INTER+i][1]*M_PI/180) + cos(htrace[index*INTER+i][1]*M_PI/180) * cos(est_frame_vp[index*INTER+i][1]*M_PI/180) * cos(abs(est_frame_vp[index*INTER+i][0] - htrace[index*INTER+i][0])*M_PI/180)) / M_PI * 180;
            //printf("(%d,%d) - (%d, %d) %.2f\n", est_frame_vp[index * INTER + i][0], est_frame_vp[index * INTER + i][1], htrace[index * INTER + i][0],htrace[index * INTER + i][1], est_err_ang[index*INTER+i]);
          est_frame_vp_2[index*INTER + i][0] = (int)(est_frame_vp[index*INTER+i][0] + est_err[index][0] * i * 1.0 / (INTER-1));
          est_frame_vp_2[index*INTER + i][1] = (int)(est_frame_vp[index*INTER+i][1] + est_err[index][1] * i * 1.0 / (INTER-1));
          // phi
          while(est_frame_vp_2[index * INTER + i][0] >= 180)
            est_frame_vp_2[index * INTER + i][0] -= 360;
          while(est_frame_vp_2[index * INTER + i][0] < -180)
            est_frame_vp_2[index * INTER + i][0] += 360;
          // theta
          while(est_frame_vp_2[index * INTER + i][1] >= 90)
            est_frame_vp_2[index * INTER + i][1] -= 90;
          while(est_frame_vp_2[index * INTER + i][1] <= -90)//
            est_frame_vp_2[index * INTER + i][1] += 90;
          est_err_ang_2[index*INTER+i] = acos(sin(htrace[index*INTER+i][1]*M_PI/180) * sin(est_frame_vp_2[index*INTER+i][1]*M_PI/180) + cos(htrace[index*INTER+i][1]*M_PI/180) * cos(est_frame_vp_2[index*INTER+i][1]*M_PI/180) * cos(abs(est_frame_vp_2[index*INTER+i][0] - htrace[index*INTER+i][0])*M_PI/180)) / M_PI * 180;
        }
        est_vp[index] = est_frame_vp[index * INTER];
        ang_speed[index] = sqrt(speed[index][0] * speed[index][0] + speed[index][1] * speed[index][1]);
        // printf("#[vp_estimator] est_vp[%d] = (%d, %d)\n", index, est_vp[index][0], est_vp[index][1]);
      }
      break;
  }
  for(i=0; i < INTER; i++){
   est_err_frame[index*INTER + i][0] = htrace[index*INTER + i][0] - est_frame_vp[index*INTER + i][0];
   if(est_err_frame[index*INTER + i][0] > 180)
     est_err_frame[index*INTER + i][0] -= 360;
   if(est_err_frame[index*INTER + i][0] < -180)
     est_err_frame[index*INTER + i][0] += 360;
   // 
   est_err_frame[index*INTER + i][1] = htrace[index*INTER + i][1] - est_frame_vp[index*INTER + i][1];
  }
  //
  est_vp[index] = est_frame_vp[index * INTER];
}
int* DecisionEngine::DASH(int index){
  int* tile = new int[No_tile];
  int selectVer;
  selectVer = 0;
  // printf("#[DASH] index=%d\n est_thrp=%.2f\n", index, est_seg_thrp[index]);
  for(int i=0; i < NO_VER; i++)
    // printf("#[DASH] BR[%d][%d] = %.2f\n", index, i, DASH_SEG_BR[index][i]);
    while(selectVer < NO_VER && DASH_SEG_BR[index][selectVer] < est_seg_thrp[index]) selectVer ++;
  //
  selectVer = (selectVer > 0)?(selectVer-1):selectVer;
  // printf("#[DASH] selected version: %d\n", selectVer);
  //
  for(int i=0; i < No_tile; i++){
    tile[i] = selectVer;
  }
  return tile;
}
int* DecisionEngine::EXT_ALL(int index, int ext_width){
  int* vmask = NULL;
  decide_width[index] = ext_width;
  if(index >=  BUFF / INTER){
    vmask = get_visible_tile(est_vp[index]);
    vmask = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_width);
    if(index==5 && ext_width == 0){
      showTileVersion(vmask, No_tile_h, No_tile_v);
    }
  }
  return ROI(TILE_SEG_BR[index], No_tile, NO_VER, vmask, est_seg_thrp[index]);
}
int* DecisionEngine::EXT_ALL_same_ver(int index, int ext_width){
  int* vmask = NULL;
  decide_width[index] = ext_width;
  if(index >  BUFF / INTER){
    vmask = get_visible_tile(est_vp[index]);
    vmask = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_width);
  }
  return ROI_same_ver(TILE_SEG_BR[index], No_tile, NO_VER, vmask, est_seg_thrp[index]);
}
int* DecisionEngine::EQUAL(int index){
  int* tileVer = new int[No_tile];
  int i,j;
  double totalBR = est_seg_thrp[index];
  double tmpSum;
  decide_width[index] = 0;
  i=0;
  while(true){
    tmpSum = 0;
    for(j=0; j < No_tile; j++)
      tmpSum += TILE_SEG_BR[index][j][i];
    if(tmpSum > totalBR || i == NO_VER -1) break;
    i++;
  }
  i = (i==0)?0:(i-1);
  for(j=0; j < No_tile; j++)
    tileVer[j] = i;
  return tileVer;
}
/*
* - Assign lowest version for invisible tiles
* - Assign highest possible (same) version for visible tiles
* - Using remaining bandwidth to improve tile quality
*   + calculate amount of BW required to increase the version a tile to the next level
*   + starting from the tile with the lowest required BW, increase the version to the next higher version,  till all BW is used.  
*
*/
int* DecisionEngine::ROI(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp){
  int i,j;
  int select_ver;
  int* tile = new int[NO_TILE];
  double sum;
  double delta_R;
  typedef std::pair<int,double> pair;
  std::map<int,double> candiTile;
  std::vector<pair> vec;
  if(vmask == NULL){ // should be changed to EQUAL later
    for(i=0; i < NO_TILE; i++){
      tile[i] = 0;
    }
  }else{
    // allocate the minimum version for visible tiles
    sum = 0;
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] == 0){
        tile[i] = 0;
        sum += TILE_BR[i][tile[i]];
      }
    }
    // showTileVersion(tile, No_tile_h, No_tile_v);
    delta_R = est_thrp - sum;
    // assign highest possible version for visible tiles
    for(j=0; j < NO_VER; j++){
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          sum += TILE_BR[i][j];
        }
      }
      // printf("j=%d sum=%.2f delta_R=%.2f\n", j, sum, delta_R);
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
    // showTileVersion(tile, No_tile_h, No_tile_v);
    //
    bool FLAG = true;
    /* enhance quality of visible tiles */
    //while(delta_R > 0){
      //FLAG = true;
      // search for candidate tiles
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 1)
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R){
            candiTile[i] = TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]];
          }
      }
      // sort candidate tiles according to the required additional bandwidth
      if(candiTile.size() > 0){
        std::copy(candiTile.begin(), candiTile.end(), std::back_inserter<std::vector<pair>>(vec));
        std::sort(vec.begin(), vec.end(),
            [](const pair& l, const pair& r) {
            if (l.second != r.second)
            return l.second < r.second;

            return l.first < r.first;
            });
        for (auto const &pair: vec) {
          std::cout << '{' << pair.first << "," << pair.second << '}' << '\n';
          if(pair.second < delta_R){
            delta_R -= pair.second;
            tile[pair.first] += 1;
            FLAG = false;
          }
        }
      }
      /*
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] == 1){
          if(tile[i] < NO_VER -1 && (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]) <= delta_R){
            delta_R -= (TILE_BR[i][tile[i] + 1] - TILE_BR[i][tile[i]]);
            tile[i] += 1;
            FLAG = false;
          }
        }
      }
      */
      //if(FLAG) break; // break if no improvement is done in the last round.
    //}
    /* enhance quality of invisible tiles */
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
    // showTileVersion(tile, No_tile_h, No_tile_v);
  }
  return tile;
}
/*
* - Assign lowest version for invisible tiles
* - Assign highest possible (same) version for visible tiles
* - NOT use remaining BW to improve tiles' versions (avoid bias when comparing to other methods)
*/
int* DecisionEngine::ROI_same_ver(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp){
  int i,j;
  int select_ver;
  int* tile = new int[NO_TILE];
  double sum;
  double delta_R;
  typedef std::pair<int,double> pair;
  std::map<int,double> candiTile;
  std::vector<pair> vec;
  if(vmask == NULL){ // should be changed to EQUAL later
    for(i=0; i < NO_TILE; i++){
      tile[i] = 0;
    }
  }else{
    // allocate the lowest version for every tile
    sum = 0;
    for(i=0; i < NO_TILE; i++){
        tile[i] = 0;
        sum += TILE_BR[i][tile[i]];
    }
    // showTileVersion(tile, No_tile_h, No_tile_v);
    delta_R = est_thrp - sum;
    if(delta_R <= 0)
      return tile;
    // assign highest possible version for visible tiles
    for(j=1; j < NO_VER; j++){
      sum = 0;
      for(i=0; i < NO_TILE; i++){
        if(vmask[i] > 0){
          sum += (TILE_BR[i][j] - TILE_BR[i][0]);
        }
      }
      // printf("j=%d sum=%.2f delta_R=%.2f\n", j, sum, delta_R);
      if(sum > delta_R)
        break;
    }
    for(i=0; i < NO_TILE; i++){
      if(vmask[i] > 0){
        tile[i] = j - 1;
      }
    }
  }
  return tile;
}
int* DecisionEngine::ISM(int index){
  int width;
  int**tileVer = new int*[3];
  int* selectTileVer = NULL;
  double max_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  int INTERVAL = INTER;
  int NO_TILE = No_tile;
  if(index <= BUFF/INTER){
    return EXT_ALL(index, 0);
  }
  for(width = 0; width <= 2; width ++){
    // printf("[#ISM] index=%d with=%d est_vp: (%d, %d)\n", index, width, est_vp[index][0], est_vp[index][1]);
    if(width == 0){
      tileVer[width] = EXT_ALL(index, 0); // ROI
    }else{
      tileVer[width] = ISM(index, width);
    }
    // showTileVersion(tileVer[width], No_tile_h, No_tile_v);
    tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
    tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr = est_vp_psnr(TILE_SEG_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
    tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr_last = est_vp_psnr(TILE_SEG_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    //
    // printf("[#ISM] index=%d width=%d: vppsnr=%.2f\n",index, width, (vp_psnr +vp_psnr_last)/2);
    if((vp_psnr + vp_psnr_last) > max_psnr){
      selectTileVer = tileVer[width];
      max_psnr = vp_psnr + vp_psnr_last;
      decide_width[index] = width;
    }
  }
  // printf("index=%d select: max_psnr=%.2f w=%d\n", index, max_psnr, decide_width[index]);
  return selectTileVer;
}
int* DecisionEngine::ISM(int index, int ext_width){
  int i, j, ii, jj;
  int tiles = No_tile;
  int *tileVer = NULL;
  int mean_ver = NO_VER;
  int mean_ver_2 = NO_VER;
  int vp_ver;
  int ext_1_ver;
  double BW = est_seg_thrp[index];
  double remainBW;
  bool FLAG;
  int *selectTileVer = new int[tiles];
  double max_vp_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  //
  int width = 0;
  int INTERVAL = INTER;
  double** TILE_MSE = TILE_SEG_MSE[index];
  double** TILE_BR = TILE_SEG_BR[index];
  int NO_TILE = No_tile;
  // int NO_VER = NO_VER;
  //
  // printf("#[ISM][seg #%d]: ext_width:%d BW=%.2f\n", index, ext_width, BW);
  //
  int *vmask = NULL;
  if(index >=  BUFF / INTER){
    vmask = get_visible_tile(est_vp[index]);
    vmask = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_width);
  }
  // showTileVersion(vmask, No_tile_h, No_tile_v);
  // get the selected version when all tiles are of equal quality
  tileVer = ROI(TILE_SEG_BR[index], NO_TILE, NO_VER, vmask, est_seg_thrp[index]);
  // showTileVersion(tileVer, No_tile_h, No_tile_v);
  if(vmask == NULL)
    return tileVer;
  // calculate viewport-psnr given 'tileVer'
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  // first frame & last frame
  tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
  tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
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
  tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
  tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
  // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
  //
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
        tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
        tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
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
        tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
        tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
        // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
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
      tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
      tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
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
      tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
      // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
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
      tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
      tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
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
      tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
      // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
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
        tmp_vp[0] = est_frame_vp[index * INTERVAL][0];
        tmp_vp[1] = est_frame_vp[index * INTERVAL][1];
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
        tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
        tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
        // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
        //
        ext_1_ver ++;
      }
      vp_ver++;
    }
  }
  return selectTileVer;
}
int* DecisionEngine::ISM_same_ver(int index){
  int width;
  int**tileVer = new int*[3];
  int* selectTileVer = NULL;
  double max_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  int INTERVAL = INTER;
  int NO_TILE = No_tile;
  if(index <= BUFF/INTER){
    return EXT_ALL(index, 0);
  }
  for(width = 0; width <= 2; width ++){
    // printf("[#ISM] index=%d with=%d est_vp: (%d, %d)\n", index, width, est_vp[index][0], est_vp[index][1]);
    if(width == 0){
      tileVer[width] = EXT_ALL_same_ver(index, 0); // ROI
    }else{
      tileVer[width] = ISM_same_ver(index, width);
    }
    // showTileVersion(tileVer[width], No_tile_h, No_tile_v);
    tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
    tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr = est_vp_psnr(TILE_SEG_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
    tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
    if(tmp_vp[0] >= 180)
      tmp_vp[0] -= 360;
    if(tmp_vp[0] < -180)
      tmp_vp[0] += 360;
    //
    if(tmp_vp[1] >= 90)
      tmp_vp[1] -= 180;
    if(tmp_vp[1] <= -90)//
      tmp_vp[1] += 180;
    vp_psnr_last = est_vp_psnr(TILE_SEG_MSE[index], NO_TILE, tileVer[width], tmp_vp);
    //
    // printf("[#ISM] index=%d width=%d: vppsnr=%.2f\n",index, width, (vp_psnr +vp_psnr_last)/2);
    if((vp_psnr + vp_psnr_last) > max_psnr){
      selectTileVer = tileVer[width];
      max_psnr = vp_psnr + vp_psnr_last;
      decide_width[index] = width;
    }
  }
  // printf("index=%d select: max_psnr=%.2f w=%d\n", index, max_psnr, decide_width[index]);
  return selectTileVer;
}
int* DecisionEngine::ISM_same_ver(int index, int ext_width){
  int i, j, ii, jj;
  int tiles = No_tile;
  int *tileVer = NULL;
  int mean_ver = NO_VER;
  int mean_ver_2 = NO_VER;
  int vp_ver;
  int ext_1_ver;
  double BW = est_seg_thrp[index];
  double remainBW;
  bool FLAG;
  int *selectTileVer = new int[tiles];
  double max_vp_psnr = 0;
  double vp_psnr;
  double vp_psnr_last;
  //
  int width = 0;
  int INTERVAL = INTER;
  double** TILE_MSE = TILE_SEG_MSE[index];
  double** TILE_BR = TILE_SEG_BR[index];
  int NO_TILE = No_tile;
  double tmp_sum;
  // int NO_VER = NO_VER;
  //
  // printf("#[ISM][seg #%d]: ext_width:%d BW=%.2f\n", index, ext_width, BW);
  //
  int *vmask = NULL;
  if(index >=  BUFF / INTER){
    vmask = get_visible_tile(est_vp[index]);
    vmask = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_width);
  }
  // showTileVersion(vmask, No_tile_h, No_tile_v);
  // get the selected version when all tiles are of equal quality
  tileVer = ROI_same_ver(TILE_SEG_BR[index], NO_TILE, NO_VER, vmask, est_seg_thrp[index]);
  // showTileVersion(tileVer, No_tile_h, No_tile_v);
  if(vmask == NULL)
    return tileVer;
  // calculate viewport-psnr given 'tileVer'
  int* tmp_vp = (int*) malloc(2 * sizeof(int));
  // first frame & last frame
  tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
  tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
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
  tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
  tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
  // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
  // CASE 1: EXT_WIDTH=1
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
      if(remainBW < 0)
        break;
      for(j=1; j < NO_VER; j++){
          tmp_sum = 0;
          for(i=0; i < No_tile; i++){
            if(vmask[i] == 2)
              tmp_sum += (TILE_BR[i][j] - TILE_BR[i][tileVer[i]]);  
          }
          if(tmp_sum > remainBW) break;
        }
      for(i=0; i < No_tile; i++)
          if(vmask[i] == 2)
            tileVer[i] = j-1;
      //
      // if(remainBW <= 0){
      //   // reduce vp tiles's version to meet the bandwidth constraint
      //   for(i=0; i < tiles; i++){
      //     if(vmask[i] == 1){
      //       tileVer[i] -= 1; // assume that all tiles have same priority
      //       remainBW += TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
      //       if(remainBW > 0)
      //         break;
      //     }
      //   }
      //   // calculate viewport-psnr given 'tileVer'
      //   int tmp_vp[2];
      //   // first frame & last frame
      //   tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
      //   tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
      //   if(tmp_vp[0] >= 180)
      //     tmp_vp[0] -= 360;
      //   if(tmp_vp[0] < -180)
      //     tmp_vp[0] += 360;
      //   //
      //   if(tmp_vp[1] >= 90)
      //     tmp_vp[1] -= 180;
      //   if(tmp_vp[1] <= -90)//
      //     tmp_vp[1] += 180;
      //   vp_psnr = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      //   tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
      //   tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
      //   if(tmp_vp[0] >= 180)
      //     tmp_vp[0] -= 360;
      //   if(tmp_vp[0] < -180)
      //     tmp_vp[0] += 360;
      //   //
      //   if(tmp_vp[1] >= 90)
      //     tmp_vp[1] -= 180;
      //   if(tmp_vp[1] <= -90)//
      //     tmp_vp[1] += 180;
      //   vp_psnr_last = est_vp_psnr(TILE_MSE, NO_TILE, tileVer, tmp_vp);
      //   // update selection
      //   if((vp_psnr + vp_psnr_last)/2 > max_vp_psnr){
      //     std::copy (tileVer, tileVer + tiles, selectTileVer);
      //     max_vp_psnr = (vp_psnr + vp_psnr_last)/2;
      //   }
      //   // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
      //   break; // not enough bandwidth
      // }
      // // assign remaining bandwidth for ext_1 tiles
      // FLAG = false;
      // while(remainBW > 0 && !FLAG){
      //   FLAG = true;
      //   for(i=0; i < NO_TILE; i++){
      //     if(vmask[i] == 2 && remainBW > 0){
      //       if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
      //         remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
      //         tileVer[i] ++;
      //         FLAG = false;
      //       }
      //     }
      //   }
      // }
      // // assign remaining bandwidth to non-visible tiles
      // if(remainBW > 0){
      //   FLAG = false;
      //   while(remainBW > 0 && !FLAG){
      //     FLAG = true;
      //     for(i=0; i < NO_TILE; i++){
      //       if(vmask[i] == 0 && remainBW > 0){
      //         if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
      //           remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
      //           tileVer[i] ++;
      //           FLAG = false;
      //         }
      //       }
      //     }
      //   }
      // }
      // calculate viewport-psnr given 'tileVer'
      // first frame & last frame
      tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
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
      tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
      // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
      //
      vp_ver ++;
    }//endwhile
  }
  else{ // CASE 2: EXT_WITH=2
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
      // determine quality if 'ext1' and 'ext2' are assigned same version // dung o day, ve thoaj
      // equally assign remaining bandwidth for 'ext_1' and 'ext_2'
     for(j=1; j < NO_VER; j++){
        tmp_sum = 0;
        for(i=0; i < No_tile; i++){
          if(vmask[i] == 2 || vmask[i] == 3)
            tmp_sum += (TILE_BR[i][j] - TILE_BR[i][tileVer[i]]);  
        }
        if(tmp_sum > remainBW) break;
      }
      mean_ver_2 = j - 1;
      for(i=0; i < No_tile; i++)
        if(vmask[i] == 2 || vmask[i] == 3)
          tileVer[i] = mean_ver_2;
      // FLAG = false;
      // while(remainBW > 0 && !FLAG){
      //   FLAG = true;
      //   for(i=0; i < NO_TILE; i++){
      //     if((vmask[i] == 2 || vmask[i] == 3) && remainBW > 0){
      //       if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
      //         remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
      //         tileVer[i] ++;
      //         FLAG = false;
      //       }
      //     }
      //   }
      // }
      // find version
      // for(i=0; i < tiles; i++)
      //   if((vmask[i] == 2 || vmask[i] == 3) && tileVer[i] < mean_ver_2)
      //     mean_ver_2 = tileVer[i];
      // assign remaining bandwidth to non-visible tiles
      // if(remainBW > 0){
      //   FLAG = false;
      //   while(remainBW > 0 && !FLAG){
      //     FLAG = true;
      //     for(i=0; i < NO_TILE; i++){
      //       if(vmask[i] == 0 && remainBW > 0){
      //         if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
      //           remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
      //           tileVer[i] ++;
      //           FLAG = false;
      //         }
      //       }
      //     }
      //   }
      // }
      //
      // calculate viewport-psnr given 'tileVer'
      // first frame & last frame
      tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
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
      tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
      tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
      // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
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
        // assign highest possible version for 'ext_2' tiles
        for(j=1; j < NO_VER; j++){
          tmp_sum = 0;
          for(i=0; i < No_tile; i++){
            if(vmask[i] == 3)
              tmp_sum += (TILE_BR[i][j] - TILE_BR[i][tileVer[i]]);  
          }
          if(tmp_sum > remainBW) break;
        }
      for(i=0; i < No_tile; i++)
          if(vmask[i] == 3)
            tileVer[i] = j-1;
        // FLAG = false;
        // while(remainBW > 0 && !FLAG){
        //   FLAG = true;
        //   for(i=0; i < NO_TILE; i++){
        //     if(vmask[i] == 3 && remainBW > 0){
        //       if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
        //         remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
        //         tileVer[i] ++;
        //         FLAG = false;
        //       }
        //     }
        //   }
        // }
        // // assign remaining bandwidth to non-visible tiles
        // if(remainBW > 0){
        //   FLAG = false;
        //   while(remainBW > 0 && !FLAG){
        //     FLAG = true;
        //     for(i=0; i < NO_TILE; i++){
        //       if(vmask[i] == 0 && remainBW > 0){
        //         if(tileVer[i] < NO_VER - 1 && (TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]]) < remainBW){
        //           remainBW -= TILE_BR[i][tileVer[i]+1] - TILE_BR[i][tileVer[i]];
        //           tileVer[i] ++;
        //           FLAG = false;
        //         }
        //       }
        //     }
        //   }
        // }
        // calculate viewport-psnr given 'tileVer'
        // first frame & last frame
        tmp_vp[0] = est_frame_vp[index * INTERVAL][0] + est_err[index][0];
        tmp_vp[1] = est_frame_vp[index * INTERVAL][1] + est_err[index][1];
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
        tmp_vp[0] = est_frame_vp[(index + 1) * INTERVAL - 1][0] + est_err[index][0];
        tmp_vp[1] = est_frame_vp[(index + 1) * INTERVAL - 1][1] + est_err[index][1];
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
        // showTileVersion(selectTileVer, No_tile_h, No_tile_v);
        //
        ext_1_ver ++;
      }
      vp_ver++;
    }
  }
  return selectTileVer;
}
int* DecisionEngine::Ireland(int index){
  int i,j, ext_w, k;
  int select_ver;
  int* tile_ver = new int[No_tile];
  double* tile_w = new double[No_tile];
  double* tile_br = new double[No_tile];
  double* tile_d = new double[No_tile];
  double* k_i = new double[No_tile];
  int* vmask, *vmask_cur, *pixel;
  int* vmask_ext;
  double BR_budget = est_seg_thrp[index];
  double BR_budget_visi;
  double BR_budget_invi;
  double BR_used = 0;
  int NT = No_tile;
  int NT_W = No_tile_v;
  int NT_H = No_tile_h;
  int W_T = tile_W;
  int W_H = tile_H;
  int VP_T = 960;
  int VP_H = 960;
  int W = 3840;
  int H = 1920;
  double lamda = 0.8, u, v, m, n;
  int tid;
  double tmp, tmp2;
  decide_width[index] = 0;
  if(index <= BUFF/INTER){
    for(i=0; i < NT; i++){
      tile_ver[i] = 0;
    }
    return tile_ver;
  }
  /* marking */
  vmask = get_visible_tile(est_vp[index]);
  pixel = get_visible_pixel(est_vp[index]);
  // printf("#[Ireland][index=%d] visible budget=%.2f (%d, %d)\n", index, BR_budget, est_vp[index][0], est_vp[index][1]);
  // printf("#[Ireland][index=%d] vmask:\n", index);
  //showArrayInt(vmask, NT);
  // printf("#[Ireland][index=%d] pixel:\n", index);
  //showArrayInt(pixel, NT);
  // check if the available bw is enough for the lowest version
  tmp = 0;
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    tmp += TILE_SEG_BR[index][i][tile_ver[i]];
  }
  if(tmp >= BR_budget)
    return tile_ver;
  /* assign bitrate for each visible/invisible areas */
  tmp = 0;
  tmp2 = 0;
  for(i=0; i < NT; i++){
    if(vmask[i] == 0)
      tmp += TILE_SEG_BR[index][i][0];
    if(vmask[i] == 1)
      tmp2 += TILE_SEG_BR[index][i][0];
  }
  if(lamda * BR_budget >= tmp2 && (1-lamda)*BR_budget >= tmp){
    BR_budget_invi = (1-lamda) * BR_budget;
    BR_budget_visi = lamda * BR_budget;
  }else{
    BR_budget_invi = tmp;
    BR_budget_visi = BR_budget - tmp;
  }

  // if(tmp > (1-lamda) * BR_budget){
  //   BR_budget_invi = tmp;
  //   BR_budget_visi = BR_budget - BR_budget_invi;
  // }else{
  //   BR_budget_invi = (1-lamda) * BR_budget;
  //   BR_budget_visi = lamda * BR_budget;
  // }
  /* assign bitrate for visible tiles */
  // assign the lowest version to each tile inadvance
  for(i=0; i < NT; i++){
    if(vmask[i] == 1){
      /* calculate tile weight and allocated bitrate*/
      tile_w[i] = pixel[i] * 1.0 / (VP_T * VP_H);
      tile_br[i] = BR_budget_visi * tile_w[i];
      /* decide the version of each tile */
      j = 0;
      while(j < NO_VER && TILE_SEG_BR[index][i][j] < tile_br[i]) j++;
      if(j==0){ /*allocated bitrate is not enough even for the lowest version */
        // printf("#[Ireland] tile #%d\n", i);
        tile_ver[i] = 0;//
      }else{
        tile_ver[i] = j-1;
      }
      tmp+=tile_w[i];
      BR_used += TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  printf("sum_weight=%.2f BR_used=%.2f\n", tmp, BR_used);
  showTileVersion(tile_ver, No_tile_h, No_tile_v);
  // printf("#[Ireland][index=%d] Visible tiles version:\n", index);
  // showArrayInt(tile_ver, NT);
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
  // u = est_vp[index][0]/(2*M_PI) + 0.5;
  // v = 0.5 - est_vp[index][1]/M_PI;
  // m = u * W  - 0.5;
  // n = v * H - 0.5;
  // for(k=0; k < No_face; k++){
  //  for(i=0; i < NT_H; i++){
  //    for(j=0; j < NT_W; j++){
  //      if(vmask[k* NT_H * NT_W + i * NT_W + j] == 0){
  //  /* calculate distace to viewport center*/
  //  tile_d[k* NT_H * NT_W + i * NT_W + j] = calc_distance(est_vp[index], k, i, j, No_face, No_tile_h, No_tile_v, face_W, face_H);
  //  /* find the max. distance */
  //  if(tile_d[i*NT_W + j] > max_d)
  //     max_d = tile_d[i*NT_W + j];//
  //      }
  //    }
  //  }
  // }
  for(tid = 0; tid < NT; tid++){
    if(vmask[tid] == 0){
      tile_d[tid] = calc_distance(est_vp[index], tid, No_face, No_tile_h, No_tile_v, face_W, face_H);
      /* find the longest distance */
      if(tile_d[tid] > max_d)
        max_d = tile_d[tid];
    }else{
      tile_d[tid] = 0; // distance is zero if it is a visible tile
    }
  }
  // printf("#[Ireland][index=%d] Invisible tiles distances:\n", index);
  // showArrayDouble(tile_d, NT);
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
      tile_br[i] = BR_budget_invi * k_i[i]/sum_d;
      j = 0;
      while(j < NO_VER && TILE_SEG_BR[index][i][j] < tile_br[i]) j++;
      if(j==0){ /*allocated bitrate is not enough even for the lowest version */
        tile_ver[i] = 0;//
      }else{
        tile_ver[i] = j-1;
      }
      BR_used += TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  printf("BR_budget=%.2f\tBR_used=%.2f\n", BR_budget, BR_used);
  /*
  for(i=0; i < NT; i++){
    if(vmask[i] == 0)
      printf("%.2f ", k_i[i]/sum_d);
    else
      printf("%.2f- ", tile_w[i]);
    if((i+1) % No_tile_h == 0)
      printf("\n");
  }
  */
  // printf("#[Ireland][index=%d] outside budget=%.2f used=%.2f\n", index,(1-lamda) * BR_budget, BR_used);
  // printf("#[Ireland][index=%d] Invisible tiles version:\n", index);
  // showArrayDouble(tile_br, NT);
  // showArrayInt(tile_ver, NT);
  // for(i=0; i < NT; i++)
  // printf("%.2f\t", TILE_SEG_BR[index][i][tile_ver[i]]);
  // printf("\n");
  // printf("#[Ireland][index=%d] Invisible tiles distances:\n", index);
  // for(i=0; i < NT_H; i++){
  //   for(j=0; j < NT_W; j++)
  //     printf("%d ", tile_ver[i*NT_W + j]);
  //   printf("\t");
  //   for(j=0; j < NT_W; j++)
  //     printf("%.2f ", TILE_SEG_BR[index][i*NT_W + j][tile_ver[i*NT_W + j]]);
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
// modify to remove the problem in which the total tile BR exceeds available BW
int* DecisionEngine::Ireland_v2(int index){
  int i,j, ext_w, k;
  int select_ver;
  int* tile_ver = new int[No_tile];
  double* tile_w = new double[No_tile];
  double* tile_br = new double[No_tile];
  double* tile_d = new double[No_tile];
  double* k_i = new double[No_tile];
  int* vmask, *vmask_cur, *pixel;
  int* vmask_ext;
  double BR_budget = est_seg_thrp[index];
  double BR_budget_visi;
  double BR_budget_invi;
  double BR_used = 0;
  int NT = No_tile;
  int NT_W = No_tile_v;
  int NT_H = No_tile_h;
  int W_T = tile_W;
  int W_H = tile_H;
  int VP_T = 960;
  int VP_H = 960;
  int W = 3840;
  int H = 1920;
  double lamda = 0.8, u, v, m, n;
  int tid;
  double tmp, tmp2;
  decide_width[index] = 0;
  if(index <= BUFF/INTER){
    for(i=0; i < NT; i++){
      tile_ver[i] = 0;
    }
    return tile_ver;
  }
  /* marking */
  vmask = get_visible_tile(est_vp[index]);
  pixel = get_visible_pixel(est_vp[index]);
  // printf("#[Ireland][index=%d] visible budget=%.2f (%d, %d)\n", index, BR_budget, est_vp[index][0], est_vp[index][1]);
  // printf("#[Ireland][index=%d] vmask:\n", index);
  //showArrayInt(vmask, NT);
  // printf("#[Ireland][index=%d] pixel:\n", index);
  //showArrayInt(pixel, NT);
  // check if the available bw is enough for the lowest version
  tmp = 0;
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    tmp += TILE_SEG_BR[index][i][tile_ver[i]];
  }
  if(tmp >= BR_budget)
    return tile_ver;
  /* assign bitrate for each visible/invisible areas */
  tmp = 0;
  tmp2 = 0;
  for(i=0; i < NT; i++){
    if(vmask[i] == 0)
      tmp += TILE_SEG_BR[index][i][0];
    if(vmask[i] == 1)
      tmp2 += TILE_SEG_BR[index][i][0];
  }
  if(lamda * BR_budget >= tmp2 && (1-lamda)*BR_budget >= tmp){
    BR_budget_invi = (1-lamda) * BR_budget;
    BR_budget_visi = lamda * BR_budget;
  }else{
    BR_budget_invi = tmp;
    BR_budget_visi = BR_budget - tmp;
  }
  // if(tmp > (1-lamda) * BR_budget){
  //   BR_budget_invi = tmp;
  //   BR_budget_visi = BR_budget - BR_budget_invi;
  // }else{
  //   BR_budget_invi = (1-lamda) * BR_budget;
  //   BR_budget_visi = lamda * BR_budget;
  // }
  /* assign bitrate for visible tiles */
  // assign the lowest version to each tile inadvance
  for(i=0; i < NT; i++){
    if(vmask[i] == 1){
      tile_ver[i] = 0;
      BR_budget_visi -= TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  for(i=0; i < NT; i++){
    if(vmask[i] == 1){
      /* calculate tile weight and allocated bitrate*/
      tile_w[i] = pixel[i] * 1.0 / (VP_T * VP_H);
      tile_br[i] = BR_budget_visi * tile_w[i];
      /* decide the version of each tile */
      j = 1;
      while(j < NO_VER && (TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][tile_ver[i]]) < tile_br[i]) j++;
      if(j > 1)
        tile_ver[i] = j-1;
      tmp+=tile_w[i];
      BR_used += TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  printf("sum_weight=%.2f BR_used=%.2f\n", tmp, BR_used);
  // showTileVersion(tile_ver, No_tile_h, No_tile_v);
  // printf("#[Ireland][index=%d] Visible tiles version:\n", index);
  // showArrayInt(tile_ver, NT);
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
  // u = est_vp[index][0]/(2*M_PI) + 0.5;
  // v = 0.5 - est_vp[index][1]/M_PI;
  // m = u * W  - 0.5;
  // n = v * H - 0.5;
  // for(k=0; k < No_face; k++){
  //  for(i=0; i < NT_H; i++){
  //    for(j=0; j < NT_W; j++){
  //      if(vmask[k* NT_H * NT_W + i * NT_W + j] == 0){
  //  /* calculate distace to viewport center*/
  //  tile_d[k* NT_H * NT_W + i * NT_W + j] = calc_distance(est_vp[index], k, i, j, No_face, No_tile_h, No_tile_v, face_W, face_H);
  //  /* find the max. distance */
  //  if(tile_d[i*NT_W + j] > max_d)
  //     max_d = tile_d[i*NT_W + j];//
  //      }
  //    }
  //  }
  // }
  for(tid = 0; tid < NT; tid++){
    if(vmask[tid] == 0){
      tile_d[tid] = calc_distance(est_vp[index], tid, No_face, No_tile_h, No_tile_v, face_W, face_H);
      /* find the longest distance */
      if(tile_d[tid] > max_d)
        max_d = tile_d[tid];
    }else{
      tile_d[tid] = 0; // distance is zero if it is a visible tile
    }
  }
  // printf("#[Ireland][index=%d] Invisible tiles distances:\n", index);
  // showArrayDouble(tile_d, NT);
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
      tile_ver[i] = 0;
      BR_budget_invi -= TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  for(i=0; i < NT; i++){
    if(vmask[i] == 0){
      tile_br[i] = BR_budget_invi * k_i[i]/sum_d;
      j = 1;
      while(j < NO_VER && (TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][tile_ver[i]]) < tile_br[i]) j++;
      if(j > 1)
        tile_ver[i] = j-1;
      BR_used += TILE_SEG_BR[index][i][tile_ver[i]];
    }
  }
  printf("BR_budget=%.2f\tBR_used=%.2f\n", BR_budget, BR_used);
  return tile_ver;
}
int* DecisionEngine::Ghent(int index){
  int i,j, ext_w;
  int select_ver;
  int* vmask, *vmask_cur;
  int* vmask_ext;
  double BR_budget = est_seg_thrp[index];
  int NT = No_tile;
  int NT_W = No_tile_v;
  int NT_H = No_tile_h;
  int* tile_ver = new int[NT];
  printf("#[Ghent] index=%d\n", index);
  decide_width[index] = 1;
  // marking tiles
  ext_w = 1;
  vmask = get_visible_tile(est_vp[index]);
  /* add visible tiles of the current viewport to 'vmask'*/
  vmask_cur = get_visible_tile(cur_vp[index]);
  for(i=0; i < NT; i++)
    if(vmask_cur[i] == 1)
      vmask[i] = 1;
  /* extend vmask by 1 tile in all directions --> adjcent tiles*/
  vmask_ext = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_w);
  //check point #1:
  // printf("#[Ghent][index=%d] BR_budget=%.2f, est_vp=(%d,%d), est_vp=(%d, %d)\n", index, BR_budget, est_vp[index][0], est_vp[index][1], cur_vp[index][0], cur_vp[index][1]);
  // for(i=0; i < NT_H; i++){
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", vmask[i*NT_W + j]);
  //  printf("\t");
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", vmask_cur[i*NT_W + j]);
  //  printf("\t");
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", vmask_ext[i*NT_W + j]);
  //  printf("\n");
  // }
  // assigning tile's bitrates
  /* assign the lowest version to all tiles*/
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    BR_budget -= TILE_SEG_BR[index][i][tile_ver[i]]; 
  }
  if(BR_budget < 0) /* return if avail. bw is not enough*/
    return tile_ver;
  /* assign highest possible version for viewport tiles */
  double sum_br;
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  //check point #2:
  // printf("#[Ghent][index=%d] viewport tiles, BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", tile_ver[i*NT_W + j]);
  //  printf("\n");
  // }
  /* assign highest possible version for adjcent tiles */
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  // printf("#[Ghent][index=%d] adjecent BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", tile_ver[i*NT_W + j]);
  //  printf("\n");
  // }
  /* assign highest possible version for outside tiles */
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  // printf("#[Ghent][index=%d] outside BR_budget=%.2f\n", index, BR_budget);
  // for(i=0; i < NT_H; i++){
  //  for(j=0; j < NT_W; j++)
  //    printf("%d ", tile_ver[i*NT_W + j]);
  //  printf("\n");
  // }
  return tile_ver;
}
int* DecisionEngine::Ghent_opt(int index){
  int i,j, ext_w;
  int select_ver;
  int* vmask, *vmask_cur;
  int* vmask_ext;
  double BR_budget = est_seg_thrp[index];
  int NT = No_tile;
  int NT_W = No_tile_v;
  int NT_H = No_tile_h;
  int* tile_ver = new int[NT];
  printf("#[Ghent] index=%d\n", index);
  decide_width[index] = 1;
  bool FLAG;
  // marking tiles
  ext_w = 1;
  vmask = get_visible_tile(est_vp[index]);
  /* add visible tiles of the current viewport to 'vmask'*/
  vmask_cur = get_visible_tile(cur_vp[index]);
  for(i=0; i < NT; i++)
    if(vmask_cur[i] == 1)
      vmask[i] = 1;
  /* extend vmask by 1 tile in all directions --> adjcent tiles*/
  vmask_ext = extVmask(vmask, No_face, No_tile_h, No_tile_v, ext_w);
  /* assign the lowest version to all tiles*/
  for(i=0; i < NT; i++){
    tile_ver[i] = 0;
    BR_budget -= TILE_SEG_BR[index][i][tile_ver[i]]; 
  }
  if(BR_budget < 0) /* return if avail. bw is not enough*/
    return tile_ver;
  /* assign highest possible version for viewport tiles */
  double sum_br;
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  /* assign highest possible version for adjcent tiles */
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  /* assign highest possible version for outside tiles */
  for(j=1; j < NO_VER; j++){
    sum_br = 0;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0)
        sum_br += TILE_SEG_BR[index][i][j] - TILE_SEG_BR[index][i][0];
    }
    if(sum_br > BR_budget) break;
  }
  if(j > 1 &&  BR_budget > 0){
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0){
        BR_budget -= (TILE_SEG_BR[index][i][j-1] - TILE_SEG_BR[index][i][tile_ver[i]]); 
        tile_ver[i] = j-1;
      }
    }
  }
  // if there is remaining BW, used them to improve tiles' versions
  FLAG = true;
  while(BR_budget > 0 && FLAG){
    FLAG = false;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 1 && tile_ver[i] < NO_VER-1 && (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]) < BR_budget){
        BR_budget -= (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]);
        tile_ver[i]++;
        FLAG = true;
        if(BR_budget <= 0) break;
      }
    }
  }
  FLAG = true;
  while(BR_budget > 0 && FLAG){
    FLAG = false;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 2 && tile_ver[i] < NO_VER-1 && (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]) < BR_budget){
        BR_budget -= (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]);
        tile_ver[i]++;
        FLAG = true;
        if(BR_budget <= 0) break;
      }
    }
  }
  FLAG = true;
  while(BR_budget > 0 && FLAG){
    FLAG = false;
    for(i=0; i < NT; i++){
      if(vmask_ext[i] == 0 && tile_ver[i] < NO_VER-1 && (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]) < BR_budget){
        BR_budget -= (TILE_SEG_BR[index][i][tile_ver[i]+1] - TILE_SEG_BR[index][i][tile_ver[i]]);
        tile_ver[i]++;
        FLAG = true;
        if(BR_budget <= 0) break;
      }
    }
  }
  return tile_ver;
}
int* DecisionEngine::BellLab(int index){
  int i,j,q;
  int rows = No_tile_v;
  int cols = No_tile_h;
  int tiles = rows * cols;
  double **vprob;
  double speed;
  //
  double used_bandwidth = 0;
  int *tileVer = new int[tiles];
  for(i=0; i < tiles; i++){ // allocate a minimum quality for all tiles
    tileVer[i] = 0;
    used_bandwidth += TILE_SEG_BR[index][i][tileVer[i]];
  }
  // first interval
  if(index <= BUFF/INTER)
    return tileVer;
  speed = ang_speed[index];
  // speed = 2; // degree/frame
  // generate gaussian filter kernel
  double sigma = (speed / 0.033) / 180.0 * 3.14  * BUFF * 1.0 / INTER;
  if(sigma > 0){
    int filter_size = 2 * ceil(sigma * 2) + 1;
    double** kernel = filter(filter_size, sigma);
    // printf("############ index:%d\n", index);
    // printf("############ KERNEL:\n");
    // for(i=0; i < filter_size; i++){
    //   for(j=0; j < filter_size; j++)
    //     printf("%.4f ", kernel[i][j]);
    //   printf("\n");
    // }
    // compute tile visible probabilities
    int *vmask = get_visible_tile(est_vp[index]);
    double** vmask_2D = new double*[rows];
    for(i=0; i < rows; i++){
      vmask_2D[i] = new double[cols];
      for(j=0; j < cols; j++)
        vmask_2D[i][j] = vmask[i * rows + j];
    }
    // printf("############ VMASK:\n");
    // for(i=0; i < rows; i++){
    //   for(j=0; j < cols; j++)
    //     printf("%.2f ", vmask_2D[i][j]);
    //   printf("\n");
    // }
    vprob = applyFilter(vmask_2D, kernel, rows, filter_size); // visible probabilities of tiles
    // printf("############ vprob:\n");
    // for(i=0; i < rows; i++){
    //   for(j=0; j < cols; j++)
    //     printf("%.4f ", vprob[i][j]);
    //   printf("\n"); 
    // }

  }else{
    int *vmask = get_visible_tile(est_vp[index]);
    double** vmask_2D = new double*[rows];
    for(i=0; i < rows; i++){
      vmask_2D[i] = new double[cols];
      for(j=0; j < cols; j++)
        vmask_2D[i][j] = vmask[i * rows + j];
    }
    vprob = vmask_2D;
  }
  //  //
  double** A = init2dArrayDouble(4, No_tile * (NO_VER - 1));
  for(i=0; i < No_tile; i++){
    // printf("#[Tile #%d]: ", i);
    for(q=1; q < NO_VER; q++){
      A[0][i * (NO_VER - 1) + q - 1] = i;
      A[1][i * (NO_VER - 1) + q - 1] = q;
      A[2][i * (NO_VER - 1) + q - 1] = Uti[index][i][q] * vprob[i/rows][i%rows]/Cost[index][i][q];
      A[3][i * (NO_VER - 1) + q - 1] = Cost[index][i][q];
      // printf("%.2f\t%.2f\t", A[2][i * (NO_VER - 1) + q - 1], A[3][i * (NO_VER - 1) + q - 1]);
    }
    // printf("\n");
  }
  // sort A[2] in ascending order
  std::vector<double> array(A[2], A[2] + No_tile * (NO_VER - 1));
  std::vector<int> indices = sort_index(array);
  //
  // printf("sorted array:\n");//
  //    for (std::vector<int>::iterator it = indices.begin(); it != indices.end(); ++it){
  //       std::cout << *it << " " << array[*it] << "\n";
  //    }
  //
  //calculate the selected version for each tile
  for (std::vector<int>::iterator it = indices.begin(); it != indices.end(); ++it){
    if(used_bandwidth + A[3][*it] <= est_seg_thrp[index] && tileVer[(int)A[0][*it]] == A[1][*it] - 1){
      tileVer[(int)A[0][*it]] = A[1][*it];
      used_bandwidth += A[3][*it];
    }
  }
  // printf("############ TILE VERSION:\n");
  // for(i=0; i < rows; i++){
  //   for(j=0; j < cols; j++)
  //     printf("%d ", tileVer[i * rows + j]);
  //   printf("\n");
  // }
  // printf("############ TILE VERSION:\n");
  return tileVer;
}
/*
Use gurobi to obtain more exact result
*/
// int* DecisionEngine::BellLab_v2(int index){
//   int i,j,q;
//   int rows = No_tile_v;
//   int cols = No_tile_h;
//   int tiles = rows * cols;
//   double **vprob;
//   double speed;
//   //
//   double used_bandwidth = 0;
//   int *tileVer = new int[tiles];
//   for(i=0; i < tiles; i++){ // allocate a minimum quality for all tiles
//     tileVer[i] = 0;
//     used_bandwidth += TILE_SEG_BR[index][i][tileVer[i]];
//   }
//   // first interval
//   if(index <= BUFF/INTER)
//     return tileVer;
//   speed = ang_speed[index];
//   // speed = 2; // degree/frame
//   // generate gaussian filter kernel
//   double sigma = (speed / 0.033) / 180.0 * 3.14  * BUFF * 1.0 / INTER;
//   if(sigma > 0){
//     int filter_size = 2 * ceil(sigma * 2) + 1;
//     double** kernel = filter(filter_size, sigma);
//     // printf("############ index:%d\n", index);
//     // printf("############ KERNEL:\n");
//     // for(i=0; i < filter_size; i++){
//     //   for(j=0; j < filter_size; j++)
//     //     printf("%.4f ", kernel[i][j]);
//     //   printf("\n");
//     // }
//     // compute tile visible probabilities
//     int *vmask = get_visible_tile(est_vp[index]);
//     double** vmask_2D = new double*[rows];
//     for(i=0; i < rows; i++){
//       vmask_2D[i] = new double[cols];
//       for(j=0; j < cols; j++)
//         vmask_2D[i][j] = vmask[i * rows + j];
//     }
//     // printf("############ VMASK:\n");
//     // for(i=0; i < rows; i++){
//     //   for(j=0; j < cols; j++)
//     //     printf("%.2f ", vmask_2D[i][j]);
//     //   printf("\n");
//     // }
//     vprob = applyFilter(vmask_2D, kernel, rows, filter_size); // visible probabilities of tiles
//     // printf("############ vprob:\n");
//     // for(i=0; i < rows; i++){
//     //   for(j=0; j < cols; j++)
//     //     printf("%.4f ", vprob[i][j]);
//     //   printf("\n"); 
//     // }

//   }else{
//     int *vmask = get_visible_tile(est_vp[index]);
//     double** vmask_2D = new double*[rows];
//     for(i=0; i < rows; i++){
//       vmask_2D[i] = new double[cols];
//       for(j=0; j < cols; j++)
//         vmask_2D[i][j] = vmask[i * rows + j];
//     }
//     vprob = vmask_2D;
//   }
//   // build objective function and solve the problem
//   GRBVar **x;
//   GRBVar **x_name;
//    try {
//     GRBEnv env = GRBEnv();
//     GRBModel model = GRBModel(env);
//     // init
//     x = new GRBVar*[No_tile];
//     x_name = new char**[No_tile];
//     for(i=0; i < No_tile; i++){
//         x[i] = new GRBVar[NO_VER];
//         x_name[i] = new char*[NO_VER];
//         for(j=0; j < NO_VER; j++)
//             x_name[i][j] = new char[LEN];
//     }
//     // Create variables
//     for(i=0; i < No_tile; i++){
//         for(j=0; j < NO_VER; j++){
//             x_name[i][j] = new char[LEN];
//             sprintf(x_name[i][j], "x_%d_%d", i, j);
//             x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x_name[i][j]);
//         }
//     }
//     // objective funtion
//     GRBQuadExpr obj;
//     GRBQuadExpr avgVPQL = 0;
//     GRBQuadExpr varVPQL = 0;
//     GRBQuadExpr muy_2 = 0;
//     GRBQuadExpr tmp = 0;
//     for(i=0; i < No_tile; i++){
//         act_w[i] = p[i] * w[i];
//         for(j=0; j < NO_VER; j++){
//             avgVPQL += act_w[i]*d[i][j]*x[i][j];
//         }
//     }
//     obj = avgVPQL;
//     // obj = avgVPQL + obj_w * varVPQL;
//     // model.setObjective(avgVPQL, GRB_MINIMIZE);
//     model.setObjective(obj, GRB_MINIMIZE);
//     // set constraints
//     // c1: total tiles' bitrates does not exceed bandwidth
//     GRBLinExpr constr_1 = 0;
//     GRBLinExpr* constr_2 = new GRBLinExpr[No_tile];
//     for(i=0; i < No_tile; i++){
//         constr_2[i] = 0;
//         for(j=0; j < NO_VER; j++){
//             constr_1 += TILE_SEG_BR[index][i][j] * x[i][j];
//             constr_2[i] += x[i][j];
//         }
//         model.addConstr(constr_2[i] == 1, "constr_2");
//     }
//     model.addConstr(constr_1 <= est_seg_thrp[index], "constr_1");
//     // Optimize model
//     model.optimize();
//     for(i=0; i < No_tile; i++){
//         for(j=0; j < NO_VER; j++){
//             if(x[i][j].get(GRB_DoubleAttr_X) == 1){
//                // cout << j << " ";
//               tileVer[i] = j;
//             }
//         }
//         // if((i+1)%No_tile_h == 0)
//         //     cout << endl;
//     }

//     // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

//   } catch(GRBException e) {
//     cout << "Error code = " << e.getErrorCode() << endl;
//     cout << e.getMessage() << endl;
//   } catch(...) {
//     cout << "Exception during optimization" << endl;
//   }

//   return tileVer;
// }

int* DecisionEngine::ProbDASH(int index){
  int* tileVer = new int[No_tile];
  int i,j,ii,jj;
  double obj_w = 0.0015;
  // double mean_alpha = -0.54, var_alpha = 7.03, mean_beta=0.18, var_beta=2.55, mean_gama=2.16, var_gama=0.15;
  double mean_alpha = -0.02, var_alpha = 51.12, mean_beta=1.05, var_beta= 21.79, mean_gama=2.16, var_gama=0.15;
  double m_a, m_b, m_g;
  m_a = speed[index][0];
  m_b = speed[index][1];
  m_g = 0;
  // double* p = calc_tile_view_prob(mean_alpha, var_alpha, mean_beta, var_beta, mean_gama, var_gama); // viewing probability of a tile
  double* p = calc_tile_view_prob_2(index, mean_alpha, var_alpha, mean_beta, var_beta, mean_gama, var_gama);
  double* w = calc_tile_weight(tile_W, tile_H, W, H); // weight of a tile (ratio between erp_area and actual area)
  double* act_w  = new double[No_tile];
  int LEN = 200;
  GRBVar** x;
  char*** x_name;
  double** d = TILE_SEG_MSE[index]; // tiles' distortions
  // double** d = TILE_SEG_PSNR[index]; // tiles' distortions
  // 
  if(index < BUFF/INTER){
    for(i=0; i < No_tile; i++)
      tileVer[i] = 0;
    return tileVer;
  }
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    // init
    x = new GRBVar*[No_tile];
    x_name = new char**[No_tile];
    for(i=0; i < No_tile; i++){
        x[i] = new GRBVar[NO_VER];
        x_name[i] = new char*[NO_VER];
        for(j=0; j < NO_VER; j++)
            x_name[i][j] = new char[LEN];
    }
    // Create variables
    for(i=0; i < No_tile; i++){
        for(j=0; j < NO_VER; j++){
            x_name[i][j] = new char[LEN];
            sprintf(x_name[i][j], "x_%d_%d", i, j);
            x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x_name[i][j]);
        }
    }
    // objective funtion
    GRBQuadExpr obj;
    GRBQuadExpr avgVPQL = 0;
    GRBQuadExpr varVPQL = 0;
    GRBQuadExpr muy_2 = 0;
    GRBQuadExpr tmp = 0;
    for(i=0; i < No_tile; i++){
        act_w[i] = p[i] * w[i];
        for(j=0; j < NO_VER; j++){
            avgVPQL += act_w[i]*d[i][j]*x[i][j];
        }
    }
    // for(i=0; i < No_tile; i++){
    //     for(j=0; j < NO_VER; j++){
    //       // first term
    //       varVPQL += act_w[i]*act_w[i]*d[i][j]*x[i][j]*x[i][j];
    //       // second term
    //       for(ii=0; ii < No_tile; ii++){
    //         for(jj=0; jj < NO_VER; jj++){
    //           varVPQL -= 2 * d[i][j] * x[i][j] * act_w[ii] * d[ii][jj] * x[ii][jj];
    //         }
    //       }
    //     }
    // }
    // // muy * muy_2                    
    // for(i=0; i < No_tile; i++){
    //   // sum(d_ij*x_ji)^2
    //   for(j=0; j < NO_VER; j++){
    //     // muy_2 += act_w[i] * act_w[i] * d[i][j] * d[i][j] * x[i][j] * x[i][j];
    //     muy_2 +=  d[i][j] * d[i][j] * x[i][j] * x[i][j];
    //     for(jj=0; jj < NO_VER && jj != j; jj++){
    //       // muy_2 += 2 * act_w[i] * act_w[i] * d[i][j] * x[i][j] * d[i][jj] * x[i][jj]; 
    //       muy_2 += 2 * d[i][j] * x[i][j] * d[i][jj] * x[i][jj]; 
    //     }
    //   }
    //   muy_2 *= act_w[i] * act_w[i];
    // }

    // for(i=0; i < No_tile; i++){
    //   for(ii=0; ii < No_tile && ii != i; ii++){
    //     for(j=0; j < NO_VER; j++){
    //       for(jj=0; jj < NO_VER; jj++){
    //         tmp += act_w[i] * act_w[i] * d[i][j] * x[i][j] * d[ii][jj] * x[ii][jj];
    //       }
    //     }
    //   }
    // }
    // muy_2 += 2 * tmp;
    // varVPQL += No_tile * NO_VER * muy_2;
    obj = avgVPQL;
    // obj = avgVPQL + obj_w * varVPQL;
    // model.setObjective(avgVPQL, GRB_MINIMIZE);
    model.setObjective(obj, GRB_MINIMIZE);
    // set constraints
    // c1: total tiles' bitrates does not exceed bandwidth
    GRBLinExpr constr_1 = 0;
    GRBLinExpr* constr_2 = new GRBLinExpr[No_tile];
    for(i=0; i < No_tile; i++){
        constr_2[i] = 0;
        for(j=0; j < NO_VER; j++){
            constr_1 += TILE_SEG_BR[index][i][j] * x[i][j];
            constr_2[i] += x[i][j];
        }
        model.addConstr(constr_2[i] == 1, "constr_2");
    }
    model.addConstr(constr_1 <= all_seg_br[index], "constr_1");
    // Optimize model
    model.optimize();
    for(i=0; i < No_tile; i++){
        for(j=0; j < NO_VER; j++){
            if(x[i][j].get(GRB_DoubleAttr_X) == 1){
               // cout << j << " ";
              tileVer[i] = j;
            }
        }
        // if((i+1)%No_tile_h == 0)
        //     cout << endl;
    }

    // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return tileVer;
}
double* DecisionEngine::calc_tile_view_prob(double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama){
  int i,j, ii, m, n, tid, cnt;
  int h,w;
  struct timeval t_start, t_end;
  double phi, theta, tmp_sum;
  double R = W / (2 * M_PI); // radius of the sphere
  double* s = new double[No_tile];
  double* weight = new double[No_tile];
  double delta_theta = tile_H * M_PI / H;
  double delta_phi = tile_W * 2 * M_PI / W;
  double* tileProb = new double[No_tile];
  for(tid=0; tid < No_tile; tid++){
    tileProb[tid] = 0;
    w = tid % No_tile_h * tile_W;
    h = tid / No_tile_h * tile_H;
    gettimeofday(&t_start, NULL);
    for(m=w; m < w+tile_W; m++){
      for(n=h; n < h+tile_H; n++){
        phi = (m + 0.5) * 2 * M_PI / W - M_PI;
        theta = M_PI/2 - (n+0.5) * M_PI / H;
        tileProb[tid] += calc_Ps(phi, theta, mean_alpha, var_alpha, mean_beta, var_beta, mean_gama, var_gama)/(tile_W * tile_H);
        // malloc_stats(); 
      }
    }
    gettimeofday(&t_end, NULL);
    printf("#[calc_tile_view_prob] %.2f\n", tvdiff_us(&t_end, &t_start)/1000.0);
  }
  return tileProb;
}
double* DecisionEngine::calc_tile_view_prob_2(int index, double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama){
  int i,j, ii, m, n, tid, cnt;
  int h,w;
  struct timeval t_start, t_end;
  double phi, theta, tmp_sum, P_alpha, P_beta, P_gama, phi_est, theta_est;
  int phi_deg, theta_deg, phi_est_deg, theta_est_deg;
  int err_phi, err_theta, est_phi, est_theta;
  double R = W / (2 * M_PI); // radius of the sphere
  int vp[2];
  // 
  // int No_tile = 72;
  // int No_tile_h = 12;
  // int No_tile_v = 6;
  // int tile_W = W / No_tile_h;
  // int tile_H = H / No_tile_v;
  // 
  double* s = new double[No_tile];
  double* weight = new double[No_tile];
  double delta_theta = tile_H * M_PI / H;
  double delta_phi = tile_W * 2 * M_PI / W;
  double* tileProb = new double[No_tile];
  double** probPoint = init2dArrayDouble(360,181);
  double totalProb;
  // 
  est_phi = est_vp[index][0];
  est_theta = est_vp[index][1];
  // for(est_phi = 0; est_phi <= 180; est_phi += 30){
    vp[0] = est_phi;
    vp[1] = est_theta;
    // showTileVersion(get_visible_tile(vp), No_tile_h, No_tile_v);
  // 
  for(phi_deg=-179; phi_deg <= 180; phi_deg++){
    for(theta_deg=-90; theta_deg <= 90; theta_deg ++){
      // printf("#[calc_orient_set] (phi,theta)=(%d,%d)\n", phi_deg, theta_deg)
      // convert to radian
      err_phi = phi_deg - est_phi;
        if(err_phi < -180)
          err_phi += 360;
        else if(err_phi > 180)
          err_phi -= 360;
        // theta
      err_theta = theta_deg - est_theta;
        // if(err_theta < -90)
        //   err_theta += 180;
        // else if(err_theta > 90)
        //   err_theta -= 180;
      P_alpha = (1.0/(var_alpha * sqrt(2*M_PI))) * exp(-(err_phi - mean_alpha) * (err_phi - mean_alpha) / (2 * var_alpha * var_alpha));
      P_beta = (1.0/(var_beta * sqrt(2*M_PI))) * exp(-(err_theta - mean_beta) * (err_theta - mean_beta) / (2 * var_beta * var_beta));
      P_gama = 1;
      probPoint[phi_deg+179][theta_deg+90] = P_alpha * P_beta * P_gama;
    }
  }
  //
  totalProb = 0;
  for(tid=0; tid < No_tile; tid++){
    tileProb[tid] = 0;
    w = tid % No_tile_h * tile_W;
    h = tid / No_tile_h * tile_H;
    gettimeofday(&t_start, NULL);
    for(m=w; m < w+tile_W; m++){
      for(n=h; n < h+tile_H; n++){
        phi = (m + 0.5) * 2 * M_PI / W - M_PI;
        theta = M_PI/2 - (n+0.5) * M_PI / H;
        phi_deg = (int) (phi * 180 / M_PI);
        // if(phi_deg < 0)
        //   phi_deg += 360;
        theta_deg = (int) (theta * 180 / M_PI);
        tileProb[tid] += probPoint[phi_deg+179][theta_deg+90]/(tile_W * tile_H);
        // malloc_stats(); 
      }
    }
    totalProb += tileProb[tid];
    gettimeofday(&t_end, NULL);
    // printf("#[calc_tile_view_prob] %.2f\n", tvdiff_us(&t_end, &t_start)/1000.0);
  }
  for(tid=0; tid < No_tile; tid++)
    tileProb[tid] = tileProb[tid]/totalProb; // normalize to range [0, 1]
  // printf("# (%d, %d)\n", est_phi, est_theta);
  // showTileInfo(tileProb, No_tile_h, No_tile_v);
 // }
  return tileProb;
}
double DecisionEngine::calc_Ps(double phi, double theta, double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama){
  int i, num_orient_set;
  struct timeval t_start, t_end;
  double Ps;
  gettimeofday(&t_start, NULL);
  struct Orient* L = calc_orient_set(phi, theta, &num_orient_set); // phi, theta are of radian scales.
  gettimeofday(&t_end, NULL);
  // printf("#[calc_Ps] %.2f\n", tvdiff_us(&t_start, &t_end)/1000.0);
  struct Orient* L_est = new Orient[num_orient_set];
  Ps = 0;
  for(i=0; i < num_orient_set; i++){
    // calculate L_est
    L_est[i].alpha = L[i].alpha;
    L_est[i].beta = L[i].beta;
    L_est[i].gama = L[i].gama;
    // calculate Ps
    Ps += calc_PE(L[i],L_est[i], mean_alpha, var_alpha, mean_beta, var_beta, mean_gama, var_gama)/num_orient_set;
  }
  free(L);
  free(L_est);
  return Ps;
}
struct Orient* DecisionEngine::calc_orient_set(double phi_test, double theta_test, int* N){
  // granularity (degree: 360 * 181)
  // for each possible orientation
  // check if (phi, theta) can be seen
  // update orientatation sets
  int phi_deg, theta_deg;
  Orient* orientSet = new Orient[360 * 181];
  double phi, theta;
  int m, n, i, cnt, m_tmp;
  double Fh = M_PI/2, Fv = M_PI/2;
  double u,v,x,y,z;
  int Fh_deg = 90, Fv_deg = 90;
  // cach 1: full-search
  // double* x_ = new double[vp_W * vp_H];
  // double* y_ = new double[vp_W * vp_H];
  // double* z_ = new double[vp_W * vp_H];
  // double* X = new double[vp_W * vp_H];
  // double* Y = new double[vp_W * vp_H];
  // double* Z = new double[vp_W * vp_H];
  // double* u_ = new double[vp_W * vp_H];
  // double* v_ = new double[vp_W * vp_H];
  // double* phi_ = new double[vp_W * vp_H];
  // double* theta_ = new double[vp_W * vp_H];
  // int* m_ = new int[vp_W * vp_H];
  // int* n_ = new int[vp_W * vp_H];
  // double R[3][3];
  // for(m=0; m < vp_W; m++){
  //   for(n=0; n < vp_H; n++){
  //     u = (m+0.5) * 2 * tan(Fh/2)/vp_W;
  //     v = (n+0.5) * 2 * tan(Fv/2)/vp_H;
  //     //
  //     x = u - tan(Fh/2);
  //     y = -v + tan(Fv/2);
  //     z = 1;
  //     // 
  //     x_[n*vp_W + m] = x / sqrt(x*x + y*y + z*z);
  //     y_[n*vp_W + m] = y / sqrt(x*x + y*y + z*z);
  //     z_[n*vp_W + m] = z / sqrt(x*x + y*y + z*z);
  //   }
  // }
  // // (phi, theta) -> (m,n)
  // u = phi_test / (2 * M_PI) + 0.5;
  // v = 0.5 - theta_test / M_PI;
  // m = (int)(u * W - 0.5);
  // if(m >= W)
  //   m -= W;
  // n = (int)(v * H - 0.5);
  // cnt = 0;
  // for(phi_deg=0; phi_deg < 360; phi_deg++){
  //   for(theta_deg=-90; theta_deg <= 90; theta_deg ++){
  //     // printf("#[calc_orient_set] (phi,theta)=(%d,%d)\n", phi_deg, theta_deg);
  //     phi = phi_deg * M_PI / 180;
  //     theta = theta_deg * M_PI / 180;
  //     // rotation matrix
  //     R[0][0] = cos(phi + M_PI/2);
  //     R[0][1] = -sin(phi + M_PI/2) * sin(theta);
  //     R[0][2] = sin(phi + M_PI/2) * cos(theta);
  //     // 
  //     R[1][0] = 0;
  //     R[1][1] = cos(theta);
  //     R[1][2] = sin(theta);
  //     // 
  //     R[2][0] = -sin(phi+M_PI/2);
  //     R[2][1] = -cos(phi+M_PI/2) * sin(theta);
  //     R[2][2] = cos(phi+M_PI/2) * cos(theta);
  //     // rotation
  //     for(i=0; i < vp_W * vp_H; i++){
  //       X[i] = R[0][0] * x_[i] + R[0][1] * y_[i] + R[0][2] * z_[i];
  //       Y[i] = R[1][0] * x_[i] + R[1][1] * y_[i] + R[1][2] * z_[i];
  //       Z[i] = R[2][0] * x_[i] + R[2][1] * y_[i] + R[2][2] * z_[i];
  //       // 
  //       phi_[i] = atan(-Z[i]/X[i]);
  //       phi_[i] += (X[i] < 0)?M_PI:0;
  //       theta_[i] = asin(Y[i]/sqrt(X[i]*X[i] + Y[i] * Y[i] + Z[i] * Z[i]));
  //       // 
  //       u_[i] = phi_[i] / (2 * M_PI) + 0.5;
  //       v_[i] = 0.5 - theta_[i]/M_PI;
  //       // 
  //       m_[i] = (int)(u_[i] * W - 0.5);
  //       m_[i] -= (m_[i] >= W)?W:0;
  //       n_[i] = (int)(v_[i] * H - 0.5);
  //       if(m_[i] == m && n_[i] == n){
  //         orientSet[cnt].alpha = phi;
  //         orientSet[cnt++].beta = theta;
  //         break;
  //       }
  //     }
  //   }
  // }
  // *N = cnt;
  // cach 2: assuming a rectangular observeble area.
  phi_deg = (int)(phi_test / M_PI * 180);
  theta_deg = (int)(theta_test / M_PI * 180);
  if(phi_deg < 0)
    phi_deg += 360;
  cnt = 0;
  // printf("#[calc_orient_set] phi: (%d->%d), theta: (%d->%d)\n", phi_deg - Fh_deg/2, phi_deg + Fh_deg/2, (theta_deg - Fv_deg/2 >= -90)?(theta_deg - Fv_deg/2):-90, (theta_deg + Fv_deg/2 <= 90)?(theta_deg + Fv_deg/2):90);
  for(m=phi_deg - Fh_deg/2; m <= phi_deg + Fh_deg/2; m++){
    for(n=((theta_deg - Fv_deg/2 >= -90)?(theta_deg - Fv_deg/2):-90); n <= ((theta_deg + Fv_deg/2 <= 90)?(theta_deg + Fv_deg/2):90); n++){
      if(m < 0) m_tmp = m + 360;
      else{
        if(m >= 360)
          m_tmp = m % 360;
        else
          m_tmp = m;
      }
      orientSet[cnt].alpha = m_tmp / 180.0 * M_PI; // phi
      orientSet[cnt].beta = n / 180.0 * M_PI; // theta
      orientSet[cnt++].gama = 0;
    }
  }
  *N = cnt;
  return orientSet;
}
double DecisionEngine::calc_PE(struct Orient L,struct Orient L_est, double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama){
  double P_alpha, P_beta, P_gama;
  P_alpha = (1.0/(var_alpha * sqrt(2*M_PI))) * exp(-(L.alpha - L_est.alpha) * (L.alpha - L_est.alpha) / (2 * var_alpha * var_alpha));
  P_beta = (1.0/(var_beta * sqrt(2*M_PI))) * exp(-(L.beta - L_est.beta) * (L.beta - L_est.beta) / (2 * var_beta * var_beta));
  P_gama = (1.0/(var_gama * sqrt(2*M_PI))) * exp(-(L.gama - L_est.gama) * (L.gama - L_est.gama) / (2 * var_gama * var_gama));
  P_gama = 1;
  return P_alpha * P_beta * P_gama;
}

double* DecisionEngine::calc_tile_weight(int tile_W, int tile_H, int W, int H){
  int i,j, ii, m, n, tid, cnt;
  int h,w;
  double phi, theta, tmp_sum;
  double R = W / (2 * M_PI); // radius of the sphere
  double* s = new double[No_tile];
  double* weight = new double[No_tile];
  double delta_theta = tile_H * M_PI / H;
  double delta_phi = tile_W * 2 * M_PI / W;
  std::vector<double> phi_arr[No_tile];
  std::vector<double> theta_arr[No_tile];
  tmp_sum = 0;
  bool FLAG;
  double* area = new double[No_tile_v/2];
  // original implementation according to 360ProbDASH papers (NOT CORRECT)
  cnt=0;
  for(theta=0; theta < M_PI/2 - delta_theta; theta += delta_theta){
    area[cnt++] = R*R*delta_phi*(sin(theta+delta_theta) - sin(theta));
    // printf("# theta=PI/%d area=%.2f\n", (int) (theta/delta_theta), area);
  }
  for(i=0; i < No_tile; i++){
    if(i < No_tile/2){
      s[i] = area[No_tile_v/2 - i/No_tile_h - 1];
    }else{
      s[i] = s[(No_tile_v-1 - i/No_tile_h) * No_tile_h];
    }
    tmp_sum += s[i];
  //   phi = M_PI/2 - h * M_PI / H;
  //   theta = w * 2 * M_PI / W;
  //   s[i] = delta_theta * R * R * (cos(phi + delta_phi) - cos(phi));
  //   tmp_sum += s[i];
  //   printf("(%.5f,%.5f)\t", theta, phi);
  //   if((i+1)%No_tile_h == 0)
  //     printf("\n");
  // }
  // DUC's implementation
  // for(tid=0; tid < No_tile; tid++){
  //   cnt = 0;
  //   s[tid] = 0;
  //   w = tid % No_tile_h * tile_W;
  //   h = tid / No_tile_h * tile_H;
  //   for(m=w; m < w+tile_W; m++){
  //     for(n=h; n < h+tile_H; n++){
  //       phi = (m + 0.5) * 2 * M_PI / W - M_PI;
  //       theta = M_PI/2 - (n+0.5) * M_PI / H;
  //       phi_arr[tid].push_back(phi);
  //       phi_arr[tid].push_back(theta);
  //       // if(phi_arr[tid].size() == 0){
  //       //   phi_arr[tid].push_back(phi);
  //       //   theta_arr[tid].push_back(theta);
  //       //   s[tid]++;
  //       // }else{
  //       //   FLAG = true;
  //       //   for(i=0; i < phi_arr[tid].size(); i++){
  //       //     if((theta_arr[tid][i] == theta && phi_arr[tid][i] == phi) || (theta_arr[tid][i] == M_PI/2 && theta == M_PI/2) || (theta_arr[tid][i] == -M_PI/2 && theta == -M_PI/2)){
  //       //       FLAG = false;
  //       //       break;
  //       //     }
  //       //   }
  //       //   if(FLAG){
  //       //     phi_arr[tid].push_back(phi);
  //       //     theta_arr[tid].push_back(theta);
  //       //     s[tid] ++;
  //       //   }
  //       // }
  //       cnt++;
  //       s[tid]++;
  //     }
  //   }
  //   tmp_sum+= s[tid];
  }
  for(i=0; i < No_tile; i++){
    weight[i] = s[i] / tmp_sum;
  }
  return weight;
}
int* DecisionEngine::get_visible_tile(int* vp){
  int phi,theta;
  if(vp[0] < 0)
    phi = vp[0] + 360;
  else
    phi = vp[0];
  // convert from [-90; 90] to [0; 180]
  theta = vp[1] + 90;
  // printf("#[get_visible_tile]: (%d, %d) -> (%d, %d)\n", vp[0], vp[1], phi, theta);
  return vmask[phi][theta]; 
}
int* DecisionEngine::get_visible_pixel(int* vp){
  int phi,theta;
  if(vp[0] < 0)
    phi = vp[0] + 360;
  else
    phi = vp[0];
  // convert from [-90; 90] to [0; 180]
  theta = vp[1] + 90;
  // printf("#[getVisibleTile]: (%d, %d) -> (%d, %d)\n", vp[0], vp[1], phi, theta);
  return pixel[phi][theta]; 
}
double DecisionEngine::est_vp_psnr(double** TILE_MSE, int NO_TILE, int* tileVer, int* est_vp){
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
  for(int i=0; i < No_tile; i++){
    if(vmask[i] > 0){
      est_vp_mse += TILE_MSE[i][tileVer[i]] * pixel[i] * 1.0 / (vp_W * vp_H);
    }
  }
  v_psnr = 10 * log10((255*255)/est_vp_mse);
  return v_psnr;
}
double DecisionEngine::calc_distance(int* vp, int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H){
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
double DecisionEngine::get_Euclide_distance(double phi, double theta, double phi2, double theta2){
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
  // printf("#[get_Euclide_distance] dis: %.2f (x,y,z)=(%.2f,%.2f,%.2f) (x2,y2,z2)=(%.2f,%.2f,%.2f)\n",dist,x,y,z,x2,y2,z2);
  return dist;
}
void DecisionEngine::get_tile_center_pos(int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H, double* phi, double *theta){
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
      n = (tid == 0)?(H - face_H)/4:((H - face_H)*3/4 + face_H);
    }else{
      tile_id_v = (tid - 2) / No_tile_h;
      tile_id_h = (tid - 2) % No_tile_h;
      m = tile_id_h * tile_W + tile_W/2;
      n = (H - face_H)/2 + tile_id_v * tile_H + tile_H/2;
    }
    u = (m+0.5)/W;
    v = (n+0.5)/H;
    *phi = (u-0.5) * 2 * M_PI;
    *theta = (0.5 - v) * M_PI;
    // printf("#[get_tile_center_pos]: tid:%d m=%.0f n=%.0f phi=%.2f theta=%.2f\n", tid, m,n,*phi, *theta);
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
      // printf("#[get_tile_center_pos]: tid:%d phi=%.2f theta=%.2f tid_h=%d tid_v=%d\n", tid, *phi * 180.0 / M_PI, *theta * 180.0 / M_PI, tile_id_h, tile_id_v);
    }
  }
}
void DecisionEngine::write_result(int NO_SEG){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/log_tile_ver_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE_ID, BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/log_tile_ver_linear_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE_ID,BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/log_dec_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE_ID,BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/log_frame_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE_ID, BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vp_psnr, htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  }else{//DASH
    fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
    for(int ii=0; ii < NO_SEG; ii++){
      for(int k=0; k < INTER; k++){
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
      }
    }
    fflush(log_frame_psnr);
  }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_0(int NO_SEG, int phi, int theta, int DASH_VER){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
   // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedVP/log_tile_ver_DASH_VER_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, DASH_VER, phi,theta,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedVP/log_tile_ver_linear_DASH_VER_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedVP/log_dec_DASH_VER_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedVP/log_frame_DASH_VER_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii*INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_1(int NO_SEG, int phi, int theta, int BW){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_tile_ver_BW_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, BW, phi,theta,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_tile_ver_linear_BW_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_dec_BW_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_frame_BW_%d_phi_%d_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii*INTER +k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_2(int NO_SEG, int phi, int theta, int err, int point_id, int BW){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_tile_ver_BW_%d_phi_%d_theta_%d_err_%d_point_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, BW, phi,theta,err, point_id,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_tile_ver_linear_BW_%d_phi_%d_theta_%d_err_%d_point_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,err, point_id, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_dec_BW_%d_phi_%d_theta_%d_err_%d_point_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,err, point_id, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedVP/log_frame_BW_%d_phi_%d_theta_%d_err_%d_point_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,err, point_id,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii*INTER +k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_3(int NO_SEG, int phi, int theta, int speed, int DASH_VER){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedspeed/log_tile_ver_DASH_VER_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, DASH_VER, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedspeed/log_tile_ver_linear_DASH_VER_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedspeed/log_dec_DASH_VER_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/fixedspeed/log_frame_DASH_VER_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,DASH_VER, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_4(int NO_SEG, int phi, int theta, int speed, int BW){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_tile_ver_BW_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, BW, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_tile_ver_linear_BW_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_dec_BW_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_frame_BW_%d_phi_%d_theta_%d_speed_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_5(int NO_SEG, int phi, int theta, int speed, int BW){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_tile_ver_BW_%d_phi_%d_theta_%d_speed_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, BW, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_tile_ver_linear_BW_%d_phi_%d_theta_%d_speed_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_dec_BW_%d_phi_%d_theta_%d_speed_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/fixedspeed/log_frame_BW_%d_phi_%d_theta_%d_speed_theta_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,BW, phi,theta,speed,  TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_6(int NO_SEG, int HTRACE, int DASH_VER){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/realtrace/log_tile_ver_TRACE_%d_DASH_VER_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, HTRACE, DASH_VER, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/realtrace/log_tile_ver_linear_TRACE_%d_DASH_VER_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, DASH_VER, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/realtrace/log_dec_TRACE_%d_DASH_VER_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, DASH_VER,TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedDASH/realtrace/log_frame_TRACE_%d_DASH_VER_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, DASH_VER, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    printf("#[write_result] Cannot open log files\n");
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\n", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err[ii][0], est_err[ii][1], decide_width[ii], seg_br[ii]);
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_8(int NO_SEG, int HTRACE, int bwtrace_id){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/realtrace/realtrace/log_tile_ver_TRACE_%d_BWTRACE_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, HTRACE, bwtrace_id, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/realtrace/realtrace/log_tile_ver_linear_TRACE_%d_BWTRACE_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, bwtrace_id, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/realtrace/realtrace/log_dec_TRACE_%d_BWTRACE_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, bwtrace_id,TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/realtrace/realtrace/log_frame_TRACE_%d_BWTRACE_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, bwtrace_id, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    cout << "#[write_result] Cannot open log files\n";
    cout << fname << endl;
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi_1\test_theta_1\terr_1\test_phi_2\test_theta_2\terr_2\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      fprintf(log_dec, "%.2f\t", seg_thrp[ii]);
      fprintf(log_dec, "%.2f\t", buff_level[ii]);
      fprintf(log_dec, "%.2f\t", seg_down_time[ii]);
      fprintf(log_dec, "%.2f\t", time_to_next_request[ii]);
      fprintf(log_dec, "%d\t", last_frame_id[ii]);
      fprintf(log_dec, "%d\t%d\t", est_err[ii][0], est_err[ii][1]);
      fprintf(log_dec, "\n");
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t", est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err_ang[ii*INTER+k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t", est_frame_vp_2[ii*INTER + k][0], est_frame_vp_2[ii*INTER + k][1], est_err_ang_2[ii*INTER+k]);
        /*
        for(int j=0; j < NO_VER; j++)
          fprintf(log_frame_psnr, "%.2f\t", lowQLPixelPercent[ii*INTER+k][j]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%.2f\t", ext_tile[ii], useful_ext_tile[ii*INTER+k], useful_ext_tile_percent[ii*INTER+k] * 100, avg_visi_tile_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", s_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", avg_ext_tile_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", ext_tile_br[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", est_err_ang[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", visiTileOut_percent[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", visiTileOut_br[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", ext_tile_br_useful[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", useful_ext_tile_br[ii*INTER+k]);
        */
        fprintf(log_frame_psnr, "\n");
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
void DecisionEngine::write_result_7(int NO_SEG, int HTRACE, int BW){
  int jj=0; // frame id
  double BR = 0;
  int tid;
  char out_dir[] = "result";
  char fname[1024];
  FILE* log_tile_ver, *log_tile_ver_linear, *log_dec, *log_frame_psnr;
  // generate log files
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/realtrace/log_tile_ver_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H, HTRACE, BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_tile_ver = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/realtrace/log_tile_ver_linear_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_tile_ver_linear = fopen(fname, "w");
  // printf("#%s\n", fname);
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/realtrace/log_dec_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, BW,TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  log_dec = fopen(fname, "w");
  sprintf(fname, "%s/video_%s/%dx%d/fixedBW/realtrace/log_frame_TRACE_%d_BW_%d_METHOD_%d_INTER_%d_BUFF_%d_EST_%d.txt",out_dir,metadata->video_info.name.c_str(),W,H,HTRACE, BW, TILE_SELECT_METHOD, INTER, BUFF, VP_EST_METHOD);
  // printf("#%s\n", fname);
  log_frame_psnr = fopen(fname, "w");
  if(log_tile_ver == NULL || log_tile_ver_linear == NULL || log_dec == NULL){
    cout << "#[write_result] Cannot open log files\n";
    cout << fname << endl;
    exit(1);
  }
  // if(TILE_SELECT_METHOD > 1){
    fprintf(log_dec, "id\testThrp(kbps)\tbitrate(kbps)\tvp-bitrate(kbps)\tvp-psnr(dB)\tcalcTime(ms)\text_width\n");
    fprintf(log_frame_psnr, "fid\tdecid\test_vp_psnr\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\text_width\n");
    for(int ii=0; ii < NO_SEG; ii++){
      fprintf(log_tile_ver, "\nseg #%d calcTime: %d(ms)\n", ii, proc_time[ii]);
      for(int i=0; i < No_tile_v; i++){
        for(int j=0; j < No_tile_h; j++){
          fprintf(log_tile_ver, "%d ",tile_ver[ii][i * No_tile_h + j]);
          fprintf(log_tile_ver_linear, "%d ",tile_ver[ii][i * No_tile_h + j]);
        }
        fprintf(log_tile_ver, "\n");
      }
      fprintf(log_tile_ver_linear, "\n");
      fprintf(log_dec, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t", ii, est_seg_thrp[ii], seg_br[ii], vp_br[ii], seg_vp_psnr[ii], proc_time[ii]/1000.0, decide_width[ii]);
      fprintf(log_dec, "%.2f\t", seg_thrp[ii]);
      fprintf(log_dec, "%.2f\t", seg_down_time[ii]);
      fprintf(log_dec, "%.2f\n", buff_level[ii]);
      // calculate viewport psnr of frames in this interval
      // decide the adaptation result of this segment
      for(int k=0; k < INTER; k++){
        // printf("#frame #%d\n", ii * INTER + k);
        // double vp_psnr = est_vp_psnr(TILE_SEG_MSE[ii], No_tile, tile_ver[ii], htrace[ii * INTER + k]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%d\t", ii * INTER + k, ii, vpsnr[ii * INTER + k], htrace[ii*INTER +k][0], htrace[ii*INTER +k][1], est_frame_vp[ii*INTER + k][0], est_frame_vp[ii*INTER + k][1], est_err_frame[ii*INTER + k][0], est_err_frame[ii*INTER +k][1], decide_width[ii], seg_br[ii], visiTileOut[ii*INTER +k]);
        for(int j=0; j < NO_VER; j++)
          fprintf(log_frame_psnr, "%.2f\t", lowQLPixelPercent[ii*INTER+k][j]);
        fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%.2f\t", ext_tile[ii], useful_ext_tile[ii*INTER+k], useful_ext_tile_percent[ii*INTER+k] * 100, avg_visi_tile_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", s_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", avg_ext_tile_ver[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", ext_tile_br[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", est_err_ang[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", visiTileOut_percent[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", visiTileOut_br[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", ext_tile_br_useful[ii*INTER+k]);
        fprintf(log_frame_psnr, "%.2f\t", useful_ext_tile_br[ii*INTER+k]);
        fprintf(log_frame_psnr, "\n");
      }
      fflush(log_tile_ver);
      fflush(log_tile_ver_linear);
      fflush(log_dec);
      fflush(log_frame_psnr);
    }
  // }else{//DASH
  //   fprintf(log_frame_psnr, "fid\tdecid\tver\tbitrate\test_vp_psnr\n");
  //   for(int ii=0; ii < NO_SEG; ii++){
  //     for(int k=0; k < INTER; k++){
  //       fprintf(log_frame_psnr, "%d\t%d\t%.2f\t%d\t%.2f\n", ii * INTER + k, ii, DASH_SEG_PSNR[ii][tile_ver[ii][0]],tile_ver[ii][0], DASH_SEG_BR[ii][tile_ver[ii][0]]);
  //     }
  //   }
  //   fflush(log_frame_psnr);
  // }
  fclose(log_frame_psnr);
  fclose(log_tile_ver);
  fclose(log_tile_ver_linear);
  fclose(log_dec);
}
int* DecisionEngine::ISM_ext(int index){
  int NUM_EXT_PART = 9;
  int i[100], j, ii;
  int EXT_W = 2;
  double remainBW,remainBW2, max_vpPSNR = 0, vpPSNR, max_remainBW;
  double tmpsum[100];
  int* tileVer = (int*) malloc(No_tile * sizeof(int));
  int* ret = (int*) malloc(No_tile * sizeof(int));
  // 
  // printf("#[ISM_ext] (%d, %d)\n", est_vp[index][0], est_vp[index][1]);
  int* vmask = get_visible_tile(est_vp[index]);
  if(vmask == NULL){
    printf("(%d, %d)\n", est_vp[index][0], est_vp[index][1]);
    exit(1);
  }
  int* adapt_area = extVmask2(vmask, No_face, No_tile_h, No_tile_v, EXT_W);
  int count = 0, count2=0;
  double sum;
  int WITH_BREAK = 1;
  if(index < BUFF/INTER){
    for(j=0; j < No_tile; j++)
      ret[j] = 0;
    return ret;
  }
  //
  // printf("#[ISM_ext]: index=%d est_thrp=%.2f(kbps)\n", index, est_seg_thrp[index]);
  /* display the adapt area */
  /*
     for(j=0; j < No_tile; j++){
     printf("%d ", adapt_area[j]);
     if((j+1) % 8 == 0) printf("\n");
     }
     */
  // select the lowest version for background tiles
  remainBW = est_seg_thrp[index];
  for(j=0; j < No_tile; j++){
    if(adapt_area[j] == 0){
      tileVer[j] = 0;
      remainBW -= TILE_SEG_BR[index][j][0];
    }
  }
  // loop to find optimal versions of tiles in other areas
  for(i[0] = 0; i[0] < NO_VER; i[0]++){//1
    //
    tmpsum[0] = 0;      
    for(j=0; j < No_tile; j++){
      if(adapt_area[j] == 1){
        tmpsum[0] += TILE_SEG_BR[index][j][i[0]];
      }
    }
    if(tmpsum[0] >= remainBW) break;
    //
    for(i[1] = 0; i[1] < NO_VER; i[1]++){//2
      //
      tmpsum[1] = 0;
      for(j=0; j < No_tile; j++){
        if(adapt_area[j] == 2){
          tmpsum[1] += TILE_SEG_BR[index][j][i[1]];
        }
      }
      if(tmpsum[0] + tmpsum[1] >= remainBW) break;
      if(i[1] > i[0]) break;
      //
      for(i[2] = 0; i[2] < NO_VER; i[2]++){//3
        //
        tmpsum[2] = 0;
        for(j=0; j < No_tile; j++){
          if(adapt_area[j] == 3){
            tmpsum[2] += TILE_SEG_BR[index][j][i[2]];
          }
        }
        if(tmpsum[0] + tmpsum[1] + tmpsum[2] >= remainBW) break;
        if(i[2] > i[0]) break;
        //
        for(i[3] = 0; i[3] < NO_VER; i[3]++){//4
          //
          tmpsum[3] = 0;
          for(j=0; j < No_tile; j++){
            if(adapt_area[j] == 4){
              tmpsum[3] += TILE_SEG_BR[index][j][i[3]];
            }
          }
          if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] >= remainBW) break;
          if(i[3] > i[0]) break;
          //
          for(i[4] = 0; i[4] < NO_VER; i[4]++){//5
            //
            tmpsum[4] = 0;
            for(j=0; j < No_tile; j++){
              if(adapt_area[j] == 5){
                tmpsum[4] += TILE_SEG_BR[index][j][i[4]];
              }
            }
            if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] + tmpsum[4] >= remainBW) break;
            if(i[4] > i[0]) break;
            //
            for(i[5] = 0; i[5] < NO_VER; i[5]++){//6
              //
              tmpsum[5] = 0;
              for(j=0; j < No_tile; j++){
                if(adapt_area[j] == 6){
                  tmpsum[5] += TILE_SEG_BR[index][j][i[5]];
                }
              }
              if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] + tmpsum[4] + tmpsum[5] >= remainBW) break;
              if(i[5] > i[0] || i[5] > i[1] || i[5] > i[2] || i[5] > i[3] || i[5] > i[4]) break;
              //if(i[5] > i[0] || i[5] > i[1]) break;
              //
              for(i[6] = 0; i[6] < NO_VER; i[6]++){//7
                //
                tmpsum[6] = 0;
                for(j=0; j < No_tile; j++){
                  if(adapt_area[j] == 7){
                    tmpsum[6] += TILE_SEG_BR[index][j][i[6]];
                  }
                }
                if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] + tmpsum[4] + tmpsum[5] + tmpsum[6] >= remainBW) break;
                if(i[6] > i[0] || i[6] > i[1] || i[6] > i[2] || i[6] > i[3] || i[6] > i[4]) break;
                //if(i[6] > i[0] || i[6] > i[2]) break;
                //
                for(i[7] = 0; i[7] < NO_VER; i[7]++){//8
                  //
                  tmpsum[7] = 0;
                  for(j=0; j < No_tile; j++){
                    if(adapt_area[j] == 8){
                      tmpsum[7] += TILE_SEG_BR[index][j][i[7]];
                    }
                  }
                  if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] + tmpsum[4] + tmpsum[5] + tmpsum[6] + tmpsum[7] >= remainBW) break;
                  if(i[7] > i[0] || i[7] > i[1] || i[7] > i[2] || i[7] > i[3] || i[7] > i[4]) break;
                  //if(i[7] > i[0] || i[7] > i[3]) break;
                  //
                  for(i[8] = 0; i[8] < NO_VER; i[8]++){
                    //
                    tmpsum[8] = 0;
                    for(j=0; j < No_tile; j++){
                      if(adapt_area[j] == 9){
                        tmpsum[8] += TILE_SEG_BR[index][j][i[8]];
                      }
                    }
                    if(tmpsum[0] + tmpsum[1] + tmpsum[2] + tmpsum[3] + tmpsum[4] + tmpsum[5] + tmpsum[6] + tmpsum[7] + tmpsum[8] >= remainBW) break;
                    if(i[8] > i[0] || i[8] > i[1] || i[8] > i[2] || i[8] > i[3] || i[8] > i[4]) break;
                    //if(i[8] > i[0] || i[8] > i[4]) break;
                    count2++;
                    // assign version for each area
                    sum = 0;
                    for(j=0; j < No_tile; j++){
                      if(adapt_area[j] > 0){
                        tileVer[j] = i[adapt_area[j]-1];
                        sum += TILE_SEG_BR[index][j][tileVer[j]];
                      }
                    }
                    //
                    if(sum <= remainBW){
                      // using remaining bandwidth (if any) to increase tile versions
                      remainBW2 = remainBW - sum;
                      int BREAK_FLAG;
                      int area_id;
                      while(remainBW2 > 0){
                        BREAK_FLAG = 1;
                        //for(area_id = 5; area_id >= 2; area_id --){
                        for(ii=No_tile-1; ii >= 0; ii--){
                          if(adapt_area[ii] == 1  && tileVer[ii] < NO_VER -1 && (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]) < remainBW2){
                            remainBW2 = remainBW2 - (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]);
                            tileVer[ii] += 1;
                            BREAK_FLAG = 0;
                          }
                        }
                        //}
                        if(BREAK_FLAG == 1)
                          break;
                      }
                      while(remainBW2 > 0){
                        BREAK_FLAG = 1;
                        //for(area_id = 5; area_id >= 2; area_id --){
                        for(ii=0; ii < No_tile; ii++){
                          if(adapt_area[ii] >=2 && adapt_area[ii] <=5 && tileVer[ii] < NO_VER -1 && (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]) < remainBW2){
                            remainBW2 = remainBW2 - (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]);
                            tileVer[ii] += 1;
                            BREAK_FLAG = 0;
                          }
                        }
                        //}
                        if(BREAK_FLAG == 1)
                          break;
                      }
                      while(remainBW2 > 0){
                        BREAK_FLAG = 1;
                        for(area_id = 6; area_id <= 9; area_id ++){
                          for(ii=0; ii < No_tile; ii++){
                            if(adapt_area[ii] == area_id && tileVer[ii] < NO_VER -1 && (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]) < remainBW2){
                              remainBW2 = remainBW2 - (TILE_SEG_BR[index][ii][tileVer[ii] +1] - TILE_SEG_BR[index][ii][tileVer[ii]]);
                              tileVer[ii] += 1;
                              BREAK_FLAG = 0;
                            }
                          }
                        }
                        if(BREAK_FLAG == 1)
                          break;
                      }
                      // compute viewport PSNR
                      vpPSNR = 0;
                      int tmp_vp[2];
                      for(j=0; j < INTER; j++){
                        tmp_vp[0] = est_frame_vp[index * INTER + j][0] + est_err[index][0] * j * 1.0 / (INTER -1);
                        tmp_vp[1] = est_frame_vp[index * INTER + j][1] + est_err[index][1] * j * 1.0 / (INTER -1);
                        if(tmp_vp[0] > 180)
                          tmp_vp[0] -= 360;
                        if(tmp_vp[0] <= -180)
                          tmp_vp[0] += 360;
                        if(tmp_vp[1] > 90)
                          tmp_vp[1] = 90;
                        if(tmp_vp[1] < -90)
                          tmp_vp[1] = -90;
                        vpPSNR += 1.0 / INTER * est_vp_psnr(TILE_SEG_MSE[index], No_tile, tileVer, tmp_vp);
                      }
                      if(vpPSNR > max_vpPSNR){
                        //for(j=0; j < No_tile; j++)
                        //  ret[j] = tileVer[j];
                        std::copy (tileVer, tileVer + No_tile, ret);
                        max_vpPSNR = vpPSNR;
                        max_remainBW = remainBW2;
                        if(index==3 && max_vpPSNR >= 42.39){
                          double total_br = 0;
                          // printf("[DUC]: remainBW=%.2f\n", remainBW2);
                          for(j=0; j < No_tile; j++){
                            total_br += TILE_SEG_BR[index][j][ret[j]];
                            // printf("%d ", tileVer[j]);
                            // if((j+1) % 8 == 0) printf("\n");
                          }
                          // printf("[DUC]: total_br=%.2f remainBW=%.2f\n", total_br, est_seg_thrp[index] - total_br);
                        }
                      }
                      count++;
                    }
                    //
                  }//9
                }//8
              }//7
            }//6
          }//5
        }//4
      }//3
    }//2
  }//1
  /*
     printf("# select_version: index=%d vpPSNR=%.2f count=%d count2=%d remainBW=%.2f\n",index, max_vpPSNR, count, count2, max_remainBW);
     double tmp_bitrate = 0;
     for(ii=0; ii < No_tile; ii++){
     tmp_bitrate += TILE_SEG_BR[index][ii][ret[ii]];
     printf("%d ", ret[ii]);
     if((ii+1) % 8 == 0) printf("\n");
     }
     decide_width[index] = 2;
     printf("# bitrate=%.2f remainBW=%.2f\n", tmp_bitrate, est_seg_thrp[index] - tmp_bitrate);
     */
  return ret;
}
int* DecisionEngine::OPTIMAL(int index){
  int* tileVer = new int[No_tile];
  int* visi_pixel;
  int i,j,ii,jj;
  double* w = new double[No_tile];
  int LEN = 200;
  GRBVar** x;
  char*** x_name;
  double** d = TILE_SEG_MSE[index]; // tiles' distortions
  // calculate tiles' weights
  for(i=0; i < No_tile; i++){
    w[i] = 0;
  }
  if(index==0){
    for(i=0; i < No_tile; i++)
      tileVer[i] = 0;
    return tileVer;
  }
  for(j=0; j < INTER; j++){
    visi_pixel = get_visible_pixel(htrace[index*INTER+j]); 
    for(i=0; i < No_tile; i++){
      w[i] += visi_pixel[i] * 1.0 / (vp_W * vp_H);
    }
  }
  // 
  try {
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);
    // init
    x = new GRBVar*[No_tile];
    x_name = new char**[No_tile];
    for(i=0; i < No_tile; i++){
      x[i] = new GRBVar[NO_VER];
      x_name[i] = new char*[NO_VER];
      for(j=0; j < NO_VER; j++)
        x_name[i][j] = new char[LEN];
    }
    // Create variables
    for(i=0; i < No_tile; i++){
      for(j=0; j < NO_VER; j++){
        x_name[i][j] = new char[LEN];
        sprintf(x_name[i][j], "x_%d_%d", i, j);
        x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, x_name[i][j]);
      }
    }
    // objective funtion
    GRBQuadExpr obj;
    GRBQuadExpr avgVPQL = 0;
    GRBQuadExpr tmp = 0;
    for(i=0; i < No_tile; i++){
      for(j=0; j < NO_VER; j++){
        avgVPQL += w[i]*d[i][j]*x[i][j];
      }
    }
    obj = avgVPQL;
    model.setObjective(obj, GRB_MINIMIZE);
    // set constraints
    // c1: total tiles' bitrates does not exceed bandwidth
    GRBLinExpr constr_1 = 0;
    GRBLinExpr* constr_2 = new GRBLinExpr[No_tile];
    for(i=0; i < No_tile; i++){
      constr_2[i] = 0;
      for(j=0; j < NO_VER; j++){
        constr_1 += TILE_SEG_BR[index][i][j] * x[i][j];
        constr_2[i] += x[i][j];
      }
      model.addConstr(constr_2[i] == 1, "constr_2");
    }
    model.addConstr(constr_1 <= est_seg_thrp[index], "constr_1");
    // Optimize model
    model.optimize();
    for(i=0; i < No_tile; i++){
      for(j=0; j < NO_VER; j++){
        if(x[i][j].get(GRB_DoubleAttr_X) == 1){
          // cout << j << " ";
          tileVer[i] = j;
        }
      }
      // if((i+1)%No_tile_h == 0)
      //     cout << endl;
    }

    // cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  } catch(GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch(...) {
    cout << "Exception during optimization" << endl;
  }

  return tileVer;
}
template <typename T>
std::vector<int> DecisionEngine::sort_index(std::vector<T> const& values) {
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));
    std::sort(
        begin(indices), end(indices),
        [&](int a, int b) { return values[a] > values[b]; }
    );
    return indices;
}
