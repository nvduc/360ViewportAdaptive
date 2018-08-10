#include "Metadata.h"
#include <stdio.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "common.h"
#include <sys/types.h>
#include <cstring>
#include <dirent.h>
using namespace std;
#ifdef DEBUG
int main(int argc, char* argv[]){
  const char* video_cfg = argv[1];
  Metadata meta(video_cfg);
  int tileVer_ROI[] = {
4, 4, 1, 1, 1, 1, 1, 4,
4, 4, 1, 1, 1, 1, 4, 4,
4, 4, 1, 1, 1, 1, 4, 4,
3, 4, 1, 0, 0, 0, 0, 4,
3, 3, 0, 0, 0, 0, 0, 4,
0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 1, 1, 1, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0
};
  int tileVer_ISM[] = {
    0, 0, 4, 4, 4, 4, 0, 0, 
    0, 0, 4, 4, 4, 4, 0, 0, 
    0, 0, 4, 4, 4, 4, 0, 0, 
    0, 0, 4, 4, 5, 5, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0};
  int tileVer_Ghent[] = {
3, 3, 1, 1, 1, 1, 1, 3,
3, 3, 1, 0, 0, 1, 3, 3,
3, 3, 1, 0, 0, 1, 3, 3,
3, 3, 1, 0, 0, 1, 3, 3,
3, 3, 1, 0, 0, 1, 1, 3,
3, 1, 1, 0, 0, 0, 1, 3,
1, 1, 0, 0, 0, 0, 1, 1,
0, 0, 0, 0, 0, 0, 0, 0
};
  int tileVer_prob360[] = {
  1, 1, 0, 0, 0, 0, 1, 0,
  5, 3, 2, 0, 0, 1, 3, 4,
  5, 4, 3, 0, 0, 2, 3, 5,
  4, 4, 3, 0, 0, 2, 3, 5,
  2, 2, 0, 0, 0, 0, 1, 2,
  0, 0, 0, 0, 0, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0
  };
  int tileVer_BellLab[] = {
  4, 5, 1, 0, 0, 1, 2, 4,
  5, 4, 0, 0, 0, 0, 4, 4,
  4, 4, 0, 0, 0, 0, 3, 4,
  3, 4, 0, 0, 0, 0, 4, 4,
  3, 3, 0, 0, 0, 0, 0, 3,
  3, 0, 0, 0, 0, 0, 0, 5,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0
  };
  int index=3, tid, No_tile = 64, i, j;
  int INTER_ID = 0;
  int TILING  = 1;
  int vp[] = {-166 + 360, 16};
  int vp_W = 960;
  int vp_H = 960;
  int FRAME = 1792;
  int INTER = 32;
  int BUFF = INTER;
  int htrace_id = 1;
  int HTRACE_NUM = 300;
  double totalBR;
  std::vector<int> tileList = {15, 22};
  int* pixel = meta.video_info.tile[TILING].pixel[vp[0]][vp[1]+90];
  double* tileWeight = new double[No_tile];
  // for(i=0; i < tileList.size(); i++){
  //   tid = tileList[i];
  //   cout << "tid: " << tid << endl;
  //   for(j = 0; j < meta.video_info.NO_VER; j++){
  //     printf("%.2f\t", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][j]);
  //   }
  //   printf("\n");
  // }
  int** est_frame_vp;
  for(htrace_id = 0; htrace_id < HTRACE_NUM; htrace_id++){
//  for(htrace_id = 94; htrace_id < 95; htrace_id++){
    est_frame_vp = meta.est_head_trace(htrace_id, FRAME, INTER, BUFF);
  }
  exit(1);
  for(tid=0; tid < No_tile; tid++)
    tileWeight[tid] = pixel[tid] * 1.0 / (vp_W * vp_H);
  // // show tiles' psnr
  printf("#vmask\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%d ", meta.video_info.tile[TILING].vmask[vp[0]][vp[1]+90][tid]);
    if((tid+1)%8==0)
      printf("\n");
  }
  exit(1);
  printf("#tiles' weights\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.4f ", meta.video_info.tile[TILING].pixel[vp[0]][vp[1]+90][tid] * 1.0 / (vp_W * vp_H));
    if((tid+1)%8==0)
      printf("\n");
  }
   printf("#tiles' pixel:\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%d ",pixel[tid]);
    if((tid+1)%8==0)
      printf("\n");
  }
     // show tiles' psnr
  printf("#ROI-mse\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][tileVer_ROI[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }  
  printf("#BellLab-mse\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][tileVer_BellLab[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
    // show tiles' psnr
  printf("#Ghent-mse\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][tileVer_Ghent[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
  // show tiles' psnr
  printf("#prob360-mse\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][tileVer_prob360[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
    // show tiles' bitrates
  printf("#ROI-bitrate\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_BR[INTER_ID][index][tid][tileVer_ROI[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
    // show tiles' bitrates
  printf("#Ghent-bitrate\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_BR[INTER_ID][index][tid][tileVer_Ghent[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
  printf("#prob360-bitrate\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_BR[INTER_ID][index][tid][tileVer_prob360[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
    // show tiles' bitrates
  printf("#BellLab-bitrate\n");
  for(tid=0; tid < No_tile; tid++){
    printf("%.2f ", meta.video_info.tile[TILING].TILE_SEG_BR[INTER_ID][index][tid][tileVer_BellLab[tid]]);
    if((tid+1)%8==0)
      printf("\n");
  }
  printf("ROI: %.2f\nBellLab: %.2f\n Ghent: %.2f\nprob360: %.2f\n", meta.calc_s_version(tileVer_ROI, tileWeight, No_tile), meta.calc_s_version(tileVer_BellLab, tileWeight, No_tile), meta.calc_s_version(tileVer_Ghent, tileWeight, No_tile), meta.calc_s_version(tileVer_prob360, tileWeight, No_tile));
  printf("ROI: %.2f\nBellLab: %.2f\n Ghent: %.2f\nprob360: %.2f\n", meta.calc_avg_version(tileVer_ROI, tileWeight, No_tile), meta.calc_avg_version(tileVer_BellLab, tileWeight, No_tile), meta.calc_avg_version(tileVer_Ghent, tileWeight, No_tile), meta.calc_avg_version(tileVer_prob360, tileWeight, No_tile));
  for(j = 0; j < meta.video_info.NO_VER; j++){
    cout << (j+1) << " ";
    for(i=0; i < No_tile; i++){
      tid = i;
      printf("%.2f\t", meta.video_info.tile[TILING].TILE_SEG_MSE[INTER_ID][index][tid][j]);
    }
    printf("\n");
  }
} 
#endif
Metadata::Metadata(const char* video_cfg){
  int ret=load_video_info(video_cfg);
  if(ret != 0){
    printf("#[Metadata] Failed\n");
    exit(-1);
  }else{
    printf("#[Metadata] [OK] Data load finished!\n");
  }
}
Metadata::~Metadata(){
}
void Metadata::print_meta(void){
}
int Metadata::load_video_info(const char* video_cfg){
  int i;
  string s;
  string delimiter = "=";
  string comment = "#";
  string key;
  string val_str;
  string token;
  string deli = ",";
  double val;
  size_t pos_deli = 0;
  size_t pos_comm;
  size_t pos = 0;
  char buff[1024];
  ifstream infile(video_cfg);
  if(infile == NULL){
    cout << "#[load_video_info] Cannot open file " << video_cfg << endl;
    return -1;
  }else{
    cout << "#[load_video_info] Reading file" << video_cfg << endl;
  }
  while(std::getline(infile, s)){
    if((pos_deli=s.find(delimiter)) != std::string::npos){
      key = s.substr(0, pos_deli-1);
      pos_comm = s.find(comment);
      if(pos_comm != std::string::npos){
        val_str = s.substr(pos_deli + 2, pos_comm-pos_deli-2);
      }else{
        val_str = s.substr(pos_deli + 2, s.length());
      }
    }
    //printf("%s\n", val_str.c_str());
    // viewport
    if(key.compare("name")==0) video_info.name = val_str;
    if(key.compare("W")==0) video_info.W = (int) std::stod(val_str);
    if(key.compare("H")==0) video_info.H = (int) std::stod(val_str);
    if(key.compare("FPS")==0) video_info.FPS = (int) std::stod(val_str);
    if(key.compare("NO_FRAME")==0) video_info.NO_FRAME = (int) std::stod(val_str);
    if(key.compare("NO_FRAME_ORIGIN")==0) video_info.NO_FRAME_ORIGIN = (int) std::stod(val_str);
    if(key.compare("NO_HEAD_TRACE")==0) video_info.htrace.No_user = (int) std::stod(val_str);
    if(key.compare("NO_TILING")==0){
      video_info.NO_TILING = (int) std::stod(val_str);
      video_info.tile = new Tiling[video_info.NO_TILING];
      std::getline(infile, s);
      for(i=0; i < video_info.NO_TILING; i++){
        std::getline(infile, s);
        deli = "\t";
        sscanf(s.c_str(), "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d", &(video_info.tile[i].No_face),&(video_info.tile[i].face_W),&(video_info.tile[i].face_H),&(video_info.tile[i].No_tile_h),&(video_info.tile[i].No_tile_v), &(video_info.tile[i].FoV), &(video_info.tile[i].vp_W), &(video_info.tile[i].vp_H));
        //
        video_info.tile[i].No_tile = video_info.tile[i].No_face * video_info.tile[i].No_tile_h * video_info.tile[i].No_tile_v;
        video_info.tile[i].tile_W = video_info.tile[i].face_W / video_info.tile[i].No_tile_h;
        video_info.tile[i].tile_H = video_info.tile[i].face_H / video_info.tile[i].No_tile_v; 
        //
        printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", video_info.tile[i].No_face, video_info.tile[i].face_W, video_info.tile[i].face_H, video_info.tile[i].No_tile_h, video_info.tile[i].No_tile_v, video_info.tile[i].No_tile, video_info.tile[i].tile_W, video_info.tile[i].tile_H);
      } 
    }
    if(key.compare("NO_INTER")==0) video_info.NO_INTER = (int) std::stod(val_str);
    if(key.compare("NO_VER")==0) video_info.NO_VER = (int) std::stod(val_str);
    if(key.compare("INTER_LIST") == 0){
      video_info.INTER_LIST = new int[video_info.NO_INTER];
      i = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        video_info.INTER_LIST[i++] = (int) std::stod(token);
        val_str.erase(0, pos + deli.length());
      }
    }
  }
  cout <<"#[load_video_info] video.name:  " << video_info.name << endl;
  cout <<"#[load_video_info] video.W:  " << video_info.W << endl;
  cout <<"#[load_video_info] video.H:  " << video_info.H << endl;
  cout <<"#[load_video_info] video.FPS:  " << video_info.FPS << endl;
  cout <<"#[load_video_info] video.NO_TILING:  " << video_info.NO_TILING << endl;
  cout <<"#[load_video_info] video.NO_INTER:  " << video_info.NO_INTER << endl;
  for(i=0; i < video_info.NO_INTER; i++)
    printf("%d,", video_info.INTER_LIST[i]);
  printf("\n");
  cout <<"#[load_video_info] video.NO_VER:  " << video_info.NO_VER << endl;
  // load other information
  // visible mask
  for(i=0; i < video_info.NO_TILING; i++){
    if(load_visible_mask(i) != 0){
      return -1; 
    }; 
  }
  // load tiles' info
  if(load_tile_info() != 0){
    printf("Failed to load tile info\n");
    return -1; 
  }
  // load headtrace' info
  if(load_headtrace_info() != 0){
    printf("Failed to load head traces\n");
    return -1;
  }
}
int Metadata::load_visible_mask(int TILING){
  int i,j,k;
  vector <double> v;
  int rows;
  int cols;
  int phi_num = 360;
  int theta_num = 181;
  int ret;
  char buff[1024];
  int No_tile = video_info.tile[TILING].No_tile;
  int No_face = video_info.tile[TILING].No_face;
  sprintf(buff,"data/vmask/vmask_%dx%d_%d_face_%dx%d_FoV_%d.txt", video_info.W, video_info.H, video_info.tile[TILING].No_face, video_info.tile[TILING].No_tile_h, video_info.tile[TILING].No_tile_v, video_info.tile[TILING].FoV);
  printf("#[load_visible_mask]:\n");
  ret = import_matrix_from_txt_file(buff, v, rows, cols);
  if(ret != 0){
    printf("#[load_visible_mask] Cannot open file !\n");
    return -1;
  }
  int NO_VMASK_SAMPLE = rows;
  video_info.tile[TILING].vmask = init3dArrayInt(phi_num, theta_num, No_tile);
  video_info.tile[TILING].pixel = init3dArrayInt(phi_num, theta_num, No_tile);
  for(i=0; i < NO_VMASK_SAMPLE; i++){
    // printf("%d %d\n", int(v[i*cols]), int(v[i*cols + 1]) + 90);
    for(j=0; j < No_tile; j++){
      if(video_info.W == 3840) // ERP
        video_info.tile[TILING].pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = v[i * cols + 2*j + 3];
      if(video_info.W == 7680)
        video_info.tile[TILING].pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = v[i * cols + j + 2];
      if(video_info.tile[TILING].pixel[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] > 0)
        video_info.tile[TILING].vmask[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = 1;
      else
        video_info.tile[TILING].vmask[int(v[i*cols])][int(v[i*cols + 1]) + 90][j] = 0;
    }
  }
  return 0;
}
int Metadata::import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols){
  ifstream file_X;
  string line;
  // erase all current elements
  v.erase(v.begin(), v.end());
  cout << "#[import_matrix_from_txt_file] open text file: " << filename_X << endl;
  file_X.open(filename_X);
  if (file_X.is_open())
  {
    int i=0;
    getline(file_X, line);
    cols = ReadNumbers( line, v );
    rows = 1;
    // cout << "cols:" << cols << endl;
    while(getline(file_X, line) != 0){
      ReadNumbers(line, v);
      rows ++;
    }        
    file_X.close();
    return 0;
  }
  else{
    cout << "#[import_matrix_from_txt_file] Cannot open file " << filename_X << endl;
    return -1;
  }
}
int Metadata::ReadNumbers( const string & s, vector <double> & v ){
  istringstream is( s );
  double n;
  while( is >> n ) {
    v.push_back( n );
  }
  return v.size();
}
int Metadata::load_tile_info(){
  int INTER,i,ii,NO_SEG,tid,f,t;
  int rows,cols;
  int sid, vid, fid;
  vector<double> v;
  char buff[1024];
  for(ii=0; ii < video_info.NO_TILING; ii++){
    video_info.tile[ii].TILE_SEG_BR = new double***[video_info.NO_INTER];
    video_info.tile[ii].TILE_SEG_SIZE = new double***[video_info.NO_INTER];
    video_info.tile[ii].TILE_SEG_MSE = new double***[video_info.NO_INTER];
    video_info.tile[ii].TILE_SEG_PSNR = new double***[video_info.NO_INTER];
    video_info.tile[ii].TILE_FR_SIZE = new double***[video_info.NO_INTER];
    /* load segment's info of correspond to each interval */
    for(i=0; i < video_info.NO_INTER; i++){
      INTER = video_info.INTER_LIST[i];
      NO_SEG = video_info.NO_FRAME/INTER;
      // allocate memory
      // segment levels
      video_info.tile[ii].TILE_SEG_BR[i] = init3dArrayDouble(NO_SEG, video_info.tile[ii].No_tile, video_info.NO_VER);  
      video_info.tile[ii].TILE_SEG_SIZE[i] = init3dArrayDouble(NO_SEG, video_info.tile[ii].No_tile, video_info.NO_VER);  
      video_info.tile[ii].TILE_SEG_PSNR[i] = init3dArrayDouble(NO_SEG, video_info.tile[ii].No_tile, video_info.NO_VER);  
      video_info.tile[ii].TILE_SEG_MSE[i] = init3dArrayDouble(NO_SEG, video_info.tile[ii].No_tile, video_info.NO_VER);  
      // frame levels
      video_info.tile[ii].TILE_FR_SIZE[i] = init3dArrayDouble(video_info.NO_FRAME, video_info.tile[ii].No_tile, video_info.NO_VER);  
      for(tid=0; tid < video_info.tile[ii].No_tile; tid++){
        get_face_tid(video_info.tile[ii].No_face, video_info.tile[ii].No_tile_h, video_info.tile[ii].No_tile_v, tid, &f, &t);
        /* segment info */
        sprintf(buff,"data/tile/%s_%dx%d_%dfr_%dfps/tile_%dx%d/%dframe/f%d_t%d.txt",video_info.name.c_str(), video_info.W, video_info.H, video_info.NO_FRAME_ORIGIN, video_info.FPS,video_info.tile[ii].No_tile_h, video_info.tile[ii].No_tile_v, INTER, f,t);
        // sprintf(buff,"data/tile/%s_%dx%d_%dfr_%dfps/tile_%dx%d/%dframe/tile_%d.txt",video_info.name.c_str(), video_info.W, video_info.H, video_info.NO_FRAME_ORIGIN, video_info.FPS,video_info.tile[ii].No_tile_h, video_info.tile[ii].No_tile_v, INTER, t);
        if(import_matrix_from_txt_file(buff,v,rows,cols) != 0)
          return -1;
        for(sid=0; sid < video_info.NO_FRAME_ORIGIN/INTER; sid++){
          for(vid=0; vid < video_info.NO_VER; vid++){
            video_info.tile[ii].TILE_SEG_BR[i][sid][tid][vid] = v[sid * cols + 3*vid];
            video_info.tile[ii].TILE_SEG_PSNR[i][sid][tid][vid] = v[sid * cols + 3*vid + 1];
            video_info.tile[ii].TILE_SEG_MSE[i][sid][tid][vid] = v[sid * cols + 3*vid + 2];
            video_info.tile[ii].TILE_SEG_SIZE[i][sid][tid][vid] = video_info.tile[ii].TILE_SEG_BR[i][sid][tid][vid] * INTER * 1.0 / video_info.FPS; 
          }
        }
        /* repeat data to match the number of frames in a session */
        for(sid=video_info.NO_FRAME_ORIGIN/INTER; sid < NO_SEG; sid++){
          for(vid=0; vid < video_info.NO_VER; vid++){
            video_info.tile[ii].TILE_SEG_BR[i][sid][tid][vid] = video_info.tile[ii].TILE_SEG_BR[i][sid%(video_info.NO_FRAME_ORIGIN/INTER)][tid][vid]; 
            video_info.tile[ii].TILE_SEG_PSNR[i][sid][tid][vid] = video_info.tile[ii].TILE_SEG_PSNR[i][sid%(video_info.NO_FRAME_ORIGIN/INTER)][tid][vid]; 
            video_info.tile[ii].TILE_SEG_MSE[i][sid][tid][vid] = video_info.tile[ii].TILE_SEG_MSE[i][sid%(video_info.NO_FRAME_ORIGIN/INTER)][tid][vid]; 
            video_info.tile[ii].TILE_SEG_SIZE[i][sid][tid][vid] = video_info.tile[ii].TILE_SEG_SIZE[i][sid%(video_info.NO_FRAME_ORIGIN/INTER)][tid][vid]; 
          }
        }
        /* Frame info */
        sprintf(buff,"data/tile/%s_%dx%d_%dfr_%dfps/tile_%dx%d/%dframe/per_frame/f%d_t%d.txt",video_info.name.c_str(), video_info.W, video_info.H, video_info.NO_FRAME_ORIGIN, video_info.FPS,video_info.tile[ii].No_tile_h, video_info.tile[ii].No_tile_v, INTER, f, t);
        // sprintf(buff,"data/tile/%s_%dx%d_%dfr_%dfps/tile_%dx%d/%dframe/per_frame/tile_%d.txt",video_info.name.c_str(), video_info.W, video_info.H, video_info.NO_FRAME_ORIGIN, video_info.FPS,video_info.tile[ii].No_tile_h, video_info.tile[ii].No_tile_v, INTER, t);
        if(import_matrix_from_txt_file(buff,v,rows,cols) != 0)
          return -1;
        for(fid=0; fid < video_info.NO_FRAME_ORIGIN; fid++){
          for(vid=0; vid < video_info.NO_VER; vid++){
            video_info.tile[ii].TILE_FR_SIZE[i][fid][tid][vid] = v[fid * cols + 2*vid];
          }
        }
        /* repeat data to match the number of frames in a session */
        for(fid=video_info.NO_FRAME_ORIGIN; fid < video_info.NO_FRAME; fid++){
          for(vid=0; vid < video_info.NO_VER; vid++){
            video_info.tile[ii].TILE_FR_SIZE[i][fid][tid][vid] = video_info.tile[ii].TILE_FR_SIZE[i][fid%video_info.NO_FRAME_ORIGIN][tid][vid];
          }
        }
      }
    }
  }
  return 0;
}
int Metadata::load_headtrace_info(){
  char buff[1024], buff2[1024];
  int i,j,trace_id;
  vector <double> v;
  int rows, cols;
  sprintf(buff, "data/head_trace/");
  DIR* dirp = opendir(buff);
  struct dirent * dp;
  //
  video_info.htrace.trace = init3dArrayInt(video_info.htrace.No_user,video_info.NO_FRAME, 2);
  trace_id = 0;
  while ((dp = readdir(dirp)) != NULL) {
    printf("%s %d\n", dp->d_name, trace_id);
    if(strstr(dp->d_name, "xyz") == NULL)
      continue;
    sprintf(buff2,"%s%s", buff, dp->d_name);
    if(import_matrix_from_txt_file(buff2, v, rows, cols) != 0)
     return -1; 
    for(i=0; i < video_info.NO_FRAME; i++){
      for(j=0; j < cols; j++){
        video_info.htrace.trace[trace_id][i][j] = v[(i+1) * cols + j];
      }
    }
    trace_id ++;
  }
  closedir(dirp);
  return 0;
}
int** Metadata::est_head_trace(int htrace_id, int FRAME, int INTER,  int BUFF){
  int** htrace = video_info.htrace.trace[htrace_id];
  int** est_frame_vp = init2dArrayInt(FRAME, 2);
  int** est_err = init2dArrayInt(FRAME, 2);
  double* est_err_ang = new double[FRAME];
  int NO_SEG = FRAME / INTER;
  double* dist_ang = new double[NO_SEG];
  double* speed_ang = new double[NO_SEG];
  double* ang_speed = new double[NO_SEG];
  double** speed = init2dArrayDouble(NO_SEG, 2);
  int** cur_vp = init2dArrayInt(NO_SEG, 2);
  int index, i, j, delta_phi, delta_theta, last_frame_id;
  int VP_EST_WIN = INTER;
  char buff[1024];
  int FPS = 30;
  // 
  sprintf(buff, "data/head_trace/est/trace_%d_INTER_%d_BUFF_%d.txt", htrace_id, INTER, BUFF);
  FILE* fout = fopen(buff,"w");
  fprintf(fout, "frame_id\tphi\ttheta\test_phi\test_theta\terr_phi\terr_theta\tspeed\terr_ang\n");
  // 
  sprintf(buff, "data/head_trace/est/speed_%d_INTER_%d_BUFF_%d.txt", htrace_id, INTER, BUFF);
  FILE* fout2 = fopen(buff, "w");
  for(index = 0; index < NO_SEG; index ++){
   if(index < BUFF/INTER){
    for(i=0; i < INTER; i++){
      est_frame_vp[index * INTER + i][0] = 0;
      est_frame_vp[index * INTER + i][1] = 0;
      est_err_ang[index*INTER+i] =  acos(sin(htrace[index*INTER+i][1]*M_PI/180) * sin(est_frame_vp[index*INTER+i][1]*M_PI/180) + cos(htrace[index*INTER+i][1]*M_PI/180) * cos(est_frame_vp[index*INTER+i][1]*M_PI/180) * cos(abs(est_frame_vp[index*INTER+i][0] - htrace[index*INTER+i][0])*M_PI/180)) / M_PI * 180;

      fprintf(fout, "%d\t%d\t%d\t%d\t%d\t%.2f\n", index*INTER + i, htrace[index*INTER + i][0], htrace[index*INTER + i][1], est_frame_vp[index*INTER + i][0],est_frame_vp[index*INTER + i][1], est_err_ang[index*INTER + i]);
    }
    // calculate angular distance between two end-points of an interval
    dist_ang[index] = acos(sin(htrace[index*INTER+INTER-1][1]*M_PI/180) * sin(htrace[index*INTER][1]*M_PI/180) + cos(htrace[index*INTER+INTER-1][1]*M_PI/180) * cos(htrace[index*INTER][1]*M_PI/180) * cos(abs(htrace[index*INTER+INTER-1][0] - htrace[index*INTER][0])*M_PI/180)) / M_PI * 180;
    // calculate angular speed during this interval
    ang_speed[index] = dist_ang[index]/(INTER*1.0/FPS);
  }else{
    last_frame_id = index * INTER - BUFF;
    cur_vp[index][0] = htrace[last_frame_id][0];
    cur_vp[index][1] = htrace[last_frame_id][1];
    if(last_frame_id >= VP_EST_WIN){
      delta_phi = htrace[last_frame_id][0] - htrace[last_frame_id-VP_EST_WIN][0];
    }
    else{
      delta_phi = htrace[last_frame_id][0] - htrace[0][0];
    }
    if(delta_phi < -180)
      delta_phi += 360;
    else
      if(delta_phi > 180)
        delta_phi -= 360;
    if(last_frame_id == 0)
      speed[index][0] = 0;
    else
      speed[index][0] = delta_phi/(1.0 * ((last_frame_id >= VP_EST_WIN)?VP_EST_WIN:last_frame_id));
    // theta
    if(last_frame_id >= VP_EST_WIN)
      delta_theta = htrace[last_frame_id][1] - htrace[last_frame_id-VP_EST_WIN][1];
    else
      delta_theta = htrace[last_frame_id][1] - htrace[0][1];
    // printf("delta_theta=%d\n", delta_theta);
    if(delta_theta < -90)
      delta_theta += 180;
    else
      if(delta_theta > 90)
        delta_theta -= 180;
    if(last_frame_id == 0)
      speed[index][1] = 0;
    else
      speed[index][1] = delta_theta/(1.0 * ((last_frame_id >= VP_EST_WIN)?VP_EST_WIN:last_frame_id));
    // 
    ang_speed[index] = sqrt(speed[index][0] * speed[index][0] + speed[index][1] * speed[index][1]);
    // estimatt viewport 
    for(i=0; i < INTER; i++){
      est_frame_vp[index * INTER + i][0] = (int) (cur_vp[index][0] + speed[index][0] * (BUFF + i));
      est_frame_vp[index * INTER + i][1] = (int) (cur_vp[index][1] + speed[index][1] * (BUFF + i));
      // phi
      while(est_frame_vp[index * INTER + i][0] >= 180)
        est_frame_vp[index * INTER + i][0] -= 360;
      while(est_frame_vp[index * INTER + i][0] < -180)
        est_frame_vp[index * INTER + i][0] += 360;
      // theta
      if(est_frame_vp[index * INTER + i][1] >= 90)
        est_frame_vp[index * INTER + i][1] = 90;
      if(est_frame_vp[index * INTER + i][1] <= -90)//
        est_frame_vp[index * INTER + i][1] = -90;
      est_err_ang[index*INTER+i] =  acos(sin(htrace[index*INTER+i][1]*M_PI/180) * sin(est_frame_vp[index*INTER+i][1]*M_PI/180) + cos(htrace[index*INTER+i][1]*M_PI/180) * cos(est_frame_vp[index*INTER+i][1]*M_PI/180) * cos(abs(est_frame_vp[index*INTER+i][0] - htrace[index*INTER+i][0])*M_PI/180)) / M_PI * 180;

      fprintf(fout, "%d\t%d\t%d\t%d\t%d\t%.2f\n", index*INTER + i, htrace[index*INTER + i][0], htrace[index*INTER + i][1], est_frame_vp[index*INTER + i][0],est_frame_vp[index*INTER + i][1], est_err_ang[index*INTER + i]);
    }
    // calculate angular distance between two end-points of an interval
    dist_ang[index] = acos(sin(htrace[index*INTER+INTER-1][1]*M_PI/180) * sin(htrace[index*INTER][1]*M_PI/180) + cos(htrace[index*INTER+INTER-1][1]*M_PI/180) * cos(htrace[index*INTER][1]*M_PI/180) * cos(abs(htrace[index*INTER+INTER-1][0] - htrace[index*INTER][0])*M_PI/180)) / M_PI * 180;
    // calculate angular speed during this interval
    ang_speed[index] = dist_ang[index]/(INTER*1.0/FPS);
     }
    // 
    fprintf(fout2, "%d\t%.2f\n", index, ang_speed[index]);
  }
  fclose(fout);
  fclose(fout2);
  return est_frame_vp;
}
double Metadata::calc_s_version(int* tileVer, double* tileWeight, int No_tile){
  int tid;
  double thres = 0.1;
  double s_ver = 0;
  for(tid=0; tid < No_tile; tid++){
    if(tileWeight[tid] >= thres){
      s_ver += tileWeight[tid] * tileVer[tid];
    }
  }
  return s_ver;
}
double Metadata::calc_avg_version(int* tileVer, double* tileWeight, int No_tile){
  int tid;
  double thres = 0.1;
  double s_ver = 0;
  for(tid=0; tid < No_tile; tid++){
    if(tileWeight[tid] >= thres){
      s_ver += tileVer[tid];
    }
  }
  return s_ver;
}
