#include "DecisionEngine.h"
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
//
int main(int argc, char* argv[]){
  if(argc != 3){
    printf("usage: ./main video_cfg run_cfg\n");
    exit(1);
  }
  struct timeval t_start, t_end;
  const char* video_cfg = argv[1];
  const char* adapt_cfg = argv[2];
  int index, inter_id, buff_id, htrace_id, bw_id, method_id, vp_est_id, speed_id, dash_ver_id, bwtrace_id;
  int i, phi, theta, phi_id, theta_id, speed, bw, err, err_id, point_id, phi_est, theta_est;
  int No_point_per_err = 4;
  int* tile_ver;
  char buf[1024];
  FILE* fout;
  int phi_num = 360/15;
  int theta_num = 180/15 + 1;
  int speed_num = 90/5 + 1;
  int htrace_num = 400;
  double*** avgVPSNR = init3dArrayDouble(phi_num, theta_num,20);
  double**** avgVPSNR_case2 = init4dArrayDouble(phi_num, theta_num, speed_num, 20);
  double** avgVPSNR_case6 = init2dArrayDouble(htrace_num, 20);
  double tmp_phi, tmp_theta;
  const char* vname;
  // 
  AdaptInfo adaptInfo;
  //Metadata meta(video_cfg);
  Metadata* meta;
  meta = new Metadata(video_cfg);
  //
  run_cfg rcfg = load_run_cfg(adapt_cfg);
  //
  vname = meta->video_info.name.c_str();
  //
  adaptInfo.TILING = 1;
  for(vp_est_id=0; vp_est_id < rcfg.VP_EST_METHOD_NUM; vp_est_id ++){
    adaptInfo.VP_EST_METHOD = rcfg.VP_EST_METHOD_LIST[vp_est_id];
    switch(rcfg.VP_MODE){
      /* fixed viewport */
      // DASH-version
      case 0:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(dash_ver_id=0; dash_ver_id < rcfg.DASH_VER_NUM; dash_ver_id++){
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    // fix viewport
                    for(i=0; i < decEngine.NO_FRAME; i++){
                      decEngine.htrace[i][0] = phi;
                      // adaptLogic.metadata.trace.frame_vp[i][0] = phi;
                      decEngine.htrace[i][1] = theta;
                    }
                    /* run adaptation */
                    for(index = 0; index < rcfg.NO_SEG; index++){
                      // decEngine.seg_thrp[index] = adaptInfo.BW;
                      decEngine.seg_thrp[index] = decEngine.DASH_SEG_BR[index][rcfg.DASH_VER[dash_ver_id]] + 10;
                      tile_ver = decEngine.get_next_segment(index);
                    }
                    /* calculate  metrics */
                    decEngine.calc_result(rcfg.NO_SEG);
                    avgVPSNR[phi_id][theta_id][method_id] = decEngine.avgVPSNR;
                    /* write adaptation results to files */
                    decEngine.write_result_0(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, rcfg.DASH_VER[dash_ver_id]);
                    gettimeofday(&t_end, NULL);
                    printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              theta_id = 0;
              for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                sprintf(buf, "result/video_%s/%dx%d/fixedDASH/fixedVP/avg_fixed_vp_DASH_VER_%d_theta_%d.txt",vname, meta->video_info.W, meta->video_info.H, rcfg.DASH_VER[dash_ver_id], theta);
                fout = fopen(buf, "w");
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  fprintf(fout, "%d\t", phi);
                  for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                    fprintf(fout, "%.2f\t", avgVPSNR[phi_id][theta_id][method_id]);
                  }
                  fprintf(fout, "\n");
                  phi_id ++;
                }
                fclose(fout);
                theta_id++;
              }
            }
          }
        }      
        break;
        // fixed bandwidth
      case 1:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bw_id = 0; bw_id < rcfg.BW_NUM; bw_id ++){
              bw = rcfg.BW_LIST[bw_id];
              adaptInfo.BW = bw;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    // fix viewport
                    for(i=0; i < decEngine.NO_FRAME; i++){
                      decEngine.htrace[i][0] = phi;
                      // adaptLogic.metadata.trace.frame_vp[i][0] = phi;
                      decEngine.htrace[i][1] = theta;
                    }
                    /* run adaptation */
                    for(index = 0; index < rcfg.NO_SEG; index++){
                      // decEngine.seg_thrp[index] = adaptInfo.BW;
                      decEngine.seg_thrp[index] = bw;
                      tile_ver = decEngine.get_next_segment(index);
                    }
                    /* calculate  metrics */
                    decEngine.calc_result(rcfg.NO_SEG);
                    avgVPSNR[phi_id][theta_id][method_id] = decEngine.avgVPSNR;
                    /* write adaptation results to files */
                    decEngine.write_result_1(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, bw);
                    gettimeofday(&t_end, NULL);
                    printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              theta_id = 0;
              for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                sprintf(buf, "result/video_%s/%dx%d/fixedBW/fixedVP/avg_fixed_vp_BW_%d_theta_%d.txt",vname, meta->video_info.W, meta->video_info.H, bw, theta);
                fout = fopen(buf, "w");
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  fprintf(fout, "%d\t", phi);
                  for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                    fprintf(fout, "%.2f\t", avgVPSNR[phi_id][theta_id][method_id]);
                  }
                  fprintf(fout, "\n");
                  phi_id ++;
                }
                fclose(fout);
                theta_id++;
              }
            }
          }
        }
        break;
        // bandwidth traces
      case 2:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bw_id = 0; bw_id < rcfg.BW_NUM; bw_id ++){
              bw = rcfg.BW_LIST[bw_id];
              adaptInfo.BW = bw;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    for(err_id = 0; err_id < rcfg.ERR_NUM; err_id ++){
                      err = rcfg.ERR[err_id];
                      for(point_id = 1; point_id <= No_point_per_err; point_id++){
                        getErrorPoint(phi, theta, err, point_id, No_point_per_err, &phi_est, &theta_est);
                        for(i=0; i < decEngine.NO_FRAME; i++){
                          decEngine.htrace[i][0] = phi;
                          // adaptLogic.metadata.trace.frame_vp[i][0] = phi;
                          decEngine.htrace[i][1] = theta;
                        }
                        for(i=0; i < decEngine.NO_SEG; i++){
                          decEngine.est_err[i][0] = phi_est - phi;
                          decEngine.est_err[i][1] = theta_est - theta;
                        }
                        /* run adaptation */
                        for(index = 0; index < rcfg.NO_SEG; index++){
                          // decEngine.seg_thrp[index] = adaptInfo.BW;
                          decEngine.seg_thrp[index] = bw;
                          tile_ver = decEngine.get_next_segment(index);
                        }
                        /* calculate  metrics */
                        // fix viewport
                        for(i=0; i < decEngine.NO_FRAME; i++){
                          decEngine.htrace[i][0] = phi_est;
                          // adaptLogic.metadata.trace.frame_vp[i][0] = phi;
                          decEngine.htrace[i][1] = theta_est;
                        }
                        decEngine.calc_result(rcfg.NO_SEG);
                        avgVPSNR[phi_id][theta_id][method_id] = decEngine.avgVPSNR;
                        /* write adaptation results to files */
                        decEngine.write_result_2(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, err, point_id, bw);
                        gettimeofday(&t_end, NULL);
                        printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                      }
                    }
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              // theta_id = 0;
              // for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
              //   sprintf(buf, "result/video_%s/%dx%d/fixedBW/fixedVP/avg_fixed_vp_BW_%d_theta_%d.txt",vname, meta->video_info.W, meta->video_info.H, bw, theta);
              //   fout = fopen(buf, "w");
              //   phi_id = 0;
              //   for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
              //       fprintf(fout, "%d\t", phi);
              //       for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
              //         fprintf(fout, "%.2f\t", avgVPSNR[phi_id][theta_id][method_id]);
              //       }
              //       fprintf(fout, "\n");
              //       phi_id ++;
              //    }
              //    fclose(fout);
              //    theta_id++;
              // }
            }
          }
        }
        break;

        /* constant speed */
        // DASH-version
      case 3:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(dash_ver_id=0; dash_ver_id < rcfg.DASH_VER_NUM; dash_ver_id++){
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                vname = decEngine.metadata->video_info.name.c_str();
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    speed_id = 0;
                    for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
                      // fix viewport
                      for(i=0; i < decEngine.NO_FRAME; i++){
                        if(i==0){
                          decEngine.htrace[i][0] = phi;
                          tmp_phi = phi;
                        }
                        else{
                          tmp_phi += 1.0/decEngine.FPS * speed;
                          if(tmp_phi > 180)
                            tmp_phi -= 360;
                          decEngine.htrace[i][0] = (int) tmp_phi;
                        }
                        decEngine.htrace[i][1] = theta;
                      }
                      /* run adaptation */
                      for(index = 0; index < rcfg.NO_SEG; index++){
                        // decEngine.seg_thrp[index] = adaptInfo.BW;
                        decEngine.seg_thrp[index] = decEngine.DASH_SEG_BR[index][rcfg.DASH_VER[dash_ver_id]] + 10;
                        tile_ver = decEngine.get_next_segment(index);
                      }
                      /* calculate  metrics */
                      decEngine.calc_result(rcfg.NO_SEG);
                      avgVPSNR_case2[phi_id][theta_id][speed_id][method_id] = decEngine.avgVPSNR;
                      /* write adaptation results to files */
                      decEngine.write_result_3(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, speed, rcfg.DASH_VER[dash_ver_id]);
                      gettimeofday(&t_end, NULL);
                      printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                      speed_id ++;
                    }  
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              theta_id = 0;
              for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  speed_id = 0;
                  sprintf(buf, "result/video_%s/%dx%d/fixedDASH/fixedspeed/avg_fixed_speed_DASH_VER_%d_theta_%d_phi_%d.txt",vname, meta->video_info.W, meta->video_info.H, rcfg.DASH_VER[dash_ver_id], phi, theta);
                  printf("#[main] vname=%s\n", vname);
                  fout = fopen(buf, "w");
                  for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
                    fprintf(fout, "%d\t", speed);
                    for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                      fprintf(fout, "%.2f\t", avgVPSNR_case2[phi_id][theta_id][speed_id][method_id]);
                    }
                    fprintf(fout, "\n");
                    speed_id++;
                  }
                  fclose(fout);
                  phi_id ++;
                }
                theta_id++;
              }
            }
          }
        }      
        break;
        // fixed bandwidth
      case 4:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bw_id = 0; bw_id < rcfg.BW_NUM; bw_id ++){
              bw = rcfg.BW_LIST[bw_id];
              adaptInfo.BW = bw;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                vname = decEngine.metadata->video_info.name.c_str();
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    speed_id = 0;
                    for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
                      // fix viewport
                      for(i=0; i < decEngine.NO_FRAME; i++){
                        if(i==0){
                          decEngine.htrace[i][0] = phi;
                          tmp_phi = phi;
                        }
                        else{
                          tmp_phi += 1.0/decEngine.FPS * speed;
                          if(tmp_phi > 180)
                            tmp_phi -= 360;
                          decEngine.htrace[i][0] = (int) tmp_phi;
                        }
                        decEngine.htrace[i][1] = theta;
                      }
                      /* run adaptation */
                      for(index = 0; index < rcfg.NO_SEG; index++){
                        // decEngine.seg_thrp[index] = adaptInfo.BW;
                        decEngine.seg_thrp[index] = bw;
                        tile_ver = decEngine.get_next_segment(index);
                      }
                      /* calculate  metrics */
                      decEngine.calc_result(rcfg.NO_SEG);
                      avgVPSNR_case2[phi_id][theta_id][speed_id][method_id] = decEngine.avgVPSNR;
                      /* write adaptation results to files */
                      decEngine.write_result_4(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, speed, bw);
                      gettimeofday(&t_end, NULL);
                      printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                      speed_id ++;
                    }  
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              // theta_id = 0;
              // for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
              //   phi_id = 0;
              //   for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
              //       speed_id = 0;
              //       sprintf(buf, "result/video_%s/%dx%d/fixedBW/fixedspeed/avg_fixed_speed_BW_%d_theta_%d_phi_%d.txt",vname, meta->video_info.W, meta->video_info.H, bw, phi, theta);
              //       printf("#[main] vname=%s\n", vname);
              //       fout = open_file(buf);
              //       for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
              //         fprintf(fout, "%d\t", speed);
              //         for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
              //           fprintf(fout, "%.2f\t", avgVPSNR_case2[phi_id][theta_id][speed_id][method_id]);
              //         }
              //         fprintf(fout, "\n");
              //         speed_id++;
              //       }
              //       fclose(fout);
              //       phi_id ++;
              //    }
              //    theta_id++;
              // }
            }
          }
        }      
        break;
        // speed theta
      case 5:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bw_id = 0; bw_id < rcfg.BW_NUM; bw_id ++){
              bw = rcfg.BW_LIST[bw_id];
              adaptInfo.BW = bw;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                adaptInfo.HTRACE_ID = 0;
                DecisionEngine decEngine(meta, adaptInfo);
                vname = decEngine.metadata->video_info.name.c_str();
                phi_id = 0;
                for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
                  theta_id = 0;
                  for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
                    speed_id = 0;
                    for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
                      // fix viewport
                      for(i=0; i < decEngine.NO_FRAME; i++){
                        if(i==0){
                          decEngine.htrace[i][1] = theta;
                          tmp_theta = theta;
                        }
                        else{
                          tmp_theta += 1.0/decEngine.FPS * speed;
                          if(tmp_theta > 90)
                            tmp_theta = 90;
                          decEngine.htrace[i][1] = (int) tmp_theta;
                        }
                        decEngine.htrace[i][0] = phi;
                      }
                      /* run adaptation */
                      for(index = 0; index < rcfg.NO_SEG; index++){
                        // decEngine.seg_thrp[index] = adaptInfo.BW;
                        decEngine.seg_thrp[index] = bw;
                        tile_ver = decEngine.get_next_segment(index);
                      }
                      /* calculate  metrics */
                      decEngine.calc_result(rcfg.NO_SEG);
                      avgVPSNR_case2[phi_id][theta_id][speed_id][method_id] = decEngine.avgVPSNR;
                      /* write adaptation results to files */
                      decEngine.write_result_5(rcfg.NO_SEG, (phi >=0)?phi:(phi+360), theta, speed, bw);
                      gettimeofday(&t_end, NULL);
                      printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                      speed_id ++;
                    }  
                    theta_id ++;
                  }
                  phi_id ++;
                }
              }
              // theta_id = 0;
              // for(theta=rcfg.THETA[0]; theta <= rcfg.THETA[1]; theta += 15){
              //   phi_id = 0;
              //   for(phi=rcfg.PHI[0]; phi <= rcfg.PHI[1]; phi+=15){
              //       speed_id = 0;
              //       sprintf(buf, "result/video_%s/%dx%d/fixedBW/fixedspeed/avg_fixed_speed_BW_%d_theta_%d_phi_%d.txt",vname, meta->video_info.W, meta->video_info.H, bw, phi, theta);
              //       printf("#[main] vname=%s\n", vname);
              //       fout = open_file(buf);
              //       for(speed = rcfg.SPEED[0]; speed <= rcfg.SPEED[1]; speed += 5){ // phi
              //         fprintf(fout, "%d\t", speed);
              //         for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
              //           fprintf(fout, "%.2f\t", avgVPSNR_case2[phi_id][theta_id][speed_id][method_id]);
              //         }
              //         fprintf(fout, "\n");
              //         speed_id++;
              //       }
              //       fclose(fout);
              //       phi_id ++;
              //    }
              //    theta_id++;
              // }
            }
          }
        } 
        break;
        /* real head traces */
        // DASH-version
      case 6:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(dash_ver_id=0; dash_ver_id < rcfg.DASH_VER_NUM; dash_ver_id++){
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                for(htrace_id=0; htrace_id < rcfg.HEADTRACE_NUM; htrace_id++){
                  adaptInfo.HTRACE_ID = rcfg.HEADTRACE_LIST[htrace_id];
                  DecisionEngine decEngine(meta, adaptInfo);
                  vname = decEngine.metadata->video_info.name.c_str();
                  /* run adaptation */
                  for(index = 0; index < rcfg.NO_SEG; index++){
                    // decEngine.seg_thrp[index] = adaptInfo.BW;
                    decEngine.seg_thrp[index] = decEngine.DASH_SEG_BR[index][rcfg.DASH_VER[dash_ver_id]] + 10;
                    tile_ver = decEngine.get_next_segment(index);
                  }
                  /* calculate  metrics */
                  decEngine.calc_result(rcfg.NO_SEG);
                  avgVPSNR_case6[htrace_id][method_id] = decEngine.avgVPSNR;
                  /* write adaptation results to files */
                  decEngine.write_result_6(rcfg.NO_SEG, rcfg.HEADTRACE_LIST[htrace_id], rcfg.DASH_VER[dash_ver_id]);
                  gettimeofday(&t_end, NULL);
                  printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                }
              }
              for(htrace_id=0; htrace_id < rcfg.HEADTRACE_NUM; htrace_id++){
                sprintf(buf, "result/video_%s/%dx%d/fixedDASH/realtrace/avg_realtrace_DASH_VER_%d_htrace_%d.txt",vname, meta->video_info.W, meta->video_info.H, rcfg.DASH_VER[dash_ver_id], rcfg.HEADTRACE_LIST[htrace_id]);
                printf("#[main] vname=%s\n", vname);
                fout = fopen(buf, "w");                    
                fprintf(fout, "%d\t", rcfg.HEADTRACE_LIST[htrace_id]);
                for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                  fprintf(fout, "%.2f\t", avgVPSNR_case6[htrace_id][method_id]);
                }
                fprintf(fout, "\n");

              }
            }
          }
        }      
        break;
        // fixed bandwidth
      case 7:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bw_id = 0; bw_id < rcfg.BW_NUM; bw_id ++){
              bw = rcfg.BW_LIST[bw_id];
              adaptInfo.BW = bw;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                for(htrace_id=0; htrace_id < rcfg.HEADTRACE_NUM; htrace_id++){
                  adaptInfo.HTRACE_ID = rcfg.HEADTRACE_LIST[htrace_id];
                  DecisionEngine decEngine(meta, adaptInfo);
                  vname = decEngine.metadata->video_info.name.c_str();
                  /* run adaptation */
                  for(index = 0; index < rcfg.NO_SEG; index++){
                    // decEngine.seg_thrp[index] = adaptInfo.BW;
                    tile_ver = decEngine.get_next_segment(index);
                    // download tiles' versions
                    decEngine.seg_thrp[index] = bw;
                    decEngine.down_next_segment(index); 
                  }
                  /* calculate  metrics */
                  decEngine.calc_result(rcfg.NO_SEG);
                  avgVPSNR_case6[htrace_id][method_id] = decEngine.avgVPSNR;
                  /* write adaptation results to files */
                  decEngine.write_result_7(rcfg.NO_SEG, rcfg.HEADTRACE_LIST[htrace_id], bw);
                  gettimeofday(&t_end, NULL);
                  printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                }
              }
            }
          }
        }
        break;
        // bandwidth traces
      case 8:
        for(inter_id = 0; inter_id < rcfg.INTER_NUM; inter_id++){
          adaptInfo.INTER = rcfg.INTER_LIST[inter_id];
          for(buff_id = 0; buff_id < rcfg.BUFF_NUM; buff_id ++){
            adaptInfo.BUFF = rcfg.BUFF_LIST[buff_id];
            for(bwtrace_id = rcfg.BWTRACE_LIST[0]; bwtrace_id <= rcfg.BWTRACE_LIST[1]; bwtrace_id ++){
              adaptInfo.BWTRACE_ID = bwtrace_id;
              for(method_id = 0; method_id < rcfg.METHOD_NUM; method_id ++){
                adaptInfo.TILE_SELECT_METHOD = rcfg.METHOD_LIST[method_id];
                gettimeofday(&t_start, NULL);
                /* init decision engine */
                //for(htrace_id = rcfg.HEADTRACE_LIST[0]; htrace_id <= rcfg.HEADTRACE_LIST[1]; htrace_id++){
                //  adaptInfo.HTRACE_ID = htrace_id;
                for(htrace_id = 0; htrace_id < rcfg.HEADTRACE_NUM; htrace_id++){
                  adaptInfo.HTRACE_ID = rcfg.HEADTRACE_LIST[htrace_id];
                  DecisionEngine decEngine(meta, adaptInfo);
                  vname = decEngine.metadata->video_info.name.c_str();
                  /* run adaptation */
                  for(index = 0; index < rcfg.NO_SEG; index++){
                    // decEngine.seg_thrp[index] = adaptInfo.BW;
                    tile_ver = decEngine.get_next_segment(index);
                    // download tiles' versions
                    decEngine.down_next_segment_2(index); 
                  }
                  /* calculate  metrics */
                  decEngine.calc_result(rcfg.NO_SEG);
                  avgVPSNR_case6[htrace_id][method_id] = decEngine.avgVPSNR;
                  /* write adaptation results to files */
                  //decEngine.write_result_8(rcfg.NO_SEG, htrace_id, bwtrace_id);
                  decEngine.write_result_8(rcfg.NO_SEG, adaptInfo.HTRACE_ID, bwtrace_id);
                  gettimeofday(&t_end, NULL);
                  printf("#[method #%d] calc time: %.2f(ms)\n", method_id, tvdiff_us(&t_end, &t_start)/1000.0);
                }
              }
              }
            }
          }
      break;
    }
  }
  return 1;
}



