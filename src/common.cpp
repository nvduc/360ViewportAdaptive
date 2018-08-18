#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#ifndef WIN32
#ifndef ANDROID
#include <execinfo.h>
#endif /* !ANDROID */
#include <signal.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/syscall.h>
#endif /* !WIN32 */
#ifdef ANDROID
#include <android/log.h>
#endif /* ANDROID */
#ifdef __APPLE__
#include <syslog.h>
#endif

#if !defined(WIN32) && !defined(__APPLE__) && !defined(ANDROID)
#include <X11/Xlib.h>
#endif

#include "common.h"

#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
int* norm_vp_range(int* vp){
  int* norm_vp = new int[2];
  norm_vp[0] = vp[0];
  norm_vp[1] = vp[1];
  while(norm_vp[0] >= 180)
    norm_vp[0] -= 360;
  while(norm_vp[0] < -180)
    norm_vp[0] += 360;
  while(norm_vp[1] >= 90)
    norm_vp[1] -= 90;
  while(norm_vp[1] <= -90)
    norm_vp[1] += 90;
  return norm_vp;
}
double max(double a, double b){
  if(a >= b) return a;
  return b;
}
double min(double a, double b){
  if(a <= b) return a;
  return b;
}
int min(int a, int b){
  if(a <= b) return a;
  return b;
}
int max(int a, int b){
  if(a >= b) return a;
  return b;
}/**
 * Compute the time difference for two \a timeval data structure, i.e.,
 * \a tv1 - \a tv2.
 *
 * @param tv1 [in] Pointer to the first \a timeval data structure.
 * @param tv2 [in] Pointer to the second \a timeval data structure.
 * @return The difference in micro seconds.
 */
long long
tvdiff_us(struct timeval *tv1, struct timeval *tv2) {
	struct timeval delta;
	delta.tv_sec = tv1->tv_sec - tv2->tv_sec;
	delta.tv_usec = tv1->tv_usec - tv2->tv_usec;
	if(delta.tv_usec < 0) {
		delta.tv_sec--;
		delta.tv_usec += 1000000;
	}
	return 1000000LL*delta.tv_sec + delta.tv_usec;
}

/**
 * Sleep and wake up at \a ptv + \a interval (micro seconds).
 *
 * @param interval [in] The expected sleeping time (in micro seconds).
 * @param ptv [in] Pointer to the baseline time.
 * @return Currently always return 0.
 *
 * This function is useful for controlling precise sleeps.
 * We usually have to process each video frame in a fixed interval.
 * Each time interval includes the processing time and the sleeping time.
 * However, the processing time could be different in each iteration, so
 * the sleeping time has to be adapted as well.
 * To achieve the goal, we have to obtain the baseline time \a ptv
 * (using \a gettimeofday function)
 * \em before the processing task and call this function \em after
 * the processing task. In this case, the \a interval is set to the total
 * length of the interval, e.g., 41667 for 24fps video.
 *
 * This function sleeps for \a interval micro seconds if the baseline
 * time is not specified.
 */
long long
ga_usleep(long long interval, struct timeval *ptv) {
	long long delta;
	struct timeval tv;
	if(ptv != NULL) {
		gettimeofday(&tv, NULL);
		delta = tvdiff_us(&tv, ptv);
		if(delta >= interval) {
			usleep(1);
			return -1;
		}
		interval -= delta;
	}
	usleep(interval);
	return 0LL;
}
double**** init4dArrayDouble(int M, int N, int P, int Q){
  double**** ret;
  int i,j,k;
  ret = new double***[M];
  for(i=0; i < M; i++){
    ret[i] = new double**[N];
    for(j=0; j < N; j++){
      ret[i][j] = new double*[P];
      for(k=0; k < P; k++)
        ret[i][j][k] = new double[Q];
    }
  }
  return ret;
}
int**** init4dArrayInt(int M, int N, int P, int Q){
  int**** ret;
  int i,j,k;
  ret = new int***[M];
  for(i=0; i < M; i++){
    ret[i] = new int**[N];
    for(j=0; j < N; j++){
      ret[i][j] = new int*[P];
      for(k=0; k < P; k++)
        ret[i][j][k] = new int[Q];
    }
  }
  return ret;
}
int*** init3dArrayInt(int M, int N, int P){
	int*** ret;
	int i,j;
	ret = new int**[M];
	for(i=0; i < M; i++){
		ret[i] = new int*[N];
		for(j=0; j < N; j++)
			ret[i][j] = new int[P];
	}
	return ret;
}
int** init2dArrayInt(int M, int N){
	int** ret;
	int i,j;
	ret = new int*[M];
	for(i=0; i < M; i++){
		ret[i] = new int[N];
	}
	return ret;
}
double*** init3dArrayDouble(int M, int N, int P){
	double*** ret;
	int i,j;
	ret = new double**[M];
	for(i=0; i < M; i++){
		ret[i] = new double*[N];
		for(j=0; j < N; j++)
			ret[i][j] = new double[P];
	}
	return ret;
}
double** init2dArrayDouble(int M, int N){
	double** ret;
	int i,j;
	ret = new double*[M];
	for(i=0; i < M; i++){
		ret[i] = new double[N];
	}
	return ret;
}
double avg(double* a, int N){
	double avg=0;
	for(int i=0; i < N; i++)
		avg += a[i] / N;
	return avg;
}
double sum(double* a, int N){
	double sum=0;
	for(int i=0; i < N; i++)
		sum += a[i];
	return sum;
}
void get_face_tid(int No_face, int No_tile_h, int No_tile_v, int tid, int* face_id, int* tile_id){
	if(No_face == 1 || No_face == 6){
		*face_id = tid / (No_tile_h * No_tile_v);
		*tile_id = tid - (*face_id) * No_tile_h * No_tile_v;
	}
	if(No_face == 2){
		*face_id = 0;
		*tile_id = tid;
	}
}
void showArrayInt(int *arr, int N){
	int i;
	for(i=0; i < N; i++){
		printf("%d, ", arr[i]);
	}
	printf("\n");
}
void showArrayDouble(double *arr, int N){
	int i;
	for(i=0; i < N; i++){
		printf("%.2f, ", arr[i]);
	}
	printf("\n");
}
FILE* open_file(char* fileName){
	FILE* f = NULL;
	f = fopen(fileName, "w");
	if(f == NULL){
		printf("#[open_file] Cannot open file %s\n", fileName);
		exit(-1);
	}
}
void showTileVersion(int* tile_ver, int No_tile_h, int No_tile_v){
	int i,j;
	for(i=0; i < No_tile_v; i++){
		for(j=0; j < No_tile_h; j++){
			printf("%d ", tile_ver[i * No_tile_h + j]);
		}
		printf("\n");
	}
}
string showTileVersionInSeg(int** tile_ver, int NO_FRAME, int No_tile_h, int No_tile_v){
  std::stringstream ss;
  std::string out;
 int i,j,k, ii;
 bool FLAG = false;
 for(i=0; i < No_tile_v; i++){
   for(j=0; j < NO_FRAME; j++){
     for(k=0; k < No_tile_h; k++){
       FLAG = false;
       if(j>0){
          if(tile_ver[0][i*No_tile_h+k] == 0 && tile_ver[j][i*No_tile_h+k] == 1)
            FLAG = true;  
       }
       if(!FLAG)
         ss << tile_ver[j][i*No_tile_h+k] << " "; 
       else
         ss << "+" << tile_ver[j][i*No_tile_h+k] << " "; 
       //printf("%d ", tile_ver[j][i*No_tile_h+k]);
     }
     ss << "\t";
     //printf("\t");
   }
   ss << "\n";
   //printf("\n");
 }
 return ss.str();
}
void showTileInfo(double* tile_ver, int No_tile_h, int No_tile_v){
  int i,j;
  for(i=0; i < No_tile_v; i++){
    for(j=0; j < No_tile_h; j++){
      printf("%.5f ", tile_ver[i * No_tile_h + j]);
    }
    printf("\n");
  }
}
run_cfg load_run_cfg(const char* cfg){
  printf("#[load_run_cfg]:\n");
  string s;
  string delimiter = "=";
  string comment = "#";
  string key;
  string val_str;
  //
  string deli = ",";
  string token;
  size_t pos = 0;
  //
  double val;
  size_t pos_deli = 0;
  size_t pos_comm;
  ifstream infile(cfg);
  run_cfg config;
  if(infile == NULL){
   	exit(-1);
  }else{
    cout << "Reading file " << cfg << endl;
  }
  while(std::getline(infile, s)){
    // cout << s << endl;
    if((pos_deli=s.find(delimiter)) != std::string::npos){
      key = s.substr(0, pos_deli-1);
      pos_comm = s.find(comment);
      if(pos_comm != std::string::npos){
        val_str = s.substr(pos_deli + 2, pos_comm-pos_deli-2);
      }else{
        val_str = s.substr(pos_deli + 2, s.length());
      }
    }
    if(key.compare("NO_SEG")==0)  config.NO_SEG = (int) std::stod(val_str);
    if(key.compare("VP_MODE")==0)  config.VP_MODE = (int) std::stod(val_str);
    // printf("%s\n", val_str.c_str());
     if(key.compare("ERR")==0){
      printf("%sa\n", val_str.c_str());
      config.ERR_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.ERR.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.ERR_NUM ++;
        val_str.erase(0, pos + deli.length());
      }
    }
    // viewport
    if(key.compare("DASH_VER")==0){
      printf("%sa\n", val_str.c_str());
      config.DASH_VER_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.DASH_VER.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.DASH_VER_NUM ++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("TILE_SELECT_METHOD")==0){
      printf("%sa\n", val_str.c_str());
      config.METHOD_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.METHOD_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.METHOD_NUM ++;
        val_str.erase(0, pos + deli.length());
      }
      printf("#[load_run_info] no_method=%d\n", config.METHOD_NUM);
    }
    if(key.compare("BWTRACE")==0){
      // printf("%sa\n", val_str.c_str());
      config.BWTRACE_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.BWTRACE_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.BWTRACE_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
if(key.compare("BANDWIDTH")==0){
      // printf("%sa\n", val_str.c_str());
      config.BW_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.BW_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.BW_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("INTERVAL")==0){
      // printf("%sa\n", val_str.c_str());
      config.INTER_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.INTER_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.INTER_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("BUFFER")==0){
      // printf("%sa\n", val_str.c_str());
      config.BUFF_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.BUFF_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.BUFF_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("HEADTRACE")==0){
      // printf("%sa\n", val_str.c_str());
      config.HEADTRACE_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.HEADTRACE_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.HEADTRACE_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("VP_EST_METHOD")==0){
      // printf("%sa\n", val_str.c_str());
      config.VP_EST_METHOD_NUM = 0;
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.VP_EST_METHOD_LIST.push_back((int) std::stod(token));
        // cout << std::stod(token) << endl;
        config.VP_EST_METHOD_NUM++;
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("PHI")==0){
      // printf("%sa\n", val_str.c_str());
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.PHI.push_back((int) std::stod(token));
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("THETA")==0){
      // printf("%sa\n", val_str.c_str());
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.THETA.push_back((int) std::stod(token));
        val_str.erase(0, pos + deli.length());
      }
    }
    if(key.compare("SPEED")==0){
      // printf("%sa\n", val_str.c_str());
      while((pos=val_str.find(deli)) != std::string::npos){
        token = val_str.substr(0, pos);
        config.SPEED.push_back((int) std::stod(token));
        val_str.erase(0, pos + deli.length());
      }
    }
  }
  return config;
}
void getErrorPoint(int phi, int theta, int err, int point_id, int No_point, int* phi_est, int *theta_est){
  int tmp_phi_est;
  int tmp_theta_est;
  switch(point_id){
    case 1:
      tmp_phi_est = phi + err;
      tmp_theta_est = theta;
    break;
    case 2:
      tmp_phi_est = phi - err;
      tmp_theta_est = theta;
    break;
    case 3:
      tmp_phi_est = phi;
      tmp_theta_est = theta - err;
    break;
    case 4:
      tmp_phi_est = phi;
      tmp_theta_est = theta + err;
    break;
  }
  // 
  *phi_est = tmp_phi_est;
  if(tmp_phi_est <= -180)
    *phi_est += 360;
  if(tmp_phi_est > 180)
    *phi_est -= 360;
  // 
  *theta_est = tmp_theta_est;
  if(tmp_theta_est > 90)
    *theta_est = 90;
  if(tmp_theta_est < -90)
    *theta_est = -90;
}
