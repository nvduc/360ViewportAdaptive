
#ifndef __GA_COMMON_H__
#define __GA_COMMON_H__
#include <vector>
using namespace std;
long long tvdiff_us(struct timeval *tv1, struct timeval *tv2);
long long ga_usleep(long long interval, struct timeval *ptv);
int*** init3dArrayInt(int, int , int);
int** init2dArrayInt(int, int);
double*** init3dArrayDouble(int, int , int);
double** init2dArrayDouble(int, int);
double**** init4dArrayDouble(int, int , int, int);
int**** init4dArrayInt(int, int , int, int);
double avg(double* , int);
double sum(double* , int);
void get_face_tid(int No_face, int No_tile_h, int No_tile_v, int tid, int* face_id, int* tile_id);
void showArrayInt(int *arr, int N);
void showArrayDouble(double *arr, int N);
void showTileVersion(int* tile_ver, int No_tile_h, int No_tile_v);
void showTileInfo(double* info, int No_tile_h, int No_tile_v);
void getErrorPoint(int phi, int theta, int err, int point_id, int No_point, int* phi_est, int *theta_est);
FILE* open_file(char*);
double min(double, double);
double max(double, double);
int min(int, int);
int max(int, int);
// 
struct run_cfg
{
  std::vector<int> INTER_LIST;
  std::vector<int> BUFF_LIST;
  std::vector<int> BW_LIST;
  std::vector<int> METHOD_LIST;
  std::vector<int> HEADTRACE_LIST;
  std::vector<int> BWTRACE_LIST;
  std::vector<int> VP_EST_METHOD_LIST;
  std::vector<int> PHI;
  std::vector<int> THETA;
  std::vector<int> SPEED;
  std::vector<int> DASH_VER;
  std::vector<int> ERR;
  int BW_NUM;
  int METHOD_NUM;
  int INTER_NUM;
  int BUFF_NUM;
  int HEADTRACE_NUM;
  int VP_EST_METHOD_NUM;
  int NO_SEG;
  int VP_MODE;
  int DASH_VER_NUM;
  int ERR_NUM;
  int BWTRACE_NUM;
};
run_cfg load_run_cfg(const char* );
#endif



