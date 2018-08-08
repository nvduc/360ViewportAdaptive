#ifndef METADATA_H
#define METADATA_H
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;
struct Tiling{
  int No_face;			// Number of faces
  int face_W;				// Face's width
  int face_H; 			// Face's height
  int No_tile_h;		// Number of vertical tiles per faces
  int No_tile_v;		// Number of horizontal tiles per faces
  int No_tile;			// Number of tiles
  int tile_W;       // tile width
  int tile_H;       // tile height
  int FoV;          // Field of View (used to load visible mask)
  int vp_W;
  int vp_H;
  int*** vmask;			// Visible mask
  int*** pixel;			// Visible pixels
  double**** TILE_SEG_BR;
  double**** TILE_SEG_SIZE;
  double**** TILE_SEG_MSE;
  double**** TILE_SEG_PSNR;
  double**** TILE_FR_SIZE;
};
struct Headtrace{
  int No_user;
  int*** trace;
};
struct Video{
  string name;
  int NO_FRAME;         // Total number of session frames
  int NO_FRAME_ORIGIN;	// Total number of origin frames
  int W;				        // ERP's width
  int H;				        // ERP's height
  int FPS;			        // Frame rate
  int NO_VER;           // Number of quality versions
  int NO_INTER;         // Number of available intervals
  int* INTER_LIST;      // List of segment durations
  int NO_TILING;        // Number of tiling schemes
  Tiling* tile;         // Tile's info
  Headtrace htrace;     // head movement traces
};
class Metadata{
  public:
    Video video_info;
    Metadata();
    Metadata(const char* video_cfg);
    ~Metadata();
    void print_meta(void);
    int load_video_info(const char* video_cfg);
    int load_visible_mask(int TILING);
    int import_matrix_from_txt_file(const char* filename_X, vector <double>& v, int& rows, int& cols);
    int ReadNumbers( const string & s, vector <double> & v );
    int load_tile_info();
    int load_headtrace_info();
    int** est_head_trace(int htrace_id, int FRAME, int INTER, int BUFF);
    double calc_s_version(int* tileVer, double* tileWeight, int No_tile);
    double calc_avg_version(int* tileVer, double* tileWeight, int No_tile);
};
#endif
