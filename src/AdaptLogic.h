#include "Metadata.h"

struct Throughput
{
	double* est_thrp_act;
	double* seg_thrp;
	double* est_seg_thrp;
};
struct Viewport2
{
	int** est_frame_vp;
	int** est_vp;
	int** est_err;
	int** cur_vp;
};
class AdaptLogic{
public:
	Metadata metadata;
	Throughput thrp;
	Viewport2 vp2;
	// Adaptation results
	int** tile_ver;		// tiles' versions
	int** tile_size;	// tiles' sizes
	int* decide_width;	// extension width
	int* proc_time;		// processing time
	double* vpsnr;		// frames' viewport psnr
	double* seg_br;		// segments' bitrates
	double* VPSNR_thres;// ATC papers
	double* avgVPSNR; 	// ATC papers
	//
	int* get_next_segment(int index);
	void thrp_estimator(int index);
	void vp_estimator(int index);
	AdaptLogic(Metadata meta);
	~AdaptLogic();
	int* DASH_ERP(int index);
	int* DASH_CMP(int index);
	int* EXT_ALL(int index, int ext_width);
	int* ROI(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp);
	int* ROI_adapt(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp, double alpha);
	int* ROI_adapt_2(double** TILE_BR, double** TILE_MSE, int NO_TILE, int NO_VER, int* vmask, double est_thrp, int** est_vp, double VPSNR_thres);
	int* ISM(int index, int ext_width);
	int* ISM(int index);
	int* Ghent(int index);
	int* Ireland(int index);
	int* test(int index);
	int* test2(int index);
	//int* test3(int index);
	// 
	int* get_visible_tile(int* vp);
	int* get_visible_pixel(int* vp);
	double est_vp_psnr(double** TILE_MSE, int NO_TILE, int* tileVer, int* est_vp);
	double calc_distance(int* vp, int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H);
	void get_tile_center_pos(int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H, double* phi, double *theta);
	double get_Euclide_distance(double phi, double theta, double phi2, double theta2);
	void calc_result();
	char* get_method_name(int method_id);
};
