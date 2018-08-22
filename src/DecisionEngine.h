#include "Metadata.h"
struct AdaptInfo{
  int INTER;
  int BUFF;
  int TILING; 
  int HTRACE_ID;
  int VP_EST_METHOD;
  int TILE_SELECT_METHOD;
  int BW;
  int BWTRACE_ID;
};
struct Orient{
    double alpha;
    double beta;
    double gama;
};
class DecisionEngine{
  public:
    Metadata* metadata;
    /* Adaptation-related */
    int BW;       // network bandwidth (incase of constant bandwidth only)
    int INTER;    // adaptation interval
    int BUFF;     // client buffer size
    int TILE_SELECT_METHOD; 
    int TILING;   // tiling scheme
    int NO_SEG;   // number of adapattion interval
    int INTER_ID;
    int NO_VER;
    int No_tile;
    int No_tile_h;
    int No_tile_v;
    int NO_FRAME;
    int NO_FRAME_ORIGIN;
    int FPS;
    int No_face;
    int tile_W;
    int tile_H; 
    int face_W;
    int face_H;
    int vp_W;
    int vp_H;
    int W;
    int H;
    int FoV;
    int*** vmask;
    int*** pixel;
    double*** TILE_SEG_BR;
    double*** TILE_SEG_SIZE;
    double*** TILE_SEG_MSE;
    double*** TILE_SEG_PSNR;
    int No_user;
    int** htrace;
    double** bw_trace;
    int HTRACE_ID;
    int VP_EST_METHOD;
    double* VER_AVG_BR;
    // DASH
    double** DASH_SEG_BR;
    double** DASH_SEG_PSNR;
    double*** DASH_TILE_PSNR; 
    double*** DASH_TILE_MSE; 
    /* thoughput-related */
    double* seg_thrp;
    double* est_thrp_act;
    double* est_seg_thrp;
    /* buffer-related */
    int* last_frame_id;
    double* buff_level;
    double* stall_time;
    double* seg_down_time;
    double* seg_down_start_time;
    double* seg_down_finis_time;
    double* time_to_next_request;
    double B_0; // initial buffer size
    // prob360DASH
    double B_target = 0.8;
    double B_max = 0.8;
    // Ghent
    double B_cri = 2.0;
    bool BUFFERING_STATE = true;
    /* viewport-related */
    int* visiTileOut; // number of visible tiles not included in the visible tiles associated to the est. vp
    double* visiTileOut_percent; // number of visible tiles not included in the visible tiles associated to the est (pixel percentage)
    double* visiTileOut_br;
    double** lowQLPixelPercent; // percentage of pixel of each quality version in the viewport
    int** est_frame_vp;
    int** est_frame_vp_no_inpo;
    int** est_frame_vp_opp;
    int** est_frame_vp_2;
    int** est_frame_vp_2_opp;
    int** est_frame_vp_ins;
    int** est_frame_vp_ins_2;
    int** est_frame_vp_err_max_left;
    int** est_frame_vp_err_max_left_half;
    int** est_frame_vp_err_max_right;
    int** est_frame_vp_err_max_right_half;
    int** est_frame_vp_err_median_left;
    int** est_frame_vp_err_median_right;
    int** est_vp;
    int** est_vp_ins;
    int** est_err;
    int** est_err_ins;
    int** est_err_frame;
    int* max_est_err;
    int* max_est_err_frame;
    int* med_est_err;
    int* med_est_err_frame;
    std::vector<int> estErrList;
    std::vector<int> estErrListFrame;
    double* avg_est_err;
    int** cur_vp;
    double** speed;
    double** speed_ins;
    double* est_err_ang;
    double* est_err_ang_2;
    double* est_err_ang_ins;
    double* est_err_ang_ins_2;
    double* avg_visi_tile_ver;
    double* avg_ext_tile_ver;
    double* ext_tile_br;
    double* ext_tile_br_useful;
    double* s_ver;
    int* useful_ext_tile;
    double* useful_ext_tile_percent;
    double* useful_ext_tile_br;
    int* ext_tile;
    /* Adaptation results */
    int** tile_ver;   // tiles' versions
    int** tile_size;  // tiles' sizes
    int* decide_width;  // extension width
    int* proc_time;   // processing time
    double* vpsnr;    // frames' viewport psnr
    double* seg_br;   // segments' bitrates
    double* all_seg_br;   // segments' bitrates
    double* vp_br; // estimated segments' viewport bitrates 
    double* seg_vp_psnr; // segments' avg. viewport psnr
    //
    double* VPSNR_thres;
    double avgVPSNR;
    // BellLab
    double* ang_speed;
    double* net_delay;
    double*** Uti;          // utility
    double*** Cost;         // Cost
    /* function declarations */
    int* get_next_segment(int index);
    void down_next_segment(int index);
    void down_next_segment_2(int index);
    double calc_seg_down_time(double t_start, double br,double SD, double** bw_trace);
    void thrp_estimator(int index);
    void vp_estimator(int index);
    DecisionEngine(Metadata* meta, AdaptInfo adaptInfo);
    ~DecisionEngine();
    int* DASH(int index);
    int* EQUAL(int index);
    int* EXT_ALL(int index, int ext_width);
    int* EXT_ALL_same_ver(int index, int ext_width);
    int* ROI(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp);
    int* ROI_same_ver(double** TILE_BR, int NO_TILE, int NO_VER, int* vmask, double est_thrp);
    int* ISM(int index, int ext_width);
    int* ISM(int index);
    int* ISM_same_ver(int index, int ext_width);
    int* ISM_same_ver(int index);
    int* ISM_ext(int index);
    int* ISM_ext_quick(int index);
    int* Ghent(int index);
    int* Ireland(int index);
    int* Ireland_v2(int index);
    int* Ghent_opt(int index);
    int* BellLab(int index);
    // int* BellLab_v2(int index);
    int* ProbDASH(int index);
    int* OPTIMAL(int index);
    int* OPT_2_ONLY_PHI_AVG(int index);
    int* OPT_2_ONLY_PHI_INS(int index);
    int* OPT_2_ONLY_PHI_COMBI(int index);
    int* OPT_2_ONLY_PHI_ALL_AVG(int index);
    int* OPT_2_ONLY_PHI_ALL_INS(int index);
    int* OPT_2_ONLY_PHI_ERR_AVG(int index);
    int* OPT_2_ONLY_PHI_ERR_INS(int index);
    int* OPT_2_DYNAMIC(int index);
    int* OPT_2_LAST_ERR(int index);
    int* OPT_2_MAX_ERR(int index);
    int* OPT_2_NO_ERR(int index);
    int* OPT_2_NO_INPO(int index);
    int* OPT_2_ADAPTIVE(int index);
    int* OPT_2_ERR_TWO_SIDE(int index, int ext_err);
    int* PROPOSED_EXT_TWO_SIDES(int index, int ext_err);
    int* new_idea(int index);
    int* new_idea_half_max(int index);
    int* new_idea_new_idea(int index);
    int* EQUAL_THETA_0(int index);
    int* EQUAL_THETA_0_EXT(int index, int ext_width);
    // 
    int* get_visible_tile(int* vp);
    int* get_visible_pixel(int* vp);
    double est_vp_psnr(double** TILE_MSE, int NO_TILE, int* tileVer, int* est_vp);
    double calc_distance(int* vp, int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H);
    void get_tile_center_pos(int tid, int No_face, int No_tile_h, int No_tile_v, int face_W, int face_H, double* phi, double *theta);
    double get_Euclide_distance(double phi, double theta, double phi2, double theta2);
    void calc_result(int);
    char* get_method_name(int method_id);
    void write_result(int );
    void write_result_0(int, int, int, int);
    void write_result_1(int, int, int, int);
    void write_result_2(int, int, int, int, int, int);
    void write_result_3(int, int, int, int, int);
    void write_result_4(int, int, int, int, int);
    void write_result_5(int, int, int, int, int);
    void write_result_6(int, int, int);
    void write_result_7(int, int, int);
    void write_result_8(int, int, int);
    template <typename T>
    std::vector<int> sort_index(std::vector<T> const& values);
    double* calc_tile_view_prob(double, double, double, double, double, double);
    double* calc_tile_view_prob_2(int, double, double, double, double, double, double);
    double* calc_tile_weight(int, int , int , int);
    double calc_Ps(double phi, double theta, double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama);
    struct Orient* calc_orient_set(double phi, double theta, int* N);
    double calc_PE(struct Orient L,struct Orient L_est, double mean_alpha, double var_alpha, double mean_beta, double var_beta, double mean_gama, double var_gama);
    int calc_visiTileOut(int* vmask, int* vmask2);
    double calc_visiTileOut_percent(int index, int frame_id);
    double calc_visiTileOut_br(int index, int frame_id);
    double calc_lowQLPixelPercent(int* vmask, int* vmask2, int* tileVer, int VER);
    int* calc_required_ext_tile(int index, int* num);
    int calc_ext_tile(int index);
    double calc_ext_tile_avg_ver(int index);
    double calc_ext_tile_br(int index);
    double calc_ext_tile_br_useful(int index, int frame_id);
    double calc_avg_visi_tile_ver(int index, int frame_id);
    double calc_sversion(int index, int frame_id);
    int calc_useful_ext_tile(int index, int frame_id);
    double calc_useful_ext_tile_percent(int index, int frame_id);
    double calc_useful_ext_tile_br(int index, int frame_id);
};
