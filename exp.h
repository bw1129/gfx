
/******* Graphics ID *******/

#define MOUSE_ID		4
#define MX_OFF			-9
#define MY_OFF			-6
#define FIXP_WIDTH		3

/****** Error Codes ********/

/****** Error Codes ********/

#define FIX_ERR			2
#define ERR_S_EARLY		4
#define ERR_S_LATE		8
#define ERR_DIRECTION  16

/****** Button Codes *******/

#define START_BUTTON	2
#define RESP_BUTTON		1

/****** Switches ***********/

#define RESET			1
#define NO_RESET		0

/****** BUFFERS ************/

#define MAX_TDOT		1000
#define MAX_COND		1024
#define MAX_STIM		1024

/******* Type Defs *********/

#ifndef EXP_TYPES
#define EXP_TYPES	1

typedef struct {
        short x;
        short y;
        short t;
        } DOTS;

typedef struct {
		int exp_number;	   // die Nummer des Experiments
		int number;	   // die Nummer der Vp
		int block;
		int mapping;
		int play;
		int task;
		char name[20];
		char info[20];
		char data_dir[50];
		char log_file[50];
		char cal_file[50];
		char edf_file[50];
		char thres_file[50];
		} ADRESS;

typedef struct {
		float lum;
		float cb;
		float tc;
		float lmean;
		float cont;
		float sigma;
		int size;
		int id;
		} GAUSSIAN;

typedef struct {
        float lmean;
        float cont;
		float lum_cont;
		float rg_cont;
		float by_cont;
        float sigma;
        float sp_f; //spatial frequency in degrees
        float te_f; //temporal frequency in frames
        float r; //orientation of grating
		float sp_phase; 
		float vel; 
		float direction; 
        int size; //size in pixel
		int tf_flag;
		int id;
        } GABOR;

typedef struct {
        float q;
        float chi2;
        float a;
        float b;
        } REGRESSION_P;

typedef struct {
        float x;
        float y;
		long time;
		short e;
        } DOT;

typedef struct {
		int type;
		long sttime, entime;
		} MY_EVT;

typedef struct {
        short target_w;
		short bull_w;
		int window_h;
		int stim_type[MAX_STIM];
		int color_type[MAX_STIM];
		int structured_bg;
        float monitor_freq;
		float ms_per_frame;
        float h_monitor_size;
        float v_monitor_size;
		float h_monitor_size_deg;
        float v_monitor_size_deg;
        float h_pixel;
        float v_pixel;
		int monitor_delay;
		int ramp_id;
		float direction;
		int velocity;
		int frames;
		float holdx0[MAX_STIM];
        float holdy0[MAX_STIM];
		float fpx;
        float fpy;
        float vp_distance;
        float v_deg2pix;
        float h_deg2pix;
		float vel[MAX_STIM];
		float te_f[MAX_STIM];
		float sp_f[MAX_STIM];
		float size[MAX_STIM];
		float wdth[MAX_STIM];
		float hght[MAX_STIM];
		float pos[MAX_STIM];
		float cont[MAX_STIM];
		float lum_cont[MAX_STIM];
		float rg_cont[MAX_STIM];
		float by_cont[MAX_STIM];
		int ori[MAX_STIM];
		int inner_radius[MAX_STIM];
		int disptime;		
		unsigned int video_flags;
		unsigned char black[3], white[3], red[3], green[3], blue[3], gray[3], tcolor[MAX_STIM][3], icolor[MAX_STIM][3], fpcolor[3], bgcolor[3];
		REGRESSION_P regx, regy;
        } GRAPHICS;

typedef struct {
		int exp_on;
        int blk;
        int trl;
		int num;
		float cont;
        int frames[MAX_TDOT];
		char img_type;
		int img_h, img_w, img_num;
		float img_cont;
		char out_file[200];
        int tx[MAX_STIM];
		int ty[MAX_STIM];
		int tf_flag[MAX_STIM];
        char msg[200];
	} TRIAL;
#endif


/************** Prototypes from gfx.c ***********************/

void GFX_enter_2Dmode(void);
void GFX_leave_2Dmode(void);
int GFX_SDL_init(GRAPHICS gp);
void GFX_clear_buffers(int r, int g, int b);
void GFX_check_for_graphics_errors(void);
void GFX_text(float x, float y, char *s, unsigned char *c, int center);
void GFX_rotating_target(int frames, int w);
void GFX_draw_target(GRAPHICS gp, float x, float y);
void GFX_clear_buffer_dkl(double lum, double rg, double by);
int GFX_refresh(void);
void GFX_triangle(int x, int y, int w, int dir, unsigned char *c);
void GFX_draw_gaussian(GRAPHICS gp, GAUSSIAN ga, int xs, int ys);
void GFX_circle_filled(int x, int y, int w, unsigned char *c);
void GFX_circle_line(int x, int y, int w, unsigned char *c);
void GFX_circle_filled_DKL(int x, int y, int w, double lum, double rg, double by);
void GFX_circle_line_DKL(int x, int y, int w, double lum, double rg, double by);
void GFX_line(int x1, int y1, int x2, int y2, unsigned char *c);
void GFX_draw_texture(float x, float y, int size, int id);
void GFX_get_gaussian(GRAPHICS gp, GAUSSIAN *ga);
void GFX_rect_filled(int x, int y, int x_w, int y_w, int ori, unsigned char *c);
void GFX_rect_filled_DKL(int x, int y, int x_w, int y_w, int ori, double lum, double rg, double by);
void GFX_rect(int x, int y, int w, int h, unsigned char *c);
int GFX_patch(GRAPHICS gp, GABOR *ga, float sp_ph, float te_ph, int sx, int sy, int sw, int sh);
int GFX_clear(GRAPHICS gp, GABOR *ga, int sx, int sy, int sw, int sh);
void GFX_lookups(void);
//unsigned char *GFX_readimage_jpg(char *fn, int *nr, int *nc, float cont, int gray_flag);
unsigned char *GFX_readimage_bmp(char *filename, int *nr, int *nc, float cont, int rgbflag);

/************** Prototypes from com.c ************************/

extern void COM_init(void) ;
extern int COM_clear(void) ;
extern int COM_get_touch(unsigned short *x, unsigned short *y) ;

/************** Prototypes from exp.c *************************/
void leave(char *str);

/************** Prototypes from util.c ************************/

void UTIL_shuffle(int condition_max, int condition[]);
void UTIL_calibrat_touch(GRAPHICS *gp, ADRESS *adr);
void UTIL_calibrate_EL2(GRAPHICS gp, ADRESS adr, int block, int exp_on);
int UTIL_time_delta(int reset);
void UTIL_init_adress(ADRESS *adr, int argc, char *argv[]);
void UTIL_print_time(FILE *fp);
void UTIL_do_color_calibration(void);
void IO_toggle_grab(void);
int IO_handle_event(SDL_Event *event, Uint8 type, SDLKey key, Uint16 *x, Uint16 *y, Uint8 *b);
void IO_check_escape(void);
int IO_check_key(SDLKey key);
void IO_flush_events(void);
void IO_wait_for_event(Uint8 type, SDLKey key, Uint16 *x, Uint16 *y, Uint8 *b);
SDLKey IO_get_key(void);
SDLKey IO_ask_key(void);


/************** Prototypes from color.c ***********************/
void COLOR_init_mon(char *m);
void COLOR_colRGB(double lum, double cb, double tc, double *c);	/* dkl [-0.5 0.5] -> rgb [0 1] */
void COLOR_igammacorrect(int r, int g, int b, int *R, int *G, int *B);
void COLOR_set_new_gamma(void);
void COLOR_set_old_gamma(void);

/************** Staircase **************************************/
void STAIR_updoml(TRIAL *tl, GRAPHICS gp, ADRESS adr, int ntrans,int (*stmrt)(TRIAL *tl, GRAPHICS gp, ADRESS adr, int i, double amp),int ncls,double start[],double ans[],int ncorrect);
