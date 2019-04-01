#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <io.h>
#include <direct.h>
#include <wtypes.h>
#include <wingdi.h>
#include <time.h>
#include <winbase.h>

#include <iostream.h>

#include <EL2/eyelink.h>
#include <EL2/w32_exptsppt2.h>
#include <EL2/w32_demo.h>
//#include <TIFF/tiffio.h>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include <SDL/SDL_net.h>


#include "exp.h"
#include "maths.h"
#include <GLcommand.h>

static unsigned char *Limg;
static DISPLAYINFO dispinfo;

#define CLEAR_PATCH_CONSTANT 1024

int
	Break_it	= 0;
	Nrefreshdelay	= 5; //refreshes


float N_SD = 3;

union {
	unsigned char c[4];
	short s[2];
	float f;
} c2f;

float byte2float(char byte1, char byte2, char byte3, char byte4) {
	c2f.c[0]=byte1;
	c2f.c[1]=byte2;
	c2f.c[2]=byte3;
	c2f.c[3]=byte4;
	
	return (c2f.f);
}
unsigned char *float2byte(float flt) {
	c2f.f=flt;	
	return (c2f.c);
}

void draw_image(int x, int y, int w, int h, float pz)
{
	glPushMatrix();
	glRasterPos3f(x, y, 0.0);
	glPixelZoom(pz,pz);
	glDrawPixels(w,h,GL_RGB,GL_UNSIGNED_BYTE,Limg);
	glPopMatrix();
}
void GFX_draw_texture(float x, float y, int size, int id)
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture (GL_TEXTURE_2D, id);
	glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(x-size/2., y-size/2., 0.0f);
    glTexCoord2f(0.0, 1.0); glVertex3f(x+size/2., y-size/2., 0.0f);
    glTexCoord2f(1.0, 1.0); glVertex3f(x+size/2., y+size/2., 0.0f);
    glTexCoord2f(1.0, 0.0); glVertex3f(x-size/2., y+size/2., 0.0f);
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glFlush();
}

void leave(char *str)
{	
	exit(0);
}

int init_img(TRIAL *tl, GRAPHICS gp)
{
	int i,j, colflag;
	char buf[200];

	if (tl->img_type>3 || tl->img_type<1) tl->img_type=1;

	if(tl->img_type==1) {
		sprintf(buf,"C:\\src\\images\\1f-lum\\p%d.bmp",tl->img_num); 
		colflag=0; 
	}
	if(tl->img_type==2) {
		sprintf(buf,"C:\\src\\images\\1f-iso\\rg%d.bmp",tl->img_num);
		colflag=1;
	}
	if(tl->img_type==3) {
		sprintf(buf,"C:\\src\\images\\natural\\n1.bmp",tl->img_num);
		colflag=1;
	}

	if(Limg) free(Limg);
	Limg = GFX_readimage_bmp(buf, &tl->img_h, &tl->img_w, tl->img_cont, colflag);
	return 0;
}

void init_gfx(GRAPHICS *gp)
{
	gp->monitor_freq    = 60;
	gp->ms_per_frame	= 16.67;
	gp->monitor_delay	= (gp->ms_per_frame * Nrefreshdelay); //HGtv buffers ~5 refreshes at 60Hz
	gp->h_monitor_size  = 120.8;//cm horizontal
	gp->v_monitor_size  = 68;//cm vertical
	gp->video_flags		= SDL_OPENGL | SDL_NOFRAME;//  SDL_FULLSCREEN;
	gp->h_pixel         = 1920.0;
	gp->v_pixel         = 1080.0;
    gp->vp_distance		= 70.0;

	gp->h_monitor_size_deg = RAD2DEG * atan((gp->h_monitor_size*0.5)/gp->vp_distance) * 2;
	gp->v_monitor_size_deg = RAD2DEG * atan((gp->v_monitor_size*0.5)/gp->vp_distance) * 2;

	gp->h_deg2pix       = 1 / ( RAD2DEG * atan( (gp->h_monitor_size / gp->h_pixel) / gp->vp_distance )) ;
	gp->v_deg2pix       = 1 / ( RAD2DEG * atan( (gp->v_monitor_size / gp->v_pixel) / gp->vp_distance )) ;

	memset(gp->black, 0, 3);
	memset(gp->white, 255, 3);
	memset(gp->red, 0, 3); gp->red[0] = 255;
	memset(gp->green, 0, 3); gp->green[1] = 255;
	memset(gp->blue, 0, 3);	gp->blue[2] = 255;
	memset(gp->gray, 128, 3);
}

void set_intensity(GRAPHICS *gp, int intensity)
{
	memset(gp->black, 0, 3);
	memset(gp->white, intensity, 3); 
	memset(gp->red, 0, 3); gp->red[0] = intensity;
	memset(gp->green, 0, 3); gp->green[1] = intensity;
	memset(gp->blue, 0, 3);	gp->blue[2] = intensity;	
}

void print_start(GRAPHICS gp)
{
	char str[400], demo_str[2][100]={"Lab-B"};
	get_display_information(&dispinfo); 
	sprintf(str,"%s\n\n\%.1f Hz\n\%d bit\n%d x %d pixel\n%.2f x %.2f deg\n%.2f x %.2f pixel/deg\n\nExit with ESC", demo_str, dispinfo.refresh, dispinfo.bits, dispinfo.right+1, dispinfo.bottom+1, gp.h_monitor_size_deg, gp.v_monitor_size_deg, gp.h_deg2pix, gp.v_deg2pix);
	GFX_text(-50,50,str, gp.black,0);	
}

void draw_fix(GRAPHICS gp)
{
	GFX_circle_filled(gp.fpx,gp.fpy, gp.size[1], gp.fpcolor);
}
void erase_fix(GRAPHICS gp)
{
	if(gp.color_type[1]==SET_STIM_COLORS) GFX_circle_filled(gp.fpx,gp.fpy, gp.size[1], gp.bgcolor);
	if(gp.color_type[1]==SET_DKL) GFX_circle_filled_DKL(gp.fpx,gp.fpy, gp.size[1], 0.0, 0.0, 0.0);
}

void draw_target(TRIAL *tl, GRAPHICS gp, int id)
{
	char c1[3], c2[3];
	int i;
	for(i=0;i<3;i++){
		c1[i]=gp.tcolor[id][i];
		c2[i]=gp.icolor[id][i];
	}
	tl->tx[id]=round(gp.holdx0[id]);
	tl->ty[id]=round(gp.holdy0[id]);

	if(gp.stim_type[id]==DRAW_ANNULUS) {
		if(gp.color_type[id]==SET_STIM_COLORS){
			GFX_circle_filled(tl->tx[id],tl->ty[id], gp.size[id], c1);
			if (gp.inner_radius[id]) GFX_circle_line(tl->tx[id],tl->ty[id], gp.inner_radius[id], c2);
		}
		if(gp.color_type[id]==SET_DKL){
			GFX_circle_filled_DKL(tl->tx[id], tl->ty[id], gp.size[id], gp.lum_cont[id], gp.rg_cont[id], gp.by_cont[id]);
			if(gp.inner_radius[id]) GFX_circle_line_DKL(tl->tx[id],tl->ty[id], gp.inner_radius[id], gp.lum_cont[id], gp.rg_cont[id], gp.by_cont[id]);
		}
	}
	if(gp.stim_type[id]==DRAW_BAR) {
		if(gp.color_type[id]==SET_STIM_COLORS){
			GFX_rect_filled(tl->tx[id],tl->ty[id], gp.wdth[id], gp.hght[id], gp.ori[id], c1);
		}
		if(gp.color_type[id]==SET_DKL){
			GFX_rect_filled_DKL(tl->tx[id],tl->ty[id], gp.wdth[id], gp.hght[id], gp.ori[id], gp.lum_cont[id], gp.rg_cont[id], gp.by_cont[id]);		
		}	
	}
}

void erase_target(TRIAL *tl, GRAPHICS gp, int id)
{
	tl->tx[id]=round(gp.holdx0[id]);
	tl->ty[id]=round(gp.holdy0[id]);

	if(gp.stim_type[id]==DRAW_ANNULUS) {
		if(gp.color_type[id]==SET_STIM_COLORS) GFX_circle_filled(tl->tx[id],tl->ty[id], gp.size[id], gp.bgcolor);
		if(gp.color_type[id]==SET_DKL) GFX_circle_filled_DKL(tl->tx[id],tl->ty[id], gp.size[id], 0.0, 0.0, 0.0);
	}
	if(gp.stim_type[id]==DRAW_BAR) {
		if(gp.color_type[id]==SET_STIM_COLORS) GFX_rect_filled(tl->tx[id],tl->ty[id], gp.wdth[id], gp.hght[id], gp.ori[id], gp.bgcolor);
		if(gp.color_type[id]==SET_DKL) GFX_rect_filled_DKL(tl->tx[id],tl->ty[id], gp.wdth[id], gp.hght[id], gp.ori[id], 0.0, 0.0, 0.0);
	}	
}

void pursuit_target(TRIAL *tl, GRAPHICS gp, int id)
{
	char c1[3], c2[3];
	int i;

	for(i=0;i<3;i++){
		c1[i]=gp.tcolor[gp.ramp_id][i];
		c2[i]=gp.icolor[gp.ramp_id][i];
	}

	if(gp.stim_type[gp.ramp_id]==DRAW_ANNULUS) {
		if(gp.color_type[gp.ramp_id]==SET_STIM_COLORS){
			GFX_circle_filled(tl->tx[id],tl->ty[id], gp.size[gp.ramp_id], c1);
			if (gp.inner_radius[gp.ramp_id]) GFX_circle_line(tl->tx[id],tl->ty[id], gp.inner_radius[gp.ramp_id], c2);
		}
		if(gp.color_type[gp.ramp_id]==SET_DKL){
			GFX_circle_filled_DKL(tl->tx[id], tl->ty[id], gp.size[gp.ramp_id], gp.lum_cont[gp.ramp_id], gp.rg_cont[gp.ramp_id], gp.by_cont[gp.ramp_id]);
			if(gp.inner_radius[gp.ramp_id]) GFX_circle_line_DKL(tl->tx[id],tl->ty[id], gp.inner_radius[gp.ramp_id], gp.lum_cont[gp.ramp_id], gp.rg_cont[gp.ramp_id], gp.by_cont[gp.ramp_id]);
		}
	}
	if(gp.stim_type[gp.ramp_id]==DRAW_BAR) {
		if(gp.color_type[gp.ramp_id]==SET_STIM_COLORS){
			GFX_rect_filled(tl->tx[id],tl->ty[id], gp.wdth[gp.ramp_id], gp.hght[gp.ramp_id], gp.ori[gp.ramp_id], c1);
		}
		if(gp.color_type[gp.ramp_id]==SET_DKL){
			GFX_rect_filled_DKL(tl->tx[id],tl->ty[id], gp.wdth[gp.ramp_id], gp.hght[gp.ramp_id], gp.ori[gp.ramp_id], gp.lum_cont[gp.ramp_id], gp.rg_cont[gp.ramp_id], gp.by_cont[gp.ramp_id]);		
		}	
	}
}

void pursuit_target_clear(TRIAL *tl, GRAPHICS gp, int id)
{
	if(gp.stim_type[gp.ramp_id]==DRAW_ANNULUS) {
		if(gp.color_type[gp.ramp_id]==SET_STIM_COLORS) GFX_circle_filled(tl->tx[id],tl->ty[id], gp.size[gp.ramp_id], gp.bgcolor);
		if(gp.color_type[gp.ramp_id]==SET_DKL) GFX_circle_filled_DKL(tl->tx[id],tl->ty[id], gp.size[gp.ramp_id], 0.0, 0.0, 0.0);
	}
	if(gp.stim_type[gp.ramp_id]==DRAW_BAR) {
		if(gp.color_type[gp.ramp_id]==SET_STIM_COLORS) GFX_rect_filled(tl->tx[id],tl->ty[id], gp.wdth[gp.ramp_id], gp.hght[gp.ramp_id], gp.ori[gp.ramp_id], gp.bgcolor);
		if(gp.color_type[gp.ramp_id]==SET_DKL) GFX_rect_filled_DKL(tl->tx[id],tl->ty[id], gp.wdth[gp.ramp_id], gp.hght[gp.ramp_id], gp.ori[gp.ramp_id], 0.0, 0.0, 0.0);
	}		
}

void diode_signal(GRAPHICS gp, int intense)
{
	SDLKey key=0;
	set_intensity(&gp, intense);
	GFX_rect_filled(-940,-520, 50,50, 0, gp.white);
}

void calculate_position(TRIAL *tl, GRAPHICS gp, GABOR gb)
{
	int i, x, y;
	float displacement_per_frame, theta, xpix_per_frame, ypix_per_frame, step;

	if(tl->tf_flag[gb.id]==1 & gb.id>1) //STATIONARY
	{
		tl->tx[gb.id] = round(gp.holdx0[gb.id]);
		tl->ty[gb.id] = round(gp.holdy0[gb.id]);	
	}

	if((tl->tf_flag[gb.id]==2 | tl->tf_flag[gb.id]==3) & gb.id>1)//FLICKER OR DRIFT
	{		
		for(i=0;i<tl->frames[gb.id];i++) {			
			tl->tx[gb.id+i]=round(gp.holdx0[gb.id]);
			tl->ty[gb.id+i]=round(gp.holdy0[gb.id]);
		}
	}

	if(tl->tf_flag[gb.id]==4 & gb.id>1)//RAMP
	{
		displacement_per_frame=gb.vel/gp.monitor_freq;
		theta=gb.direction*DEG2RAD;//convert direction to radians
		//convert polar to cartesian coordinates, then convert back to degrees
		xpix_per_frame=(displacement_per_frame*cos(theta)) * gp.h_deg2pix; //r*cos(theta)
		ypix_per_frame=(displacement_per_frame*sin(theta)) * gp.v_deg2pix; //r*cos(theta)

		tl->tx[gb.id] = round(gp.holdx0[gb.id]);
		tl->ty[gb.id] = round(gp.holdy0[gb.id]);

		for(i=1;i<tl->frames[gb.id];i++) {	
			gp.holdx0[gb.id+i]=gp.holdx0[gb.id+(i-1)] + xpix_per_frame;
			gp.holdy0[gb.id+i]=gp.holdy0[gb.id+(i-1)] + ypix_per_frame;
			tl->tx[gb.id+i] =round(gp.holdx0[gb.id+i]);
			tl->ty[gb.id+i] =round(gp.holdy0[gb.id+i]);	
		}
	}
}

void calculate_ramp(TRIAL *tl, GRAPHICS gp, GABOR *gb)
{
	int i, x, y;
	float displacement_per_frame, theta, xpix_per_frame, ypix_per_frame;

	displacement_per_frame=gp.velocity/gp.monitor_freq;
	theta=gp.direction*DEG2RAD;//convert direction to radians
	//convert polar to cartesian coordinates, then convert back to degrees
	xpix_per_frame=(displacement_per_frame*cos(theta)) * gp.h_deg2pix; //r*cos(theta)
	ypix_per_frame=(displacement_per_frame*sin(theta)) * gp.v_deg2pix; //r*cos(theta)

	//starting point of ramp
	tl->tx[gp.ramp_id] = round(gp.holdx0[gp.ramp_id]);
	tl->ty[gp.ramp_id] = round(gp.holdy0[gp.ramp_id]);


	for(i=1;i<gp.frames;i++) {	
		gp.holdx0[gp.ramp_id+i]=gp.holdx0[gp.ramp_id+(i-1)] + xpix_per_frame;
		gp.holdy0[gp.ramp_id+i]=gp.holdy0[gp.ramp_id+(i-1)] + ypix_per_frame;
		tl->tx[gp.ramp_id+i] =round(gp.holdx0[gp.ramp_id+i]);
		tl->ty[gp.ramp_id+i] =round(gp.holdy0[gp.ramp_id+i]);		
	}

	//fixation point clear patch must be computer here for pursuit stim
	if(gp.structured_bg){
		if(gb->sigma==0) gb->sigma=0.2;
		GFX_clear(gp, gb, gp.fpx, gp.fpy, tl->img_w, tl->img_h, 1);
		//remaining clear stim for pursuit
		for(i=0;i<gp.frames;i++) 
		{		
			GFX_clear(gp, gb, tl->tx[gp.ramp_id+i], tl->ty[gp.ramp_id+i], tl->img_w, tl->img_h, gp.ramp_id+i);
		}
	}
}

int get_gauss_size (int x)
{
	int i, s=2; //max patch size=2^7pix (128pix diameter), constrained by img in GFX_patch
	for(i=0;i<7;i++) {
		if(x < s) return s;
		s = s*2;		
	}
	return 0;
}

#define G_COS	10		//accuracy of cosine lookup 1/G_COS in degrees
#define G_EXP	1000	//accuracy of exponential lookup 1/G_EXP
#define N_COS	(360*G_COS)
#define N_EXP	(5*G_EXP)

static float COS[N_COS];
static float EXP[N_COS];

void GFX_lookups(void)
{
	int i;
	for (i=0;i<N_COS;i++) {
		COS[i] = cos(i*PI/(float)(180*G_COS));
	}
	for (i=0;i<N_EXP;i++) {
		EXP[i] = exp(-i/(float)G_EXP);
	}
}

int GFX_patch(GRAPHICS gp, GABOR *gb, float sp_ph, float te_ph, int sx, int sy, int sw, int sh, int Patch_id)
{
 	int i, j, p, e, k, z, m, sxo, syo;
	float y, x, y0, x0, phase, expo, tmp;
	double l, c[3], lmax, p1, p2, lum=0, rg=0, by=0;
	GLubyte img[128*128*3];//match get_gauss_size

//	printf("Patch_id %d\n",Patch_id);
	gb->lmean=0.5;
	lmax = gb->cont * gb->lmean * te_ph;
	gb->size = gb->sigma * N_SD * gp.h_deg2pix; 
	gb->size = get_gauss_size(gb->size);
	if(!gb->size) printf("texture too large\n");

	x0 = (0.5 * gb->size) / gp.h_deg2pix;
	y0 = (0.5 * gb->size) / gp.v_deg2pix;
	p1 = cos(gb->r*DEG2RAD);
	p2 = sin(gb->r*DEG2RAD);	
	for(j=0;j<gb->size;j++) { //along the x-axis
		for(i=0;i<gb->size;i++) {
			x = j/gp.h_deg2pix;
			y = i/gp.v_deg2pix;
            phase = 360*gb->sp_f * ( (x-x0)*p1 + (y-y0)*p2 ) + sp_ph;
			expo = -( (x-x0)*(x-x0) + (y-y0)*(y-y0) ) / (gb->sigma*gb->sigma);
			p = abs(phase*G_COS+0.5) % N_COS;
			e = abs(expo*G_EXP+0.5);
			if(e >= N_EXP-1) e = N_EXP-1;
			l = lmax * (COS[p] * EXP[e]);
			lum=l * gb->lum_cont;
			rg=l * gb->rg_cont;
			by=l * gb->by_cont;
			COLOR_colRGB(lum, rg, by, c);
			k = j*gb->size*3 + i*3;
			
			img[k]=round(c[0]*255+0.5);
			img[k+1]=round(c[1]*255+0.5);
			img[k+2]=round(c[2]*255+0.5);
	
			sxo = sw/2. + sx + j;
			syo = sh/2. + sy + i;
			m = (syo-gb->size/2.) * 3 * sw + (sxo-gb->size/2.) * 3;
			k = j*gb->size*3 + i*3;
			c[0] = (Limg[m]  /255.) + c[0] - 0.5;
			c[1] = (Limg[m+1]/255.) + c[1] - 0.5;
			c[2] = (Limg[m+2]/255.) + c[2] - 0.5;
			for(z=0;z<3;z++) {
				if(c[z]<0) {
					c[z]=0.0;
				//	printf("overflow \n");
				}
				if(c[z]>1) {
				 c[z]=1.0;
				//	printf("overflow \n");
				}
				img[k+z] = (255*c[z]+0.5);
			}
        }
	}	
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glBindTexture (GL_TEXTURE_2D, Patch_id);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, gb->size, gb->size, 0, GL_RGB, GL_UNSIGNED_BYTE, img);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	return 0;
}

int GFX_clear(GRAPHICS gp, GABOR *gb, int sx, int sy, int sw, int sh, int Patch_id)
{
 	int i, j, p, e, k, z, m, sxo, syo;
	float y, x;
	double c[3],lmax;
	GLubyte img[128*128*3];//match get_gauss_size

	Patch_id=Patch_id+CLEAR_PATCH_CONSTANT;

	gb->size = gb->sigma * N_SD * gp.h_deg2pix; 
	gb->size = get_gauss_size(gb->size);
	if(!gb->size) printf("texture too large\n");
	
	for(j=0;j<gb->size;j++) { //along the x-axis
		for(i=0;i<gb->size;i++) {      
			k = j*gb->size*3 + i*3;	
			sxo = sw/2. + sx + j;
			syo = sh/2. + sy + i;
			m = (syo-gb->size/2.) * 3 * sw + (sxo-gb->size/2.) * 3;
			k = j*gb->size*3 + i*3;
			c[0] = (Limg[m]  /255.);
			c[1] = (Limg[m+1]/255.);
			c[2] = (Limg[m+2]/255.);
			for(z=0;z<3;z++) {
				if(c[z]<0) c[z]=0;
				if(c[z]>1) c[z]=1;
				img[k+z] = (255*c[z]+0.5);		
			}
        }
	}
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glBindTexture (GL_TEXTURE_2D, Patch_id);
    glTexImage2D(GL_TEXTURE_2D, 0, 3, gb->size, gb->size, 0, GL_RGB, GL_UNSIGNED_BYTE, img);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
	return 0;
}

void get_fix_clear_patch(TRIAL *tl, GRAPHICS gp, GABOR *gb)
{
	//F-OFF PATCH
	if(gb->sigma==0) gb->sigma=0.2;
	GFX_clear(gp, gb, gp.fpx, gp.fpy, tl->img_w, tl->img_h, 1); //gb->id = 1 always for fixation
}

void get_patches(TRIAL *tl, GRAPHICS gp, GABOR *gb)
{
	int i;
	float te_ph, offset, phase=0;
	int sp_ph;
	int dir[2]={1,-1};

	//STATIONARY
	if(tl->tf_flag[gb->id]==1 & gb->id > 1) 
	{
		te_ph = 1;
		sp_ph = gb->sp_phase;
		GFX_patch(gp, gb, sp_ph, te_ph, tl->tx[gb->id], tl->ty[gb->id], tl->img_w, tl->img_h, gb->id);
		GFX_clear(gp, gb, tl->tx[gb->id], tl->ty[gb->id], tl->img_w, tl->img_h, gb->id);
	}

	//COUNTERPHASE FLICKERING 
	if(tl->tf_flag[gb->id]==2 & gb->id > 1)
	{
		sp_ph = gb->sp_phase;
		offset = PI;
		for(i=0;i<tl->frames[gb->id];i++) 
		{
			te_ph = sin(2*PI*gb->te_f*(i*.016667) + offset); //sec per frame
			GFX_patch(gp, gb, sp_ph, te_ph, tl->tx[gb->id], tl->ty[gb->id], tl->img_w, tl->img_h, gb->id+i);
			GFX_clear(gp, gb, tl->tx[gb->id], tl->ty[gb->id], tl->img_w, tl->img_h, gb->id+i);
		}
	}

	// DRIFTING
	if(tl->tf_flag[gb->id]==3 & gb->id > 1)
	{
		te_ph=1;
		phase=gb->sp_phase*(180/PI);//starting phase
		for(i=0;i<tl->frames[gb->id];i++) 
		{					
			phase=phase+(gb->vel*((PI*2)/gp.monitor_freq)*gb->sp_f);
			if(phase < 0) phase=phase+(2*PI);
			if(phase > (2*PI)) phase=phase-(2*PI);
			sp_ph=phase*RAD2DEG;
			GFX_patch(gp, gb, sp_ph, te_ph, tl->tx[gb->id+i], tl->ty[gb->id+i], tl->img_w, tl->img_h, gb->id+i);
			GFX_clear(gp, gb, tl->tx[gb->id+i], tl->ty[gb->id+i], tl->img_w, tl->img_h, gb->id+i);
		}
	}
	// RAMP
	if(tl->tf_flag[gb->id]==4 & gb->id > 1)
	{
		te_ph=1;
		sp_ph = gb->sp_phase;
		for(i=0;i<tl->frames[gb->id];i++) 
		{	
			GFX_patch(gp, gb, sp_ph, te_ph, tl->tx[gb->id+i], tl->ty[gb->id+i], tl->img_w, tl->img_h, gb->id+i);
			GFX_clear(gp, gb, tl->tx[gb->id+i], tl->ty[gb->id+i], tl->img_w, tl->img_h, gb->id+i);
		}
	}
}

void display_patch(TRIAL *tl, GRAPHICS gp, GABOR gb, int id)
{
	GFX_draw_texture(tl->tx[id],tl->ty[id], gb.size, id);
}

void clear_patch(TRIAL *tl, GRAPHICS gp, GABOR gb, int id)
{	
	if(id==1) GFX_draw_texture(gp.fpx, gp.fpy, gb.size, id+CLEAR_PATCH_CONSTANT);//f-off patch	
	if(id>1) GFX_draw_texture(tl->tx[id],tl->ty[id], gb.size, id+CLEAR_PATCH_CONSTANT);//t-off patch
}

int main(int argc, char *argv[])
{	
	int i, p, j, r = 0, t=0, snum=0, n=0, bg=0, d=0;
	int Diode_on = 1;
	char str[200];
	unsigned short index;
	unsigned char tmp[8], o, k;
	unsigned short bufsize=4096;
	float cont;
	unsigned char n_refresh;

	SDLKey key=0;
	UDPsocket udpsock;
	UDPpacket *packet;

	GABOR gb;
	GRAPHICS gp;
	TRIAL tl;
	memset(&gp, 0, sizeof(GRAPHICS));
	memset(&gb, 0, sizeof(GABOR));
	
	// some initializations
	init_gfx(&gp);
	GFX_SDL_init(gp);
	IO_flush_events();
	SDLNet_Init();
	set_high_priority();//set priority for this application to high
	COLOR_init_mon("monkeylabb");
	GFX_lookups();
	gb.sigma=0;

	// draw start screen
	GFX_clear_buffer_dkl(0.0, 0.0, 0.0);
	print_start(gp);
	diode_signal(gp, 0);
	wait_for_video_refresh(); 
	SDL_GL_SwapBuffers();

	
	// MAIN LOOP
	while(key==0){
		if(key==0) key = IO_ask_key();
		IO_check_escape();	

		openUDP: //OPEN UDP	
		udpsock=SDLNet_UDP_Open(9999);
		packet=SDLNet_AllocPacket(bufsize);
		while(!SDLNet_UDP_Recv(udpsock, packet)) //wait for packet
		{
			if(key==0) key = IO_ask_key();
			IO_check_escape();
			if(key==SDLK_UP) { N_SD=N_SD + 0.5; sprintf(str,"N_sd %.1f",N_SD);}
			if(key==SDLK_DOWN) { N_SD=N_SD - 0.5;sprintf(str,"N_sd %.1f",N_SD);}
			if(key==SDLK_RIGHT) { N_SD=N_SD + 0.1;sprintf(str,"N_sd %.1f",N_SD);}
			if(key==SDLK_LEFT) { N_SD=N_SD - 0.1;sprintf(str,"N_sd %.1f",N_SD);}
			if(key==SDLK_d) { 
				if(Diode_on==0) Diode_on=1;
				else Diode_on=0;
				sprintf(str,"Diode %d",Diode_on);
			}
			
			if(key){
				GFX_text(-950,520,str, gp.white,0); 
				SDL_GL_SwapBuffers();
				SDL_Delay(500);
				key=0;
			}		
		}	
		IO_flush_events();

		//	for(i=0;i<bufsize;i++) {printf("packet->data[%d] = %d\n",i, packet->data[i]);} //print entire packet	

		//LOOP FOR A GIVEN UDP PACKET
		index=0;
		while(index<=bufsize){	
			set_high_priority();
			// GET STIMULUS INFO FROM UDP PACKET
			if(packet->data[index]==ALL_OFF) {
				GFX_clear_buffers(0,0,0);
				gp.structured_bg=0;
				goto draw;
			}
			if(packet->data[index]==SET_BACK_LUM) {
				index=index+1;
				gp.bgcolor[0]=packet->data[index];
				index=index+1;
				gp.bgcolor[1]=packet->data[index];
				index=index+1;
				gp.bgcolor[2]=packet->data[index];
				GFX_clear_buffers(gp.bgcolor[0], gp.bgcolor[1], gp.bgcolor[2]);
				gp.structured_bg=0;
				goto draw;
			}	
			if(packet->data[index]==SET_DKL_BG) {
				index=index+1;
				gb.lum_cont=(packet->data[index]-100) * 0.005; 
				index=index+1;
				gb.rg_cont=(packet->data[index]-100) * 0.005;  
				index=index+1;
				gb.by_cont=(packet->data[index]-100) * 0.005; 
				GFX_clear_buffer_dkl(gb.lum_cont, gb.rg_cont, gb.by_cont);
				gp.structured_bg=0;
				goto draw;
			}
			if(packet->data[index]==PRELOAD_STIM) {	//load .bmp image
				index=index+1;
				tl.img_type=packet->data[index];
				index=index+1;
				tl.img_num=packet->data[index];
				index=index+1;
				cont=packet->data[index] * .01;				
				tl.img_cont=cont;
				init_img(&tl, gp); 
				gp.structured_bg=1;
				goto sendUDP;
			}			
			if(packet->data[index]==SET_FP_ON_CLR) {
				index=index+1;
				gp.fpcolor[0]=packet->data[index];//fp r
				index=index+1;
				gp.fpcolor[1]=packet->data[index];//fp g
				index=index+1;
				gp.fpcolor[2]=packet->data[index];//fp b
						
				if(packet->data[index+1]==SET_FP_SIZE) {
					index=index+2;
					gp.size[1]=(packet->data[index] * gp.h_deg2pix) / 20; //fixation size 20ths of a deg
				}
	
				gp.fpx=0; gp.fpy=0; //default 			
				if(packet->data[index+1]==SET_FP_LOCATION) {
					index=index+2;
					for (j=0;j<8;j++){ //8 bytes of info for given stimulus coordinates
						tmp[j]=packet->data[index];
						index=index+1;
					}
					gp.fpx=byte2float(tmp[0], tmp[1], tmp[2], tmp[3]);//fp X
					gp.fpy=byte2float(tmp[4], tmp[5], tmp[6], tmp[7]);//fp Y
					gp.fpx=(gp.fpx * gp.h_deg2pix) / 10;
					gp.fpy=(gp.fpy * gp.v_deg2pix) / 10;
				}
				goto sendUDP;		
			}
		
			if(packet->data[index]==STIM_FROM_FIX_POINT) {
				index=index+1;
				k=packet->data[index]; ///put n objects in k
				for (i=0;i<k;i++) { 
					index=index+1;//obj #i
					snum=packet->data[index];					
					for (j=0;j<8;j++){ //8 bytes of info for given stimulus coordinates
						index=index+1;
						tmp[j]=packet->data[index];
					}
					gp.holdx0[snum]=byte2float(tmp[0], tmp[1], tmp[2], tmp[3]);//target X
					gp.holdy0[snum]=byte2float(tmp[4], tmp[5], tmp[6], tmp[7]);//target Y
					gp.holdx0[snum]=(gp.holdx0[snum] * gp.h_deg2pix) / 10;
					gp.holdy0[snum]=(gp.holdy0[snum] * gp.v_deg2pix) / 10;
					gp.holdx0[snum] = gp.holdx0[snum];
					gp.holdy0[snum] = gp.holdy0[snum];
				}
				
				if(packet->data[index+1]==LOAD_RAMP) {
					index=index+2;
					gp.ramp_id=packet->data[index];//obj #i
					for (j=0;j<4;j++){index=index+1;tmp[j]=packet->data[index];}
					gp.direction = byte2float(tmp[0], tmp[1], tmp[2], tmp[3]);
					index=index+1;
					gp.velocity = packet->data[index]; 
					index=index+1;	
					gp.frames = packet->data[index];
					calculate_ramp(&tl, gp, &gb);
				} else {
					gp.ramp_id = 0;
				}

				if(packet->data[index+1]==DRAW_ANNULUS) {
					index=index+2;
					k=packet->data[index]; //put n objects in k
					for (i=0;i<k;i++) {
						index=index+1;//obj #i
						snum=packet->data[index];
						index=index+1;//stim size
						gp.size[snum]=(packet->data[index] * gp.h_deg2pix) / 10;//put size in t-size variable					
						index=index+1;//inner radius, if 0 then circle filled
						gp.inner_radius[snum]=(packet->data[index] * gp.h_deg2pix) / 10;//inner radius  
						index=index+1;//positive contrast
						gp.stim_type[snum]=DRAW_ANNULUS;
					}
				}
				if(packet->data[index+1]==DRAW_BAR) {
					index=index+2;
					k=packet->data[index]; //put n objects in k
					for (i=0;i<k;i++) {
						index=index+1;//obj #i
						snum=packet->data[index];
						index=index+1;
						gp.wdth[snum]=(packet->data[index] * gp.h_deg2pix) / 10;
						index=index+1;
						gp.hght[snum]=(packet->data[index] * gp.h_deg2pix) / 10;
						index=index+1;
						gp.ori[snum]=packet->data[index];
						gp.stim_type[snum]=DRAW_BAR;
					}
				}

				if(packet->data[index+1]==SET_STIM_COLORS) {
					index=index+2;
					k=packet->data[index]; //put n objects in k				
					for (i=0;i<k;i++) {	
						index=index+1;//obj #i
						snum=packet->data[index];
						index=index+1;//main color for r
						gp.tcolor[snum][0]=packet->data[index];//target r 
						index=index+1;//main color for g
						gp.tcolor[snum][1]=packet->data[index];//target g 
						index=index+1;//main color for b
						gp.tcolor[snum][2]=packet->data[index];//target b 
						index=index+1;//inner color for r
						gp.icolor[snum][0]=packet->data[index];
						index=index+1;//inner color for g
						gp.icolor[snum][1]=packet->data[index];
						index=index+1;//inner color for b
						gp.icolor[snum][2]=packet->data[index];
						gp.color_type[snum]=SET_STIM_COLORS;
					}
					gp.color_type[1]=SET_STIM_COLORS;//for fixation point
				}

				if(packet->data[index+1]==SET_DKL) {
					index=index+2;
					k=packet->data[index]; //put n objects in k				
					for (i=0;i<k;i++) {	
						index=index+1;//obj #i
						snum=packet->data[index];
						index=index+1;						
						gp.lum_cont[snum] = (packet->data[index]-100) * 0.005; 
						index=index+1;
						gp.rg_cont[snum] = (packet->data[index]-100) * 0.005; 
						index=index+1;
						gp.by_cont[snum] = (packet->data[index]-100) * 0.005; 
						gp.color_type[snum]=SET_DKL;
					}
					gp.color_type[1]=SET_DKL;//for fixation point
				}

				if(packet->data[index+1]==LOAD_PATTERN) {
					index=index+2;//# gabors
					k=packet->data[index]; //put n objects in k	
					for (i=0;i<k;i++) {	
						index=index+1;//gabor id
						gb.id = packet->data[index];
						index=index+1;
						gb.cont = packet->data[index] * 0.01;
						index=index+1;	
						gb.sigma = packet->data[index] * 0.1;
						for (j=0;j<4;j++){index=index+1;tmp[j]=packet->data[index];}
						gb.r = byte2float(tmp[0], tmp[1], tmp[2], tmp[3]);
						index=index+1;
						gb.sp_f = packet->data[index] * 0.1; 
						index=index+1;	
						gb.sp_phase = packet->data[index];
						index=index+1;
						gb.te_f = packet->data[index] * 0.1;
						index=index+1;
						gb.vel = packet->data[index]; 
						index=index+1;
						gb.tf_flag = packet->data[index]; 
						for (j=0;j<4;j++){index=index+1;tmp[j]=packet->data[index];}
						gb.direction = byte2float(tmp[0], tmp[1], tmp[2], tmp[3]);
						index=index+1;
						tl.frames[gb.id] = packet->data[index];
						index=index+1;						
						gb.lum_cont = (packet->data[index]-100) * 0.01;
						index=index+1;
						gb.rg_cont = (packet->data[index]-100) * 0.01;
						index=index+1;
						gb.by_cont = (packet->data[index]-100) * 0.01;

						tl.tf_flag[gb.id]=gb.tf_flag + 1;

						calculate_position(&tl, gp, gb);
						get_patches(&tl, gp, &gb);	
					}
				}
				if(gp.structured_bg) get_fix_clear_patch(&tl, gp, &gb);//F-OFF patch for structured background
				snum=0;//initialize snum count after stim settings
				goto sendUDP;
			}		
			
		
			// INITIALIZE DRAWING COMMANDS FROM UDP PACKET
			if(packet->data[index]==DRAW_TIFF_IMAGE & packet->data[index+1]==OBJ_ON) {
				draw_image(-tl.img_w/2., -tl.img_h/2., tl.img_w, tl.img_h, 1);	
				goto draw;
			}

			if(packet->data[index]==SWITCH_FIX_POINT) {
				index=index+1;				
				if(packet->data[index]==OBJ_ON) {
					draw_fix(gp);	
				}
				if(packet->data[index]==OBJ_OFF & !gp.structured_bg) {
					erase_fix(gp);
				}
				if(packet->data[index]==OBJ_OFF & gp.structured_bg) {
					clear_patch(&tl, gp, gb, 1);
				}
				goto draw;
			}
			
			if(packet->data[index]==SWITCH_STIM) {
				index=index+1;//number of objects
				k=packet->data[index]; //put n objects in k
				for(i=0;i<k;i++){
					index=index+1;//object number
					snum=packet->data[index];
					index=index+1;//do to object
					if(packet->data[index]==OBJ_ON) {
						draw_target(&tl,gp,snum);
					}
					if(packet->data[index]==OBJ_OFF & !gp.structured_bg) {
						erase_target(&tl,gp,snum);
					}
					if(packet->data[index]==OBJ_OFF & gp.structured_bg) {
						clear_patch(&tl, gp, gb, snum);
					}
				}
				index=index+1;
				if(packet->data[index]==SWITCH_FIX_POINT & packet->data[index+1]==OBJ_ON) {
					draw_fix(gp);
				}
				if(packet->data[index]==SWITCH_FIX_POINT & packet->data[index+1]==OBJ_OFF & !gp.structured_bg) {
					erase_fix(gp);
				}
				if(packet->data[index]==SWITCH_FIX_POINT & packet->data[index+1]==OBJ_OFF & gp.structured_bg) {
					clear_patch(&tl, gp, gb, 1);
				}
	
				goto draw;
			}

			if(packet->data[index]==START_RAMP) { 
				index=index+1;//number of objects
				k=packet->data[index]; //put n objects in k
				for(i=0;i<k;i++){
					index=index+1;//object number
					snum=packet->data[index];
					index=index+1;//do to object

					// always clear fixation point for pursuit
					if(packet->data[index]==OBJ_OFF & !gp.structured_bg) erase_fix(gp);
					if(packet->data[index]==OBJ_OFF & gp.structured_bg) clear_patch(&tl, gp, gb, 1);
					if(packet->data[index]==OBJ_ON) {			
						wait_for_video_refresh();//sync to vertical retrace before 1st frame
						SDL_GL_SwapBuffers(); 
						for(p=0;p<gp.frames;p++) {								
							if(p==0) diode_signal(gp, 255);	
							pursuit_target(&tl, gp, snum+p);
							wait_for_video_refresh();
							SDL_GL_SwapBuffers(); //draw							
						//	if(p==Nrefreshdelay) SDLNet_UDP_Send(udpsock, packet->channel, packet); //delay on return signal
							if(p<gp.frames-1) {
								if(!gp.structured_bg) pursuit_target_clear(&tl, gp, snum+p);
								if(gp.structured_bg) clear_patch(&tl, gp, gb, snum+p);
							}
							if(SDLNet_UDP_Recv(udpsock, packet)){
								clear_patch(&tl, gp, gb, snum+p);
								wait_for_video_refresh();
								SDL_GL_SwapBuffers(); //draw
								goto closeUDP;							
							}
						}					
						diode_signal(gp, 0);//stim off signal	
						wait_for_video_refresh();
						SDL_GL_SwapBuffers(); //draw
						//SDL_Delay(gp.monitor_delay);//delay on return signal
						goto sendUDP;
					}
				}
			}
			
			if(packet->data[index]==DRAW_USER_PATTERN) {
				index=index+1;//number of objects
				k=packet->data[index]; //put n objects in k				
				for(i=0;i<k;i++){
					index=index+1;//object number
					snum=packet->data[index]; //printf("snum %d\n", snum);
					index=index+1;//do to object
					//if (snum!=index)

					if(packet->data[index]==OBJ_OFF) clear_patch(&tl, gp, gb, snum);						

					if(packet->data[index]==OBJ_ON) {
						if(tl.tf_flag[snum]==1) display_patch(&tl, gp, gb, snum);

						if(tl.tf_flag[snum]>1) {		
							wait_for_video_refresh();//sync to vertical retrace before 1st frame
							SDL_GL_SwapBuffers(); 
							for(p=0;p<tl.frames[snum];p++) {								
								if(p==0) diode_signal(gp, 255);								
								display_patch(&tl, gp, gb, snum+p);								
								wait_for_video_refresh();
								SDL_GL_SwapBuffers(); //draw							
								//if(p==Nrefreshdelay) SDLNet_UDP_Send(udpsock, packet->channel, packet); //delay on return signal
								if(p<tl.frames[snum]-1) clear_patch(&tl, gp, gb, snum+p);								
								if(SDLNet_UDP_Recv(udpsock, packet)){
									clear_patch(&tl, gp, gb, snum+p);
									wait_for_video_refresh();
									SDL_GL_SwapBuffers(); //draw
									goto closeUDP;							
								}
							}
							if(packet->data[index+1]==OBJ_OFF) clear_patch(&tl, gp, gb, snum);								
							diode_signal(gp, 0);//stim off signal	
							wait_for_video_refresh();
							SDL_GL_SwapBuffers(); //draw
						//	SDL_Delay(gp.monitor_delay);//delay on return signal
							goto sendUDP;
						}
					}
				}
				goto draw;
			}
			//same as previous except image is drawn simultaneous with patch
			if(packet->data[index]==DRAW_TIFF_IMAGE & packet->data[index+1]==DRAW_USER_PATTERN) {
				index=index+2;//number of objects
				k=packet->data[index]; //put n objects in k				
				for(i=0;i<k;i++){
					index=index+1;//object number
					snum=packet->data[index]; //printf("snum %d\n", snum);
					index=index+1;//do to object

					draw_image(-tl.img_w/2., -tl.img_h/2., tl.img_w, tl.img_h, 1);

					if(packet->data[index]==OBJ_OFF) clear_patch(&tl, gp, gb, snum);						

					if(packet->data[index]==OBJ_ON) {
						if(tl.tf_flag[snum]==1) display_patch(&tl, gp, gb, snum);

						if(tl.tf_flag[snum]>1) {		
							wait_for_video_refresh();//sync to vertical retrace before 1st frame
							SDL_GL_SwapBuffers(); 
							for(p=0;p<tl.frames[snum];p++) {								
								if(p==0) diode_signal(gp, 255);								
								display_patch(&tl, gp, gb, snum+p);								
								wait_for_video_refresh();
								SDL_GL_SwapBuffers(); //draw							
								//if(p==Nrefreshdelay) SDLNet_UDP_Send(udpsock, packet->channel, packet); //delay on return signal
								if(p<tl.frames[snum]-1) clear_patch(&tl, gp, gb, snum+p);								
								if(SDLNet_UDP_Recv(udpsock, packet)){
									clear_patch(&tl, gp, gb, snum+p);
									wait_for_video_refresh();
									SDL_GL_SwapBuffers(); //draw
									goto closeUDP;							
								}
							}
							if(packet->data[index+1]==OBJ_OFF) clear_patch(&tl, gp, gb, snum);								
							diode_signal(gp, 0);//stim off signal	
							wait_for_video_refresh();
							SDL_GL_SwapBuffers(); //draw
						//	SDL_Delay(gp.monitor_delay);//delay on return signal
							goto sendUDP;
						}
					}
				}
				goto draw;
			}


		index=index+1;
		}

		// SWAP BUFFERS					
		draw: //DRAW	
			//if not ENABLE_REX_FRAME_SYNC
			if(packet->data[index+1] != ENABLE_REX_FRAME_SYNC){	
				if(Diode_on) diode_signal(gp, 255);
				wait_for_video_refresh(); //return from here is syncronized to vertical retrace
				SDL_GL_SwapBuffers();//draw		
			//	SDL_Delay(gp.monitor_delay);//delay on return signal
			//	SDLNet_UDP_Send(udpsock, packet->channel, packet);//return signal for REX				
				if(Diode_on) {
					diode_signal(gp, 0);
					wait_for_video_refresh(); //return from here is syncronized to vertical retrace
					SDL_GL_SwapBuffers();//draw
				}
			}

			//if ENABLE_REX_FRAME_SYNC 
			if(packet->data[index+1] == ENABLE_REX_FRAME_SYNC) {
				index=n;
				index=index+1;
				n_refresh=packet->data[index];//number of frames
				for(i=1;i<=n_refresh;i++){
					if(i==1) {// first frame
						if(Diode_on) diode_signal(gp, 255);//init white pulse on first refresh	
						wait_for_video_refresh(); //return from here is syncronized to vertical retrace
						SDL_GL_SwapBuffers();//draw
					//t=UTIL_time_delta(RESET);
					//	SDL_Delay(gp.monitor_delay);//delay on return signal
					//	SDLNet_UDP_Send(udpsock, packet->channel, packet);//return signal for REX
						if(Diode_on) diode_signal(gp, 0);//init black pulse for 2nd refresh	
					}
					if(i>1) {// subsequent refreshes											
						wait_for_video_refresh();//sync 
						SDL_GL_SwapBuffers();//draw				
					}
				}
				//t=UTIL_time_delta(NO_RESET);
				//printf("t = %d ms\n",t);
				SDLNet_UDP_Send(udpsock, packet->channel, packet); //return signal for REX
			}		
			goto closeUDP;

		sendUDP:
			while(SDLNet_UDP_Send(udpsock, packet->channel, packet)==NULL);//return signal for REX	
			goto closeUDP;

		closeUDP: //CLOSE AND FREE UDP	
			for(i=0;i<bufsize;i++) packet->data[i]=0;
			SDLNet_FreePacket(packet);
			SDLNet_UDP_Close(udpsock);
			IO_flush_events();
			goto openUDP;				
	}	
	IO_flush_events();
	leave("regular exit");
	return 0;
}