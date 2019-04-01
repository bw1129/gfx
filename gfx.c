#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glut.h>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

#include <EL2/eyelink.h>
#include <EL2/w32_exptsppt2.h>
#include <EL2/w32_demo.h>
//#include <jpeg/jpeglib.h>

#include "exp.h"
#include "maths.h"

static const char *arrow[] = {
  /* width height num_colors chars_per_pixel */
  "    32    32        3            1",
  /* colors */
  ". c #000000",
  "X c #ffffff",
  "  c None",
  /* pixels -9,-5*/ 
 /*0123456789*/         
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "    ...........                 ",
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "         .                      ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "0,0"
};
static const char *invisible[] = {
  /* width height num_colors chars_per_pixel */
  "    32    32        3            1",
  /* colors */
  ". c #000000",
  "X c #ffffff",
  "  c None",
  /* pixels -9,-5*/ 
 /*0123456789*/         
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "                                ",
  "0,0"
};

static SDL_Cursor *init_system_cursor(const char *image[])
{
  int i, row, col;
  Uint8 data[4*32];
  Uint8 mask[4*32];
  int hot_x, hot_y;

  i = -1;
  for ( row=0; row<32; ++row ) {
    for ( col=0; col<32; ++col ) {
      if ( col % 8 ) {
        data[i] <<= 1;
        mask[i] <<= 1;
      } else {
        ++i;
        data[i] = mask[i] = 0;
      }
      switch (image[4+row][col]) {
        case 'X':
          data[i] |= 0x01;
          //k[i] |= 0x01;
          break;
        case '.':
          mask[i] |= 0x01;
          break;
        case ' ':
          break;
      }
    }
  }
  sscanf(image[4+row], "%d,%d", &hot_x, &hot_y);
  return SDL_CreateCursor(data, mask, 32, 32, hot_x, hot_y);
}

void GFX_enter_2Dmode(void)
{
	SDL_Surface *screen;
	screen = SDL_GetVideoSurface();

	glViewport(0, 0, screen->w, screen->h);

	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();

	glOrtho((GLdouble)screen->w/-2., (GLdouble)screen->w/2., (GLdouble)screen->h/-2.0, (GLdouble)screen->h/2.0, 0.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}

void GFX_leave_2Dmode(void)
{
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();

	glPopAttrib();
}

int GFX_SDL_init(GRAPHICS gp)
{
	int value, bpp=16;
	float gamma=0.0;
	Uint32 video_flags=0;
	SDL_Cursor *cursor;

	if( SDL_Init( SDL_INIT_VIDEO ) < 0 ) {
		fprintf(stderr,"Couldn't initialize SDL: %s\n",SDL_GetError());
		exit( 1 );
	}

	SDL_GL_SetAttribute( SDL_GL_RED_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_GREEN_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_BLUE_SIZE, 5 );
	SDL_GL_SetAttribute( SDL_GL_DEPTH_SIZE, bpp );
	SDL_GL_SetAttribute( SDL_GL_DOUBLEBUFFER, 1 );
	
	if ( SDL_SetVideoMode( gp.h_pixel, gp.v_pixel,  bpp, gp.video_flags ) == NULL ) {
		fprintf(stderr, "Couldn't set GL mode: %s\n", SDL_GetError());
		SDL_Quit();
		exit(1);
	}

	SDL_GL_GetAttribute( SDL_GL_RED_SIZE, &value );
	SDL_GL_GetAttribute( SDL_GL_GREEN_SIZE, &value );
	SDL_GL_GetAttribute( SDL_GL_BLUE_SIZE, &value );
	SDL_GL_GetAttribute( SDL_GL_DEPTH_SIZE, &value );
	SDL_GL_GetAttribute( SDL_GL_DOUBLEBUFFER, &value );


	SDL_WM_SetCaption( "SDL GL test", "GRAPHICS" );
	if ( gamma != 0.0 ) SDL_SetGamma(gamma, gamma, gamma);
	
	SDL_WM_GrabInput(SDL_GRAB_ON);
	cursor = init_system_cursor(arrow);
	SDL_SetCursor(cursor);
	SDL_ShowCursor(SDL_DISABLE);
    atexit(SDL_Quit);
	GFX_enter_2Dmode();
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	return(0);
}

void GFX_check_for_graphics_errors(void)		/* Check for error conditions. */
{
	GLenum gl_error;
	char* sdl_error;
	gl_error = glGetError( );
	if( gl_error != GL_NO_ERROR ) {
		fprintf( stderr, "testgl: OpenGL error: %d\n", gl_error );
	}
	sdl_error = SDL_GetError();
	if( sdl_error[0] != '\0' ) {
		fprintf(stderr, "testgl: SDL error '%s'\n", sdl_error);
		SDL_ClearError();
	}

}
void GFX_rect_filled(int x, int y, int x_w, int y_w, int ori, unsigned char *c)
{
	float newx,newy;
	float r, th;

	r=sqrt((x*x)+(y*y));
	th=atan2(y,x);
	th=th+(ori*DEG2RAD);
	if(th<0) th=th+2*PI;
	newx=r*cos(th);
	newy=r*sin(th);
	newx=x-newx;
	newy=y-newy;

	glPushMatrix();
	glTranslatef(newx, newy, 0.0);
	glRotatef(ori, 0.0, 0.0, 1.0);
	glColor3f(c[0]/255.,c[1]/255.,c[2]/255.);
	glRectf(x-x_w/2, y-y_w/2, x+x_w/2, y+y_w/2);
	glPopMatrix();
}

void GFX_rect_filled_DKL(int x, int y, int x_w, int y_w, int ori, double lum, double rg, double by)
{
	double c[3];
	float newx,newy;
	float r, th;
	
	r=sqrt((x*x)+(y*y));
	th=atan2(y,x);
	th=th+(ori*DEG2RAD);
	if(th<0) th=th+2*PI;
	//	printf("th %1f\n", th);
	newx=r*cos(th);
	newy=r*sin(th);
	newx=x-newx;
	newy=y-newy;
	
	COLOR_colRGB(lum, rg, by, c);
	glPushMatrix();	
	glTranslatef(newx, newy, 0.0);
	glRotatef(ori, 0.0, 0.0, 1.0);
	glColor3f(c[0],c[1],c[2]);
	glRectf(x-x_w/2, y-y_w/2, x+x_w/2, y+y_w/2);
	glPopMatrix();
}

void GFX_rect(int x, int y, int w, int h, unsigned char *c)
{
	glPushMatrix();
	glBegin(GL_LINES); {
		glColor3f(c[0]/255.,c[1]/255.,c[2]/255.);
		glVertex2f(x-w/2., y-h/2.);			
		glVertex2f(x-w/2., y+h/2.);
		glVertex2f(x-w/2.-1, y-h/2.);			
		glVertex2f(x-w/2.-1, y+h/2.);
		glVertex2f(x-w/2., y+h/2.);
		glVertex2f(x+w/2., y+h/2.);
		glVertex2f(x+w/2., y+h/2.);
		glVertex2f(x+w/2., y-h/2.);
		glVertex2f(x+w/2.+1, y+h/2.);
		glVertex2f(x+w/2.+1, y-h/2.);
		glVertex2f(x+w/2., y-h/2.);
		glVertex2f(x-w/2., y-h/2.);			
	} glEnd();
	glPopMatrix();
}
void GFX_text(float x, float y, char *s, unsigned char *c, int center) 
{
    int lines, length=0, max=0;
    char* p;
    for(p = s; *p; p++) {
		if (*p == '\n') {
			if(length>max) max=length;
			length=0;
		}
		length += glutBitmapWidth(GLUT_BITMAP_HELVETICA_18, *p);
    }
	if(max) length = max;
    glPushMatrix();
    glColor3f(c[0]/255.,c[1]/255.,c[2]/255.);
	if(center) x=x-length/2.;
	glRasterPos2i(x, y);
    for(p = s, lines = 0; *p; p++) {
		if (*p == '\n') {
			lines++;
			glRasterPos2i(x, y-(lines*18));
		}
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *p);
    }
    glPopMatrix();
}
void GFX_triangle(int x, int y, int w, int dir, unsigned char *c)
{
	int cp=w*6;
	glPushMatrix();
	glBegin(GL_POLYGON); {
		glColor3f(c[0]/255.,c[1]/255.,c[2]/255.);
		if(dir) {
			glVertex2f(x-w, y-w);			
			glVertex2f(x,y+w);
			glVertex2f(x+w, y-w);
		} else {
			glVertex2f(x+w, y+w);			
			glVertex2f(x,y-w);
			glVertex2f(x-w, y+w);
		}
	} glEnd();
	glPopMatrix();
}
//unsigned char *GFX_readimage_jpg(char *fn, int *nr, int *nc, float cont, int gray_flag)
//{
//	FILE * file = fopen (fn, "rb");
//	unsigned long l;
//	unsigned char *rp, *im, gray;
//	double rgb;
//	int cnt=0,min=1000,max=0,j=0;
//	printf("noise contrast=%.2f\n",cont);
//	if (file != NULL) {
//		struct jpeg_decompress_struct cinfo;
//		struct jpeg_error_mgr jerr;
//		JSAMPARRAY buffer;
//		int row_stride;
//		cinfo.err = jpeg_std_error(&jerr);
//		jpeg_create_decompress(&cinfo);
//		jpeg_stdio_src(&cinfo, file);
//		jpeg_read_header(&cinfo, TRUE);
//		jpeg_start_decompress(&cinfo);
//
//		row_stride = cinfo.output_width * cinfo.output_components;
//		buffer = (*cinfo.mem->alloc_sarray) ((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);
//		
//		*nr = cinfo.output_height; //1024
//		*nc = cinfo.output_width;  //1024
//		rp = im = calloc(cinfo.output_height*cinfo.output_width*3, 1);
//		//fprintf(stderr,"row_stride %d\n",row_stride);
//		//fprintf(stderr,"memory load: %d\n",cinfo.output_height*cinfo.output_width*3);
//		
//		if (im != 0) {
//			/* read */
//			while (cinfo.output_scanline < cinfo.output_height) {
//				/*fprintf(stderr, "cinfo.output_components=%d, cinfo.output_width=%d, cinfo.output_height=%d while cinfo.output_scanline=%d\n",
//					cinfo.output_components,
//					cinfo.output_width,
//					cinfo.output_height,
//					cinfo.output_scanline);*/
//				jpeg_read_scanlines(&cinfo, buffer, 1);
//				l = 0;
//				rp = im + (cinfo.output_height - cinfo.output_scanline) * cinfo.output_width*3;
//				cnt = 0;
//				while (l < row_stride) {
//					//fprintf(stderr,"l=%d while row_stride=%d\n", l,row_stride);					
//					rgb = ((float)(buffer[0][l])/255.-0.5) * cont + 0.5;
//					if(gray_flag) gray = rgb*255+0.5;
//					else gray = 128;// no image read, just neutral gray 
//					//gray = buffer[0][l];
//					//gray = cnt%255;						
//			/* if grayscale image it runs 3 times with same value, if rgb, runs once for each value */
//					for(j=4-cinfo.output_components; j>=1 ;j--) {
//						*rp++ = gray;
//					}					
//					l++;
//					gray = buffer[0][l];
//					if(gray < min) min = gray;
//					if(gray > max) max = gray;
//					cnt++;
//				}
//			}
//		}
//		jpeg_finish_decompress(&cinfo);
//		jpeg_destroy_decompress(&cinfo);
//		fclose(file);
//		fprintf(stderr,"min=%d max=%d\n",min,max);
//		return im;
//	} else {
//		fprintf (stderr, "File %s not found\n", fn);
//		return 0;
//	}
//}



unsigned char *GFX_readimage_bmp(char *filename, int *nr, int *nc, float cont, int rgbflag)
{
	// reads bitmap into SDL surface and returns rgb values.
	unsigned char *image, *bytepointer, *r, *g, *b;	
	Uint32 pixelvalue;
	SDL_Color* theKey;
	double c[3];
	Uint32 c2[3];
	int x,y,i;

	SDL_Surface* image_original; // create SDL surface to draw image on	
	SDL_Rect srect;
	image_original=SDL_LoadBMP(filename); //load image, put it in this variable
	*nr=image_original->h;
	*nc=image_original->w;
	image = (unsigned char*) calloc((image_original->h)*(image_original->w)*3, 1);

//	converts SDL surface to standard RGB pixel values
	for (y=0; y<image_original->h; y++){
		for (x=0; x<image_original->w; x++){	
			if(rgbflag) {
				theKey=getPixel(image_original,x,y);//unsigned char	

				theKey->r=((float)(theKey->r/255.-0.5) * cont + 0.5) * 255 + 0.5;//convert to float to adjust image contrast
				theKey->g=((float)(theKey->g/255.-0.5) * cont + 0.5) * 255 + 0.5;//convert to float to adjust image contrast
				theKey->b=((float)(theKey->b/255.-0.5) * cont + 0.5) * 255 + 0.5;//convert to float to adjust image contrast

				r=(unsigned char*) &theKey->r;//convert back to unsigned char
				g=(unsigned char*) &theKey->g;//convert back to unsigned char
				b=(unsigned char*) &theKey->b;//convert back to unsigned char

				image[x*3+(image_original->h - y-1)*image_original->w*3]=*r;
				image[x*3+(image_original->h - y-1)*image_original->w*3+1]=*g;
				image[x*3+(image_original->h - y-1)*image_original->w*3+2]=*b;
			}

			if(rgbflag==0) {
				pixelvalue=getPixel(image_original,x,y);//unsigned char	
				pixelvalue=((float)(pixelvalue/255.-0.5) * cont + 0.5) * 255 + 0.5;//convert to float to adjust image contrast
				bytepointer=(unsigned char*) &pixelvalue;//convert back to unsigned char

				image[x*3+(image_original->h - y-1)*image_original->w*3]=*bytepointer;
				image[x*3+(image_original->h - y-1)*image_original->w*3+1]=*bytepointer;
				image[x*3+(image_original->h - y-1)*image_original->w*3+2]=*bytepointer;
			}
					
		}
	}
	return image;
}


Uint32 getPixel(SDL_Surface *screen, int x, int y)
{		
    int bpp = screen->format->BytesPerPixel;	
	SDL_Color theKey;
	Uint32 pixel=0;    
    /* Here p is the address to the pixel we want to set */
    Uint8 *p = (Uint8 *)screen->pixels + y * screen->pitch + x * bpp;

    switch(bpp) {
    case 1:
        pixel = *p;
		return pixel;
        break;

    case 2:
        pixel = *(Uint16 *)p;
		return pixel;
        break;

    case 3:
	    //	fprintf(stderr,"ERROR: gfx.c->getPixel(...) failed because 24bit color mode is unsupported yet (use 32 instead)!\n");
	    // no three-byte-compatibility yet
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN) 
			pixel = (p[0] << 16 | p[1] << 8 | p[2]);
		else
			pixel = (p[0] | p[1] << 8 | p[2] << 16);
		
		SDL_GetRGB(pixel, screen->format, &theKey.r, &theKey.g, &theKey.b);
		return &theKey;
		break;

    case 4:
        pixel = *(Uint32 *)p;
		return pixel;
        break;
	}
}
void putPixel(SDL_Surface *screen, int x, int y, Uint32 pixel)
{
	int bpp = screen->format->BytesPerPixel;    
    /* Here p is the address to the pixel we want to set */
    Uint8 *p = (Uint8 *)screen->pixels + y * screen->pitch + x * bpp;

    switch(bpp) {
    case 1:
        *p = pixel;
        break;

    case 2:
        *(Uint16 *)p = pixel;
        break;

    case 3:
        if(SDL_BYTEORDER == SDL_BIG_ENDIAN) {
            p[0] = (pixel >> 16) & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = pixel & 0xff;
        } else {
            p[0] = pixel & 0xff;
            p[1] = (pixel >> 8) & 0xff;
            p[2] = (pixel >> 16) & 0xff;
        }
        break;

    case 4:
        *(Uint32 *)p = pixel;
        break;
    }
}

void GFX_draw_target(GRAPHICS gp, float x, float y)
{
	GFX_circle_filled(x, y, gp.bull_w, gp.black);
}

void GFX_circle_filled(int x, int y, int w, unsigned char *c)
{
	int i, cp=w*6;
	float angle;
	glPushMatrix();
	glBegin(GL_POLYGON); {
		glColor3f(c[0]/255., c[1]/255., c[2]/255.);
		for(i = 0; i < cp; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
	} glEnd();
	glPopMatrix();
}


void GFX_circle_filled_DKL(int x, int y, int w, double lum, double rg, double by)
{
	int i, cp=w*6;
	float angle;
	double c[3];

	COLOR_colRGB(lum, rg, by, c);

	glPushMatrix();
	glBegin(GL_POLYGON); {
		glColor3f(c[0], c[1], c[2]);
		for(i = 0; i < cp; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
	} glEnd();
	glPopMatrix();
}

void GFX_circle_line(int x, int y, int w, unsigned char *c)
{
	int i, cp=w*6;
	float angle;
	glPushMatrix();
	glBegin(GL_LINES); {
		glColor3f(c[0]/255., c[1]/255., c[2]/255.);
		for(i = 0; i < cp; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
		for(i = 1; i < cp+1; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
	} glEnd();
	glPopMatrix(); 
}

void GFX_circle_line_DKL(int x, int y, int w, double lum, double rg, double by)
{
	int i, cp=w*6;
	float angle;
	double c[3];

	COLOR_colRGB(lum, rg, by, c);

	glPushMatrix();
	glBegin(GL_LINES); {
		glColor3f(c[0], c[1], c[2]);
		for(i = 0; i < cp; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
		for(i = 1; i < cp+1; i++) {
			angle = 2*PI*i/(float)(cp);
			glVertex2f(w * (float)cos(angle) + x, w * (float)sin(angle) + y);			
		}
	} glEnd();
	glPopMatrix(); 
}

void GFX_line(int x1, int y1, int x2, int y2, unsigned char *c)
{
	glPushMatrix();
	glBegin(GL_LINES); {
		glColor3f(c[0]/255., c[1]/255., c[2]/255.);
		glVertex2f(x1, y1);			
		glVertex2f(x2, y2);
	} glEnd();
	glPopMatrix();
}
void GFX_clear_buffer_dkl(double lum, double rg, double by)
{
	double b[3];
	COLOR_colRGB(lum, rg, by, b);
	glClearColor(b[0],b[1],b[2],0);
//	glClearColor(.5,.5,.5,0);
	glClear( GL_COLOR_BUFFER_BIT);
}
void GFX_clear_buffers(int r, int g, int b)
{
//	printf("r %d, g %d, b %d\n",r,g,b);
	glClearColor(r/255., g/255., b/255.,0);
	glClear( GL_COLOR_BUFFER_BIT);
}
int GFX_refresh(void)
{
	int i,t,f=10;
	char buf[200], r[3]={255,0,0};
	SDL_Delay(200);
	UTIL_time_delta(RESET);
	for(i=0;i<f;i++) 
		SDL_GL_SwapBuffers();
	t = UTIL_time_delta(NO_RESET);
	sprintf(buf,"refresh %.2f",1000/(t/(float)f));
	GFX_text(50,50,buf, r,0);
	SDL_GL_SwapBuffers();
	SDL_Delay(200);
	return(1);
}

