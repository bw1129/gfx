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

#include "exp.h"
#include "maths.h"

static Uint16 Rgamma[256], Ggamma[256], Bgamma[256];

static double dkl2rgb[3][3], rgb2dkl[3][3];
static double rgb2lms[3][3], lms2rgb[3][3];

static double Xmon[3][3];
static double Xmat[3][3];

void COLOR_colRGB(double lum, double cb, double tc, double *c)	/* dkl [-0.5 0.5] -> rgb [0 1] */
{
	c[0] = 0.5 + dkl2rgb[0][0] * lum + dkl2rgb[0][1] * cb + dkl2rgb[0][2] * tc;
	c[1] = 0.5 + dkl2rgb[1][0] * lum + dkl2rgb[1][1] * cb + dkl2rgb[1][2] * tc;
	c[2] = 0.5 + dkl2rgb[2][0] * lum + dkl2rgb[2][1] * cb + dkl2rgb[2][2] * tc;
	//printf("dkl2rgb: %.4f %.4f %.4f\n",dkl2rgb[0][0],dkl2rgb[1][0], dkl2rgb[2][0]);
}

rgbCOL(r, g, b, lum, cb, tc)	/* rgb [0 1] -> dkl [-0.5 0.5] */
double *lum, *cb, *tc, r, g, b;
{
	r -= 0.5; g -= 0.5; b -= 0.5;
	*lum = rgb2dkl[0][0] * r + rgb2dkl[0][1] * g + rgb2dkl[0][2] * b;
	*cb =  rgb2dkl[1][0] * r + rgb2dkl[1][1] * g + rgb2dkl[1][2] * b;
	*tc =  rgb2dkl[2][0] * r + rgb2dkl[2][1] * g + rgb2dkl[2][2] * b;
}

xyYCOL(x, y, Y, r, g, b)	/* xyY -> rgb [0 1] */
double x, y, Y, *r, *g, *b;
{
	double X, Z, z;

	z = 1.0 - x - y; if (z < 0.0) z = 0.0;
	X = (x/y) * Y;
	Z = (z/y) * Y;

	*r = X * Xmat[0][0] + Y * Xmat[0][1] + Z * Xmat[0][2];
	*g = X * Xmat[1][0] + Y * Xmat[1][1] + Z * Xmat[1][2];
	*b = X * Xmat[2][0] + Y * Xmat[2][1] + Z * Xmat[2][2];
}

igammacorrect(r, g, b, R, G, B)
int *R, *G, *B;
{
	*R = Rgamma[r];
	*G = Ggamma[g];
	*B = Bgamma[b];
}

double pow();
#define GAMMA 2.4

gammacorrect(r, g, b)
double *r, *g, *b;
{
	*r = (*r < 0.0)? 0.0: pow(*r, GAMMA);
	*g = (*g < 0.0)? 0.0: pow(*g, GAMMA);
	*b = (*b < 0.0)? 0.0: pow(*b, GAMMA);
}

checkcol(r, g, b)
double *r, *g, *b;
{
	int over = 0;

	if (*r<0) { *r = 0; over++; }
	if (*g<0) { *g = 0; over++; }
	if (*b<0) { *b = 0; over++; }
	if (*r>1) { *r = 1; over++; }
	if (*g>1) { *g = 1; over++; }
	if (*b>1) { *b = 1; over++; }

	return over;
}

scalecol(r, g, b, ri, gi, bi)
double r, g, b;
int *ri, *gi, *bi;
{
	*ri = (int) (r * 255.999);
	*gi = (int) (g * 255.999);
	*bi = (int) (b * 255.999);
}

static double moncie[4][3];
static double monCIE[4][3];
static double monRGB[4][3];
static double monrgb[4][3];
static double white[3];
static double det, wr, wg, wb;

double matinv();

void COLOR_init_mon(m)
char *m;
{
	FILE *fp;
	int i, j, k;
	double sum;
	char buf[256];
	char *getenv();

	if ( ((m == 0) || (*m == 0)) && 
             ( ((m = getenv("MONITOR")) == 0) ||
	       ((m = getenv("HOST")) == 0) ) )
		m = "tomato";

	if (*m != '/')
#if defined(sgi) || defined (linux)
		sprintf(buf, "/homes/karl/gsp/lut/%s", m);
#else
		sprintf(buf, "c:\\src\\lut\\%s", m);
#endif
	else
		strcpy(buf, m);

	fprintf(stdout, "loading settings from %s.\n", buf);

	loadlut(buf, "r", Rgamma, 256);
	loadlut(buf, "g", Ggamma, 256);
	loadlut(buf, "b", Bgamma, 256);

	for(i=0;i<256;i++) {
		//printf("%d. %d %d %d >> ",i,Rgamma[i],Ggamma[i],Bgamma[i]);
		Rgamma[i] = Rgamma[i]<<8;
		Ggamma[i] = Ggamma[i]<<8;
		Bgamma[i] = Bgamma[i]<<8;
		//printf("%d %d %d\n",Rgamma[i],Ggamma[i],Bgamma[i]);
	}
	SDL_SetGammaRamp(Rgamma, Ggamma, Bgamma);

	strcat(buf, ".xyY");
	if ((fp = fopen(buf, "r")) == 0) {
		fprintf(stderr, "cannot open %s.\n", buf);
		exit(1);
	}

	for (i=0; i<3; i++)
		if (fscanf(fp, "%lg %lg %lg", &moncie[i][0], &moncie[i][1], &monCIE[i][1]) != 3) {
			fprintf(stderr, "error reading %s.\n", buf);
			exit(1);
		}

	fclose(fp);
	getdkl(moncie[0][0], moncie[0][1], monCIE[0][1],
	       moncie[1][0], moncie[1][1], monCIE[1][1],
	       moncie[2][0], moncie[2][1], monCIE[2][1], dkl2rgb);
	matinv(dkl2rgb, rgb2dkl);

	for (i=0; i<3; i++) {
		moncie[i][2] = 1.0 - moncie[i][0] - moncie[i][1];
		monCIE[i][0] = (moncie[i][0]/moncie[i][1])*monCIE[i][1];
		monCIE[i][2] = (moncie[i][2]/moncie[i][1])*monCIE[i][1];
		monRGB[i][0] = 0.15514 * monCIE[i][0] + 
				0.54312 * monCIE[i][1] - 0.03286 * monCIE[i][2];
		monRGB[i][1] = -0.15514 * monCIE[i][0] + 
				0.45684 * monCIE[i][1] + 0.03286 * monCIE[i][2];
		monRGB[i][2] = 0.01608 * monCIE[i][2];
		sum = monRGB[i][0] + monRGB[i][1];
		monrgb[i][0] = monRGB[i][0] / sum;
		monrgb[i][1] = monRGB[i][1] / sum;
		monrgb[i][2] = monRGB[i][2] / sum;
	}
	
        for (i=0; i<3; i++)
                for (j=0; j<3; j++)
                        Xmon[i][j] = monCIE[j][i];
	det = matinv(Xmon, Xmat);
	for (i=0; i<3; i++) for (j=0; j<3; j++)
		rgb2lms[i][j] = monRGB[i][j];
	matinv(rgb2lms, lms2rgb);
	mon2cones(0.5, 0.5, 0.5, &wr, &wg, &wb);
}	

mon2cones(r, g, b, R, G, B)
double r, g, b, *R, *G, *B;
{
	double c[3];
	int i;

	*R = r * monRGB[0][0] + g * monRGB[1][0] + b * monRGB[2][0];
	*G = r * monRGB[0][1] + g * monRGB[1][1] + b * monRGB[2][1];
	*B = r * monRGB[0][2] + g * monRGB[1][2] + b * monRGB[2][2];
}

cones2mon(R, G, B, r, g, b)
double R, G, B, *r, *g, *b;
{
	*r = R * lms2rgb[0][0] + G * lms2rgb[1][0] + B * lms2rgb[2][0];
	*g = R * lms2rgb[0][1] + G * lms2rgb[1][1] + B * lms2rgb[2][1];
	*b = R * lms2rgb[0][2] + G * lms2rgb[1][2] + B * lms2rgb[2][2];
fprintf(stderr, "cones2mon: %g %g %g (white %g %g %g) -> %g %g %g\n",
	R, G, B, wr, wg, wb, *r, *g, *b);
}

coco2mon(DR, DG, DB, r, g, b)
double DR, DG, DB, *r, *g, *b;
{
	cones2mon((DR*wr)+wr, (DG*wg)+wg, (DB*wb)+wb, r, g, b);
}

getgrays(r, g, b)
double *r, *g, *b;
{
        double lum;

        lum = monCIE[0][1] + monCIE[1][1] + monCIE[2][1];
        *r = monCIE[0][1]/lum;
        *g = monCIE[1][1]/lum;
        *b = monCIE[2][1]/lum;
}


static loadlut(rfn, ext, addr, len)
char *rfn, *ext;
short *addr;
{
	char buf[256];
	FILE *fp;
	int i;
	long v;
	
	sprintf(buf, "%s.%s", rfn, ext);

	if ((fp = fopen(buf, "r")) == NULL) {
		fprintf(stderr, "cannot open %s.\n", buf);
		exit(1);
	}
	
	for (i=0; i<len; i++)
		if (fscanf(fp, "%hd", &addr[i]) != 1) {
			fprintf(stderr, "error reading %s.\n", buf);
			exit(1);
		}

	fclose(fp);
}

static loadlong(rfn, ext, addr, len)
char *rfn, *ext;
long *addr;
{
	char buf[256];
	FILE *fp;
	int i;
	long v;
	
	sprintf(buf, "%s.%s", rfn, ext);

	if ((fp = fopen(buf, "r")) == NULL) {
		fprintf(stderr, "cannot open %s.\n", buf);
		exit(1);
	}
	
	for (i=0; i<len; i++)
		if (fscanf(fp, "%ld", &addr[i]) != 1) {
			fprintf(stderr, "error reading %s.\n", buf);
			exit(1);
		}

	fclose(fp);
}

#define N 3
#define abs(x) ((x) > 0.0 ? (x) : 0.0)
double matinv(b, a)
double a[N][N], b[N][N];
    {
    register i,j,k;
    double det, biga, hold;
    int l[N], m[N];		/* Row and column permutation vectors */

    for (i=0; i<N; i++)
    	for (j=0; j<N; j++)
		a[i][j] = b[i][j];

    det = 1.0;
    for (k=0;  k<N;  k++)
	{
	l[k] = k;  m[k] = k;
	biga = a[k][k];

	/* Find the biggest element in the submatrix */
	for (i=k;  i<N;  i++)
	    for (j=k;  j<N;  j++)
		if (abs(a[i][j]) > abs(biga))
		    {
		    biga = a[i][j];
		    l[k] = i;
		    m[k] = j;
		    }

	/* Interchange rows */
	i = l[k];
	if (i>k)
	    for (j=0;  j<N;  j++)
		{
	        hold = -a[k][j];
		a[k][j] = a[i][j];
		a[i][j] = hold;
		}

	/* Interchange columns */
	j = m[k];
	if (j>k)
	    for (i=0;  i<N;  i++)
		{
		hold = -a[i][k];
		a[i][k] = a[i][j];
		a[i][j] = hold;
		}

	/* Divide column by minus pivot
	    (value of pivot element is contained in biga). */
	if (biga == 0.0) return(0.0);
	for (i=0;  i<N;  i++)
	    if (i != k) a[i][k] /= -biga;

	/* Reduce matrix */
	for (i=0;  i<N;  i++)
	    if (i != k)
		{
		hold = a[i][k];
		for (j=0;  j<N;  j++)
		    if (j != k)
			a[i][j] += hold * a[k][j];
		}

	/* Divide row by pivot */
	for (j=0;  j<N;  j++)
	    if (j != k) a[k][j] /= biga;
	
	det *= biga;	/* Product of pivots */
	a[k][k] = 1.0/biga;
	}	/* K loop */

    /* Final row & column interchanges */
    for (k=N-1;  k>=0;  k--)
	{
	i = l[k];
	if (i>k)
	    for (j=0;  j<N;  j++)
		{
		hold = a[j][k];
		a[j][k] = -a[j][i];
		a[j][i] = hold;
		}
	j = m[k];
	if (j>k)
	    for (i=0;  i<N;  i++)
		{
		hold = a[k][i];
		a[k][i] = -a[j][i];
		a[j][i] = hold;
		}
	}
    return det;
    }


#define TWOPI (3.14159265*2.0)

static cart_polar(x, y, azi, amp)
double x, y, *azi, *amp;
{
	*amp = sqrt(x*x+y*y);
	*azi = atan2(y, x) / TWOPI * 360.0;
}

polar_cart(azi, ampl, x, y)
double azi, ampl, *x, *y;
{
	*y = sin(azi/360.0*TWOPI) * ampl;
	*x = cos(azi/360.0*TWOPI) * ampl;
}

cart_polar3d(lum, cb, tc, az, el, amp)
double lum, cb, tc, *az, *el, *amp;
{
	double isolen;

	*amp = sqrt(lum*lum + cb*cb + tc*tc);
	*az = atan2(tc, cb) / TWOPI * 360.0;
	isolen = sqrt(cb*cb + tc*tc);
	*el = atan2(lum, isolen) / TWOPI * 360.0;
}

polar_cart3d(az, el, amp, lum, cb, tc)
double az, el, amp, *lum, *cb, *tc;
{
        double costh = cos(el*TWOPI/360.0),
            sinphi = sin(az*TWOPI/360.0),
            cosphi = cos(az*TWOPI/360.0);

        if (amp == 0) {
                *lum = *cb = *tc = 0.0;
        } else {
                *lum = sin(el*TWOPI/360.0) * amp;
                *cb = costh * cosphi * amp;
                *tc = costh * sinphi * amp;
        }
}

static getdkl(rx, ry, rY, gx, gy, gY, bx, by, bY, dkl2rgb)
double rx, ry, rY, gx, gy, gY, bx, by, bY, dkl2rgb[3][3];
{
	double solvex(), solvey();
	double white[4];
	double delta[4];
	double r[3], g[3], b[3];
	double x, y;
	int i;
	double dRtc, dGcb, dGtc, dBcb;

	dox(&r[0], &b[0], rx, ry, rY); white[0] = rY/2.0; g[0] = 1.0 - r[0];
	dox(&r[1], &b[1], gx, gy, gY); white[1] = gY/2.0; g[1] = 1.0 - r[1];
	dox(&r[2], &b[2], bx, by, bY); white[2] = bY/2.0; g[2] = 1.0 - r[2];

	/* CONSTANT BLUE AXIS */
	dGcb = delta[1] = solvex(white[0]*b[0], white[0]*(r[0]+g[0]),
                          b[1],b[2],r[1]+g[1],r[2]+g[2]);
	dBcb = delta[2] = solvey(white[0]*b[0], white[0]*(r[0]+g[0]),
                          b[1],b[2],r[1]+g[1],r[2]+g[2]);

	/* TRITAN CONFUSION AXIS */
	dRtc = delta[0] = solvex(white[2]*r[2],white[2]*g[2],r[0],r[1],g[0],g[1]);
	dGtc = delta[1] = solvey(white[2]*r[2],white[2]*g[2],r[0],r[1],g[0],g[1]);

#define IMAX 1.0

	dkl2rgb[0][0] = IMAX;
	dkl2rgb[0][1] = IMAX;
	dkl2rgb[0][2] = dRtc*IMAX/white[0];
	dkl2rgb[1][0] = IMAX;
	dkl2rgb[1][1] = -dGcb*IMAX/white[1];
	dkl2rgb[1][2] = dGtc*IMAX/white[1];
	dkl2rgb[2][0] = IMAX;
	dkl2rgb[2][1] = -dBcb*IMAX/white[2];
	dkl2rgb[2][2] = -IMAX;
}

static double solvex(a,b,c,d,e,f)
double a,b,c,d,e,f;
{
	double x;
	x = (a*f/d - b)/(c*f/d - e);
	return(x);
}

static double solvey(a,b,c,d,e,f)
double a,b,c,d,e,f;
{
	double y;
	y = (a*e/c -b)/(d*e/c - f);
	return(y);
}

static dox(r, b, x, y)
double *r, *b, x, y;
{
	double z = 1.0 - x - y;
	double g, sum;

	if (y == 0.0)
		return 0;

	*r = 0.15514 * x + 0.54312 * y - 0.03286 * z;
	g = -0.15514 * x + 0.45684 * y + 0.03286 * z;
	*b = 0.01608 * z;

	sum = *r + g;
	*r /= sum;
	*b /= sum;

	if (*r > 1.0 || *b > 1.0)
		return 0;
	
	return 1;
}
