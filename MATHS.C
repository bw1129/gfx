#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>

#include "exp.h" 
#include "maths.h"

float round_fl(float number)
{
   if((number - floor(number)) > .5)
   return ((float)(ceil(number)));
   else
   return ((float)(floor(number)));
}
float abs_fl(float number)
{
   if(number<0) return -1*number;
   else return number;
}
int round(float number)
{
   if((number - floor(number)) > .5)
   return ((int)(ceil(number)));
   else
   return ((int)(floor(number)));
}
static int power_of_two(int input)
{
	int value = 1;

	while ( value < input ) {
		value <<= 1;
	}
	return value;
}
float p2(float number)
{
   return(number * number);
}

double mean(unsigned short *data, short n )
{
   short i;
   double mean = 0.0;

   for(i=0;i<n;i++)
      mean += (double)data[i];
   mean /= ((double)n);
   return mean;
}

float mean_fl(float *data, short n )
{
   short i;
   float mean = 0.0;

   for(i=0;i<n;i++)
      mean += data[i];
   mean /= n;
   return mean;
}
double standard_deviation(unsigned short *data, double mean, short n )
{
  short i;
  double val = 0.0;

  for(i=0;i<n;i++)
      val += ( ((double)data[i] - mean) * ((double)data[i] - mean) );
  val /= ((double)n);
  val = sqrt( val );
  return val;
}
float gammln(float xx)
{
   double x,tmp,ser;
   static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
   int j;
   x=xx-1.0;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.0;
   for (j=0;j<=5;j++) {
       x += 1.0;
       ser += cof[j]/x;
   }
   return -tmp+log(2.50662827465*ser);
}
void gcf(float *gammcf, float a, float x, float *gln)
{
   int n;
   float gold=0.0,g,fac=1.0,b1=1.0;
   float b0=0.0,anf,ana,an,a1,a0=1.0;
   *gln=gammln(a);
   a1=x;
   for (n=1;n<=ITMAX;n++) {
       an=(float) n;
       ana=an-a;
       a0=(a1+a0*ana)*fac;
       b0=(b1+b0*ana)*fac;
       anf=an*fac;
       a1=x*a0+anf*a1;
       b1=x*b0+anf*b1;
       if (a1) {
           fac=1.0/a1;
	   g=b1*fac;
	   if (fabs((g-gold)/g) < EPS) {
	       *gammcf=exp(-x+a*log(x)-(*gln))*g;
	       return;
           }
	   gold=g;
       }
   }
   fprintf(stderr,"a too large, ITMAX too small in routine GCF\n");
   printf("\a\n");
}
void gser(float *gamser, float a, float x, float *gln)
{
   int n;
   double sum,del,ap;
   *gln=gammln(a);
   if (x <= 0.0) {
       if (x < 0.0) { fprintf(stderr,"x less than 0 in routine GSER\n"); printf("\a\n");}
       *gamser=0.0;
       return;
   } else {
       ap=a;
       del=sum=1.0/a;
       for (n=1;n<=ITMAX;n++) {
           ap += 1.0;
	   del *= x/ap;
	   sum += del;
	   if (fabs(del) < fabs(sum)*EPS) {
	       *gamser=sum*exp(-x+a*log(x)-(*gln));
	       return;
	   }
       }
       printf("a too large, ITMAX too small in routine GSER\n");
       printf("\a\n");
       return;
   }
}
float gammq(float a, float x)
{
   float gamser,gammcf,gln;
   if (x < 0.0 || a <= 0.0) { fprintf(stderr,"Invalid arguments in routine GAMMQ\n"); printf("\a\n"); }
   if (x < (a+1.0)) {
       gser(&gamser,a,x,&gln);
       return 1.0-gamser;
   } else {
       gcf(&gammcf,a,x,&gln);
       return gammcf;
   }
}
void fit(float *x, float *y, int ndata, float sig[], int mwt, float *a,float *b, float *siga, float *sigb, float *chi2, float *q)
{
   int i;
   float wt,t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss,sigdat;
   *b=0.0;
   if (mwt) {
       ss=0.0;
       for (i=0;i<ndata;i++) {
           wt=1.0/SQR(sig[i]);
	   ss += wt;
	   sx += x[i]*wt;
           sy += y[i]*wt;
       }
   } else {
       for (i=0;i<ndata;i++) {
           sx += x[i];
	   sy += y[i];
       }
       ss=ndata;
   }
   sxoss=sx/ss;
   if (mwt) {
       for (i=0;i<ndata;i++) {
	   t=(x[i]-sxoss)/sig[i];
	   st2 += t*t;
	   *b += t*y[i]/sig[i];
       }
   } else {
       for (i=0;i<ndata;i++) {
	   t=x[i]-sxoss;
	   st2 += t*t;
	   *b += t*y[i];
       }
   }
   *b /= st2;
   *a=(sy-sx*(*b))/ss;
   *siga=sqrt((1.0+sx*sx/(ss*st2))/ss);
   *sigb=sqrt(1.0/st2);
   *chi2=0.0;
   if (mwt == 0) {
       for (i=0;i<ndata;i++)
	   *chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
       *q=1.0;
       sigdat=sqrt((*chi2)/(ndata-2));
       *siga *= sigdat;
       *sigb *= sigdat;
   } else {
       for (i=0;i<ndata;i++)
	   *chi2 += SQR((y[i]-(*a)-(*b)*x[i])/sig[i]);
       *q=gammq(0.5*(ndata-2),0.5*(*chi2));
   }
}

int find_index(char x[], int max, int min, char value)
{
	int index;
	for(index=min; index < max; index++)
	{
		if(x[index] == value)
		{
			return (index);
		}
	}
	return (-1);
}


