#define PI		(float)3.1415926535897
#define DEG2RAD (float) (PI / 180.)
#define RAD2DEG (float) (180. / PI)

static float sqrarg;

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define ITMAX  100
#define EPS    3.0e-7
#define FPMIN  1.0e-30

float round_fl(float number);
float p2(float number);
int round(float number);
double mean(unsigned short *data, short n );
double standard_deviation(unsigned short *data, double mean, short n );
void fit(float *x, float *y, int ndata, float sig[], int mwt, float *a,float *b, float *siga, float *sigb, float *chi2, float *q);
float abs_fl(float number);
float mean_fl(float *data, short n);
