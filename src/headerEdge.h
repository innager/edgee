#define EPSILON sqrt(DBL_EPSILON)
#define STREQ(a, b) (strcmp((a), (b)) == 0)

void getMuBias(double *smp, double n, double *mu);
void getMuUnb( double *smp, double n, double *mu);
void getm(     double *smp, double n, double *mu);
void getMuBias2(double *smp, int n_x, int n_y, double *mu);
void getMuUnb2( double *smp, int n_x, int n_y, double *mu);
void getm2(     double *smp, int n_x, int n_y, double *mu);
void getLambda(double *mu, double *lam);

void calculateK1smp(double *mu, double A, double B, double *k);
void calculateK2smp(double *mu_x, double *mu_y, double A, double B_x, double B_y,
		     double b_x, double b_y, double *k);
void tailDiag(char *qfun, double *args, double *xleft, double *xright,
	      double *yl0, double *yr0, double *npow, int nch, double *df,
	      int ldiag, int rdiag, int *lthick, int *rthick);

/* Edgeworth coefficients */
double q1k(double x, double *args);
double q2k(double x, double *args);
double q3k(double x, double *args);
double q4k(double x, double *args);
double q1s(double x, double *args);
double q2s(double x, double *args);
double q3s(double x, double *args);
double q4s(double x, double *args);

double M3two(double m3, int n_x, int n_y);
double M4two(double m2, double m4, int n_x, int n_y);
double M5two(double m2, double m3, double m5, int n_x, int n_y);
double M6two(double m2, double m3, double m4, double m6, int n_x, int n_y);
