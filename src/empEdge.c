#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "headerEdge.h"

/*-------------------------------------------------------------------------------*/
/* returns p-values for 0 - 4 term expansions, ordinary and moderated t-statistic; 
   one- and two-sample design, *one_smp has to be either 1 or 0; 
   if two-sided test, *alpha passed to C function is already divided by 2;
   for two-sample: dat is a vector c(Y, X), cbind() for matrices; nc = c(ny, nx); 
   switch X and Y because Edgeworth expansions are set up for Xbar - Ybar, 
                                      while regression is for Ybar - Xbar        */
/*-------------------------------------------------------------------------------*/

void empEdge(double *dat, int *nc, int *nr, double *d0, double *s20, double *alpha,
	     int *one_smp, char **side, int *ncheck, double *lim, int *unb_mom,
	     double *varpost, double *ts, double *tMs, double *pval, double *pvalM)
{
  // possibly a little overhead for VLA but good to have an option ncheck from R. 
  int i, j, nx, ny, m, ldiag, rdiag;
  int nch = *ncheck + 1;
  int lthick[4], rthick[4], lthickM[4], rthickM[4];
  double n, Cxy, A, AM, BM, bx, by, Bx, By, BxM, ByM, r, rM, t, tM, incr, xrtmp;
  double mu[5], df[5], dfM[5], p[5], pM[5], k[11], kM[11], npow[5], 
    xleft[nch], xright[nch], yl0[nch], yr0[nch], ylM0[nch], yrM0[nch];
  double smp[*nc + *(nc + 1)];  // for one_smp nc[1] = 0
  double (*q1)(double, double *), (*q2)(double, double *), (*q3)(double, double *),
    (*q4)(double, double *);
  void (*getMu)(double *, double, double *),
       (*getMu2)(double *, int, int, double *);  
  
  m = *nr;
  if (*one_smp) {
    n = (double) *nc;  
    df[0] = n - 1;                   /* degrees of freedom */    
    BM = n/(*d0 + df[0]);            /* for moderated t */ 
    r = sqrt(df[0]/n);
  } else {
    nx = *nc;
    ny = *(nc + 1);
    n = (nx + ny)/2.;  
    df[0] = (double) (nx + ny - 2);  /* degrees of freedom */
    Cxy = (nx + ny)/df[0];
    bx  = n/nx;
    by  = n/ny;
    Bx  = Cxy*by;
    By  = Cxy*bx;
    BxM = Bx*df[0]/(*d0 + df[0]);
    ByM = By*df[0]/(*d0 + df[0]);
    r = sqrt(df[0]/(nx + ny));
  }
  
  df[1] = df[0] + 2;
  df[2] = df[1] + 3;
  df[3] = df[2] + 3;
  df[4] = df[3] + 3;
  for (j = 0; j < 5; j++) {
    dfM[j] = df[j] + *d0;
  }
  for (j = 1; j < 5; j++) { // start with second element to match term number 
    npow[j] = pow(n, -j/2.);
  }

  /* vectors x for tail diagnostic, length = ncheck + 1 */
  incr = (lim[1] - lim[0])/(double) *ncheck;  
  xrtmp = lim[0];                             
  for (i = 0; i < nch; i++) {
    xright[i] = xrtmp;
    xleft[i]  = -xrtmp;  /* note both tails go out from center */
    xrtmp += incr;
  }
  /* yl0, yr0, yl0M, yr0M for tailDiag() */
  for (i = 0; i < nch; i++) {
    yl0[i]  = pt(xleft[i],  df[0],  1, 0);    /* use same x for all */
    yr0[i]  = pt(xright[i], df[0],  1, 0);
    ylM0[i] = pt(xleft[i],  dfM[0], 1, 0);  
    yrM0[i] = pt(xright[i], dfM[0], 1, 0);
  }

   /* for ordinary t - qs for one-smp, qk for two-smp; for moderated - qk for both */
  if (*one_smp) {                             
    q1 = q1s;
    q2 = q2s;
    q3 = q3s;
    q4 = q4s;
  } else {
    q1 = q1k;
    q2 = q2k;
    q3 = q3k;
    q4 = q4k;
  }

  if (*unb_mom) {
    getMu  = getMuUnb;
    getMu2 = getMuUnb2;
  } else {
    getMu  = getMuBias;
    getMu2 = getMuBias2;
  }
  if (STREQ(side[0], "two-sided")) {
    *alpha /= 2;
  }

  /*-------------------------------------------*/
  /* the main loop going through features/rows */
  /*-------------------------------------------*/
  for (i = 0; i < m; i++) {
    t  = ts[i];
    tM = tMs[i];
    /* rM for either one- or two-sample, relies on *one_smp being 1 or 0 */
    rM = sqrt(dfM[0]*varpost[i]/((*d0)*(*s20) + n*(2 - *one_smp)*varpost[i]));

    p[0]  = pt(t/r,   df[0],  1, 0);
    pM[0] = pt(tM/rM, dfM[0], 1, 0);

    /* if alpha < p[0] < 1 - alpha, don't use Edgeworth */
    if (*alpha < pM[0] && pM[0] < 1 - *alpha) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = fmin(p[0],  1 - p[0]);  /* fill out pval */
	pvalM[j*m + i] = fmin(pM[0], 1 - pM[0]);
      }
      continue;                                 /* proceed with the next row of data */
    }

    /* if one-sided test and t on the other side, don't need Edgeworth */
                           /* 0.5 < p-value < 1 */
    if (STREQ(side[0], "left") && t > 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = p[0];                 
	pvalM[j*m + i] = pM[0];
      }
      continue;
    }
    if (STREQ(side[0], "right") && t < 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = 1 - p[0];
	pvalM[j*m + i] = 1 - pM[0];
      }
      continue;                                
    }

    /*----------------------*/
    /* Edgeworth expansions */
    /*----------------------*/

    if (t < 0) {
      ldiag = 1;                               /* left  tail diagnostic only */
      rdiag = 0;
    } else {        
      ldiag = 0;                               /* right tail diagnostic only */
      rdiag = 1;
    }

    if (*one_smp) {                            /* ONE-SAMPLE t */
      for (j = 0; j < n; j++) {
	smp[j] = dat[j*m + i];                 
      }                                        /* a row of dat */
      getMu(smp, n, mu);
      getLambda(mu, k);                        /* k will hold lambda's */
      mu[0] = varpost[i];                      /* posterior var in k for moder. */
      AM = ((*d0)*(*s20) + n*mu[0])/dfM[0];    /* match numerator in k */
      calculateK1smp(mu, AM, BM, kM);
      tailDiag("qs", k,  xleft, xright, yl0,  yr0,  npow, nch, df,  ldiag, rdiag,
	       lthick,  rthick);
      tailDiag("qk", kM, xleft, xright, ylM0, yrM0, npow, nch, dfM, ldiag, rdiag,
	       lthickM, rthickM);
    } else {                                  /* TWO-SAMPLE t */
      for (j = 0; j < nx + ny; j++) {
	smp[j] = dat[j*m + i];                 
      }                                       /* a row of dat */
      getMu2(smp, nx, ny, mu);
      A = Cxy*(bx + by)*mu[0];               
      calculateK2smp(mu, mu, A,  Bx,  By,  bx, by, k);
      mu[0] = varpost[i];                     /* posterior var */
      AM = (bx + by)*((*d0)*(*s20) + Cxy*df[0]*mu[0])/dfM[0]; 
      calculateK2smp(mu, mu, AM, BxM, ByM, bx, by, kM);
      tailDiag("qk", k,  xleft, xright, yl0,  yr0,  npow, nch, df,  ldiag, rdiag,
	       lthick,  rthick);
      tailDiag("qk", kM, xleft, xright, ylM0, yrM0, npow, nch, dfM, ldiag, rdiag,
	       lthickM, rthickM);
    }

    /* probabilities based on tail diagnostic */
    if ((t < 0 && lthick[0]) || (t > 0 && rthick[0])) {      /* thick, term 1 */
      p[1] =  p[0]  + npow[1] * q1(t/r,    k)  * dt(t/r,   df[1],  0);
    } else {                                                /* thin */
      p[1]  = p[0];
    }
    if ((tM < 0 && lthickM[0]) || (tM > 0 && rthickM[0])) {  /* thickM */
      pM[1] = pM[0] + npow[1] * q1k(tM/rM, kM) * dt(tM/rM, dfM[1], 0);
    } else {                                               
      pM[1] = pM[0];
    }
    if ((t < 0 && lthick[1]) || (t > 0 && rthick[1])) {      /* thick, term 2 */
      p[2] =  p[1]  + npow[2] * q2(t/r,    k)  * dt(t/r,   df[2], 0);
    } else {
      p[2]  = p[1];
    }
    if ((tM < 0 && lthickM[1]) || (tM > 0 && rthickM[1])) {  /* thickM */
      pM[2] = pM[1] + npow[2] * q2k(tM/rM, kM) * dt(tM/rM, dfM[2], 0);
    } else {
      pM[2] = pM[1];
    }
    if ((t < 0 && lthick[2]) || (t > 0 && rthick[2])) {      /* thick, term 3 */
      p[3] =  p[2]  + npow[3] * q3(t/r,    k)  * dt(t/r,   df[3], 0);
    } else {
      p[3]  = p[2];
    }
    if ((tM < 0 && lthickM[2]) || (tM > 0 && rthickM[2])) {   /* thickM */
      pM[3] = pM[2] + npow[3] * q3k(tM/rM, kM) * dt(tM/rM, dfM[3], 0);
    } else {
      pM[3] = pM[2];
    }
    if ((t < 0 && lthick[3]) || (t > 0 && rthick[3])) {       /* thick, term 4 */
      p[4] =  p[3]  + npow[4] * q4(t/r,    k)  * dt(t/r,   df[4], 0);
    } else {
      p[4]  = p[3];
    }
    if ((tM < 0 && lthickM[3]) || (tM > 0 && rthickM[3])) {   /* thickM */
      pM[4] = pM[3] + npow[4] * q4k(tM/rM, kM) * dt(tM/rM, dfM[4], 0);
    } else {
      pM[4] = pM[3];
    }

    if (t < 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = p[j];                   /* fill out pval */
	pvalM[j*m + i] = pM[j];
      }
    } else {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = 1 - p[j];  
	pvalM[j*m + i] = 1 - pM[j];
      }
    }
  }
}


/* Central moment estimates - naive biased */
void getm(double *smp, double n, double *mu)
{
  int i, j;
  double m1 = 0;

  for (i = 0; i < n; i++) {
    m1 += smp[i];
  }
  m1 /= n;
  for (j = 0; j < 5; j++) {              /* central moments */
    mu[j] = 0.;
    for (i = 0; i < n; i++) {
      mu[j] += pow(smp[i] - m1, j + 2);
    }
    mu[j] /= n;
  }
}

/*
void getm2(double *smp, int n_x, int n_y, double *mu)
{
  int i, j;
  double smpx[n_x], smpy[n_y], mx[5], my[5];

  for (i = 0; i < n_x; i++) {
    smpx[i] = smp[i];
  }
  for (i = 0; i < n_y; i++) {
    smpy[i] = smp[n_x + i];
  }
  getm(smpx, n_x, mx);
  getm(smpy, n_y, my);

  for (j = 0; j < 5; j++) {
    mu[j] = (n_x*mx[j] + n_y*my[j])/(n_x + n_y);
  }
  } */
    
void getm2(double *smp, int n_x, int n_y, double *mu)
{
  int i, j;
  double sumx, sumy, mx1, my1;
  double smpx[n_x], smpy[n_y];
  
  mx1 = 0;
  my1 = 0;
  for (i = 0; i < n_x; i++) {
    smpx[i] = smp[i];
    mx1 += smpx[i];
  }
  for (i = 0; i < n_y; i++) {
    smpy[i] = smp[n_x + i];
    my1 += smpy[i];
  }
  mx1 /= n_x;
  my1 /= n_y;

  for (j = 0; j < 5; j++) {
    sumx = 0;
    sumy = 0;
    for (i = 0; i < n_x; i++) {
      sumx += pow(smpx[i] - mx1, j + 2);
    }
    for (i = 0; i < n_y; i++) {
      sumy += pow(smpy[i] - my1, j + 2);
    }
    mu[j] = (sumx + sumy)/(n_x + n_y);
  }
}

/* only second moment unbiased */
void getMuBias(double *smp, double n, double *mu)
{
  getm(smp, n, mu);
  mu[0] *= n/(n - 1);
}

/* only second moment unbiased - pooled variance */
void getMuBias2(double *smp, int n_x, int n_y, double *mu)
{
  getm2(smp, n_x, n_y, mu);
  mu[0] *= (n_x + n_y)/(n_x + n_y - 2);
}
  
void getMuUnb(double *smp, double n, double *mu)
{
  double m2, m3, m4, m5, m6;

  getm(smp, n, mu);
  m2 = mu[0];
  m3 = mu[1];
  m4 = mu[2];
  m5 = mu[3];
  m6 = mu[4];

  mu[0] = m2*n/(n - 1);
  mu[1] = m3*pow(n, 2)/((n - 1)*(n - 2));
  mu[2] = -3*pow(m2, 2)*(2*n - 3)*n/((n - 1)*(n - 2)*(n - 3)) + (pow(n, 2) - 2*n + 3)*m4*n/((n - 1)*(n - 2)*(n - 3));
  mu[3] = -10*m2*m3*pow(n, 2)/((n - 1)*(n - 3)*(n - 4)) + (pow(n, 2) - 5*n + 10)*m5*pow(n, 2)/((n - 1)*(n - 2)*(n - 3)*(n - 4));
  mu[4] = 15*pow(m2, 3)*(3*n - 10)*pow(n, 2)/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) - 40*(pow(n, 2) - 6*n + 10)*pow(m3, 2)*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) - 15*(pow(n, 3) - 8*pow(n, 2) + 29*n - 40)*m2*m4*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5)) + (pow(n, 4) - 9*pow(n, 3) + 31*pow(n, 2) - 39*n + 40)*m6*n/((n - 1)*(n - 2)*(n - 3)*(n - 4)*(n - 5));
}

void getMuUnb2(double *smp, int n_x, int n_y, double *mu)
{
  double m2, m3, m4, m5, m6;

  getm2(smp, n_x, n_y, mu);
  m2 = mu[0];
  m3 = mu[1];
  m4 = mu[2];
  m5 = mu[3];
  m6 = mu[4];

  mu[0] = m2*(n_x + n_y)/(n_x + n_y - 2);
  mu[1] = M3two(m3, n_x, n_y);
  mu[2] = M4two(m2, m4, n_x, n_y);
  mu[3] = M5two(m2, m3, m5, n_x, n_y);
  mu[4] = M6two(m2, m3, m4, m6, n_x, n_y);
}
  
void getLambda(double *mu, double *lam) 
{
  lam[0] = mu[1]/pow(mu[0], 1.5);                           /* lambda_3 */
  lam[1] = mu[2]/pow(mu[0], 2) - 3;                         /* lambda_4 */
  lam[2] = mu[3]/pow(mu[0], 2.5) - 10*lam[0];               /* lambda_5 */
  lam[3] = (mu[4] - 15*mu[2]*mu[0] - 10*pow(mu[1], 2))/pow(mu[0], 3) + 30;
}


