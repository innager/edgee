#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "headerEdge.h"

/*-------------------------------------------------------------------------------*/
/* returns p-values for 0 - 4 term expansions, ordinary t-statistic only; 
   one- and two-sample design, *one_smp has to be either 1 or 0; 
   for two-sample: dat is a vector c(Treat, Control); nc = c(ntreat, ncontrol); 
   switch X and Y because Edgeworth expansions are set up for Xbar - Ybar, 
                                      while regression is for Ybar - Xbar        */
/*-------------------------------------------------------------------------------*/
 
void empEdgeOrd(double *dat, int *nc, int *nr, double *alpha, int *one_smp, char **side,
		int *ncheck, double *lim, int *unb_mom, double *ts, double *pval)
{
  int i, j, nx, ny, m, ldiag, rdiag;
  int nch = *ncheck + 1;
  int lthick[4], rthick[4];
  double n, Cxy, A, bx, by, Bx, By, r, t, incr, xrtmp;
  double mu[5], df[5], p[5], k[11], npow[5], xleft[nch], xright[nch], yl0[nch], yr0[nch];
  double smp[*nc + *(nc + 1)];  // for one_smp nc[1] = 0
  double (*q1)(double, double *), (*q2)(double, double *), (*q3)(double, double *),
    (*q4)(double, double *);
  void (*getMu)(double *, double, double *),
       (*getMu2)(double *, double, double, double *);  
  
  m = *nr;
  if (*one_smp) {
    n = (double) *nc;  
    df[0] = n - 1;                   /* degrees of freedom */    
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
    r = sqrt(df[0]/(nx + ny));
  }
  
  df[1] = df[0] + 2;
  df[2] = df[1] + 3;
  df[3] = df[2] + 3;
  df[4] = df[3] + 3;
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
    p[0]  = pt(t/r,   df[0],  1, 0);

    /* if alpha < p[0] < 1 - alpha, don't use Edgeworth */
    if (*alpha < p[0] && p[0] < 1 - *alpha) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = fmin(p[0],  1 - p[0]); /* fill out pval */
      }
      continue;                                 /* proceed with the next row of data */
    }

    /* if one-sided test and t on the other side, don't need Edgeworth */
                           /* 0.5 < p-value < 1 */
    if (STREQ(side[0], "left") && t > 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = p[0];                 
      }
      continue;
    }
    if (STREQ(side[0], "right") && t < 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = 1 - p[0];
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
      tailDiag("qs", k,  xleft, xright, yl0,  yr0,  npow, nch, df,  ldiag, rdiag,
	       lthick,  rthick);
    } else {                                  /* TWO-SAMPLE t */
      for (j = 0; j < nx + ny; j++) {
	smp[j] = dat[j*m + i];                 
      }                                       /* a row of dat */
      getMu2(smp, nx, ny, mu);
      A = Cxy*(bx + by)*mu[0];               
      calculateK2smp(mu, mu, A,  Bx,  By,  bx, by, k);
      tailDiag("qk", k,  xleft, xright, yl0,  yr0,  npow, nch, df,  ldiag, rdiag,
	       lthick,  rthick);
    }

    /* probabilities based on tail diagnostic */
    if ((t < 0 && lthick[0]) || (t > 0 && rthick[0])) {      /* thick, term 1 */
      p[1] =  p[0]  + npow[1] * q1(t/r,    k)  * dt(t/r,   df[1],  0);
    } else {                                                /* thin */
      p[1]  = p[0];
    }
    if ((t < 0 && lthick[1]) || (t > 0 && rthick[1])) {      /* thick, term 2 */
      p[2] =  p[1]  + npow[2] * q2(t/r,    k)  * dt(t/r,   df[2], 0);
    } else {
      p[2]  = p[1];
    }
    if ((t < 0 && lthick[2]) || (t > 0 && rthick[2])) {      /* thick, term 3 */
      p[3] =  p[2]  + npow[3] * q3(t/r,    k)  * dt(t/r,   df[3], 0);
    } else {
      p[3]  = p[2];
    }
    if ((t < 0 && lthick[3]) || (t > 0 && rthick[3])) {       /* thick, term 4 */
      p[4] =  p[3]  + npow[4] * q4(t/r,    k)  * dt(t/r,   df[4], 0);
    } else {
      p[4]  = p[3];
    }

    if (t < 0) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = p[j];                   /* fill out pval */
      }
    } else {
      for (j = 0; j < 5; j++) {
	pval[j*m + i]  = 1 - p[j];  
      }
    }
  }
}

