#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "headerEdge.h"

/*-------------------------------------------------------------------------------*/
/* returns p-values for 0 - 4 term expansions for Welch t-test;
   for two-sample: dat is a vector c(Treat, X), nc = c(ntreat, ncontrol); 
   switch X and Y because Edgeworth expansions are set up for Xbar - Ybar, 
                                      while regression is for Ybar - Xbar        */
/*-------------------------------------------------------------------------------*/

void empEdgeWelch(double *dat, int *nc, int *nr, double *alpha, char **side,
		int *ncheck, double *lim, int *unb_mom, double *pval)
{
  int i, j, nx, ny, m, ldiag, rdiag;
  int nch = *ncheck + 1;
  int lthick[4], rthick[4];
  double n, Cx, Cy, A, bx, by, Bx, By, r, t, incr, xrtmp;
  double mux[5], muy[5], df[5], p[5], k[11], npow[5], xleft[nch], xright[nch],
    yl0[nch], yr0[nch], tAr[3];
  double smpx[*nc], smpy[*(nc + 1)]; 
  void (*getMu)(double *, double, double *);
  
  m = *nr;
  nx = *nc;
  ny = *(nc + 1);
  n = (nx + ny)/2.;
  df[0] = (double) (nx + ny - 2);  /* degrees of freedom */
  Cx = (double) nx/(nx - 1);
  Cy = (double) ny/(ny - 1);
  bx  = n/nx;
  by  = n/ny;
  Bx  = Cx*bx;
  By  = Cy*by;
  
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

  if (*unb_mom) {
    getMu  = getMuUnb;
  } else {
    getMu  = getMuBias;
  }
  if (STREQ(side[0], "two-sided")) {
    *alpha /= 2;
  }

  /*-------------------------------------------*/
  /* the main loop going through features/rows */
  /*-------------------------------------------*/
  for (i = 0; i < m; i++) {
    for (j = 0; j < nx; j++) {
      smpx[j] = dat[j*m + i];                  /* X (treatment) */
    }
    for (j = nx; j < nx + ny; j++) {
      smpy[j - nx] = dat[j*m + i];             /* Y (control)   */
    }
    get_tArWelch(smpx, smpy, nx, ny, bx, by, Bx, By, tAr);
    t = tAr[0];
    r = tAr[2];
    p[0] = pt(t/r, df[0], 1, 0);

    /* if alpha < p[0] < 1 - alpha, don't use Edgeworth */
    if (*alpha < p[0] && p[0] < 1 - *alpha) {
      for (j = 0; j < 5; j++) {
	pval[j*m + i] = fmin(p[0], 1 - p[0]);  /* fill out pval */
      }
      continue;                                /* proceed with the next row of data */
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

    getMu(smpx, nx, mux);
    getMu(smpy, ny, muy);
    A = tAr[1];
    calculateK2smp(mux, muy, A, Bx, By, bx, by, k);
    tailDiag("qk", k,  xleft, xright, yl0,  yr0,  npow, nch, df,  ldiag, rdiag,
	     lthick,  rthick);
 
     /* probabilities based on tail diagnostic */
    if ((t < 0 && lthick[0]) || (t > 0 && rthick[0])) {      /* thick, term 1 */
      p[1] =  p[0]  + npow[1] * q1k(t/r,    k)  * dt(t/r,   df[1],  0);
    } else {                                                /* thin */
      p[1]  = p[0];
    }
    if ((t < 0 && lthick[1]) || (t > 0 && rthick[1])) {      /* thick, term 2 */
      p[2] =  p[1]  + npow[2] * q2k(t/r,    k)  * dt(t/r,   df[2], 0);
    } else {
      p[2]  = p[1];
    }
    if ((t < 0 && lthick[2]) || (t > 0 && rthick[2])) {      /* thick, term 3 */
      p[3] =  p[2]  + npow[3] * q3k(t/r,    k)  * dt(t/r,   df[3], 0);
    } else {
      p[3]  = p[2];
    }
    if ((t < 0 && lthick[3]) || (t > 0 && rthick[3])) {       /* thick, term 4 */
      p[4] =  p[3]  + npow[4] * q4k(t/r,    k)  * dt(t/r,   df[4], 0);
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

void get_tArWelch(double *smpx, double *smpy, double nx, double ny,
	          double bx, double by, double Bx, double By, double *tAr) {
  int i;
  double mx1, mx2, my1, my2;

  mx1 = 0;
  mx2 = 0;
  my1 = 0;
  my2 = 0;
  for (i = 0; i < nx; i++) {
    mx1 += smpx[i];
    mx2 += pow(smpx[i], 2);
  }
  for (i = 0; i < ny; i++) {
    my1 += smpy[i];
    my2 += pow(smpy[i], 2);
  }
  mx1 /= nx;
  mx2 = (mx2 - nx*pow(mx1, 2))/(nx - 1);
  my1 /= ny;
  my2 = (my2 - ny*pow(my1, 2))/(ny - 1);
  tAr[0] = (mx1 - my1)/sqrt(mx2/nx + my2/ny);  /* t-stat */
  tAr[1] = Bx*mx2 + By*my2;                    /* A */
  tAr[2] = sqrt((bx*mx2 + by*my2)/tAr[1]);     /* r */
}
    
  
  

