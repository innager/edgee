#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>
#include "headerEdge.h"

void tailDiag(char *qfun, double *args, double *xleft, double *xright,
	      double *yl0, double *yr0, double *npow, int nch, double *df,
	      int ldiag, int rdiag, int *lthick, int *rthick) 
/* xleft and xright going away from center, so xleft is decreasing */
{
  int i, j;
  double yl1old, yr1old, yl1new, yr1new, yl2old, yr2old, yl2new, yr2new,
    yl3old, yr3old, yl3new, yr3new, yl4old, yr4old, yl4new, yr4new;
  double (*q1)(double, double *), (*q2)(double, double *), (*q3)(double, double *),
    (*q4)(double, double *);

  if        (STREQ(qfun, "qk")) {
    q1 = q1k;
    q2 = q2k;
    q3 = q3k;
    q4 = q4k;
  } else if (STREQ(qfun, "qs")) {
    q1 = q1s;
    q2 = q2s;
    q3 = q3s;
    q4 = q4s;
  } 
  for (j = 0; j < 4; j++) {
    lthick[j] = 1;
    rthick[j] = 1;
  }

  /* right tail */  
  if (rdiag) {
    yr1old = yr0[0] + npow[1] * q1(xright[0], args) * dt(xright[0], df[1], 0);
    yr2old = yr1old + npow[2] * q2(xright[0], args) * dt(xright[0], df[2], 0);
    yr3old = yr2old + npow[3] * q3(xright[0], args) * dt(xright[0], df[3], 0);
    yr4old = yr3old + npow[4] * q4(xright[0], args) * dt(xright[0], df[4], 0);    

    for (i = 1; i < nch; i++) { /* start from second element */
      yr1new = yr0[i] + npow[1] * q1(xright[i], args) * dt(xright[i], df[1], 0);
      if (yr1new > 1 + EPSILON      ||    /* not bounded */
	  yr1old - yr1new > EPSILON ||    /* not monotonic */
	  yr1new > yr0[i]) {              /* not conservative */
	for (j = 0; j < 4; j++) rthick[j] = 0;
	break;
      }
      yr1old = yr1new;
      if (rthick[1]) {
	yr2new = yr1old + npow[2] * q2(xright[i], args) * dt(xright[i], df[2], 0);
	if (yr2new > 1 + EPSILON || yr2old - yr2new > EPSILON || yr2new > yr0[i]) {
	  for (j = 1; j < 4; j++) rthick[j] = 0;
	}
	yr2old = yr2new;
      }
      if (rthick[2]) {
	yr3new = yr2old + npow[3] * q3(xright[i], args) * dt(xright[i], df[3], 0);
	if (yr3new > 1 + EPSILON || yr3old - yr3new > EPSILON || yr3new > yr0[i]) {
	  for (j = 2; j < 4; j++) rthick[j] = 0;
	}
	yr3old = yr3new;
      }
      if (rthick[3]) {
	yr4new = yr3old + npow[4] * q4(xright[i], args) * dt(xright[i], df[4], 0);
	if (yr4new > 1 + EPSILON || yr4old - yr4new > EPSILON || yr4new > yr0[i]) {
	  rthick[3] = 0;
	}
	yr4old = yr4new;
      }
    }
  }
    
  /* left tail */
  if (ldiag) {
    yl1old = yl0[0] + npow[1] * q1(xleft[0], args) * dt(xleft[0], df[1], 0);
    yl2old = yl1old + npow[2] * q2(xleft[0], args) * dt(xleft[0], df[2], 0);
    yl3old = yl2old + npow[3] * q3(xleft[0], args) * dt(xleft[0], df[3], 0);
    yl4old = yl3old + npow[4] * q4(xleft[0], args) * dt(xleft[0], df[4], 0);

    for (i = 1; i < nch; i++) {
      yl1new = yl0[i] + npow[1] * q1(xleft[i], args) * dt(xleft[i], df[1], 0);
      if (yl1new < 0 - EPSILON      ||
	  yl1new - yl1old > EPSILON || /* should be decreasing as xleft decreasing */
	  yl1new < yl0[i]) {
	for (j = 0; j < 4; j++) lthick[j] = 0;
	break;
      }
      yl1old = yl1new;
      if (lthick[1]) {
	yl2new = yl1old + npow[2] * q2(xleft[i], args) * dt(xleft[i], df[2], 0);
	if (yl2new < 0 - EPSILON || yl2new - yl2old > EPSILON || yl2new < yl0[i]) {
	  for (j = 1; j < 4; j++) lthick[j] = 0;
	}
	yl2old = yl2new;
      }
      if (lthick[2]) {
	yl3new = yl2old + npow[3] * q3(xleft[i], args) * dt(xleft[i], df[3], 0);
	if (yl3new < 0 - EPSILON || yl3new - yl3old > EPSILON || yl3new < yl0[i]) {
	  for (j = 2; j < 4; j++) lthick[j] = 0;
	}
	yl3old = yl3new;
      }
      if (lthick[3]) {
	yl4new = yl3old + npow[4] * q4(xleft[i], args) * dt(xleft[i], df[4], 0);
	if (yl4new < 0 - EPSILON || yl4new - yl4old > EPSILON || yl4new < yl0[i]) {
	  lthick[3] = 0;
	}
	yl4old = yl4new;
      }
    }
  }
}


/* wrapper for R function */
void tailDiagR(char **type, double *qargs, double *xleft, double *xright,
	       double *n, int *nch, double *df, int *lthick, int *rthick)
{
  int i;
  char qfun[] = "qk";
  double npow[5], yl0[*nch], yr0[*nch];

  for (i = 0; i < *nch; i++) {
    yl0[i] = pt(xleft[i],  df[0], 1, 0);
    yr0[i] = pt(xright[i], df[0], 1, 0);
  }
  if (STREQ(type[0], "short")) {
    qfun[1] = 's';
  }
  for (i = 1; i < 5; i++) { // start with second element to match term number 
    npow[i] = pow(*n, -i/2.);
  }  
  tailDiag(qfun, qargs, xleft, xright, yl0, yr0, npow, *nch, df, 1, 1, lthick, rthick);
}
	       


