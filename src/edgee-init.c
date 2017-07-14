#include <R_ext/Rdynload.h>
#include <R.h>
#include "headerEdge.h"

void empEdge(double *dat, int *nc, int *nr, double *d0, double *s20, double *alpha,
	     int *one_smp, char **side, int *ncheck, double *lim, int *unb_mom,
	     double *varpost, double *ts, double *tMs, double *pval, double *pvalM);
void empEdgeOrd(double *dat, int *nc, int *nr, double *alpha, int *one_smp, char **side,
		int *ncheck, double *lim, int *unb_mom, double *ts, double *pval);
void empEdgeOrd0(double *dat, int *nc, int *nr, double *alpha, int *one_smp, char **side,
		int *ncheck, double *lim, int *unb_mom, double *pval);
void empEdgeWelch(double *dat, int *nc, int *nr, double *alpha, char **side,
		  int *ncheck, double *lim, int *unb_mom, double *pval);
void tailDiagR(char **type, double *qargs, double *xleft, double *xright,
	       double *n, int *nch, double *df, int *lthick, int *rthick);

static const R_CMethodDef cMethods[] = {
  {"empEdge",      (DL_FUNC) &empEdge,      16},
  {"empEdgeOrd",   (DL_FUNC) &empEdgeOrd,   11},
  {"empEdgeOrd0",  (DL_FUNC) &empEdgeOrd0,  10},
  {"empEdgeWelch", (DL_FUNC) &empEdgeWelch, 9},
  {"tailDiagR",    (DL_FUNC) &tailDiagR,    9},
  {NULL}
};
  
void R_init_edgee(DllInfo *info)
{
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

