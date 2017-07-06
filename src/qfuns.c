#include <math.h>

double He0(double x);
double He1(double x);
double He2(double x);
double He3(double x);
double He4(double x);
double He5(double x);
double He6(double x);
double He7(double x);
double He8(double x);
double He9(double x);
double He10(double x);
double He11(double x);

double q1k(double x, double *args) // args: k's, k21 = r^2
{
  double k12, k31, r2, r;
  k12 = args[0];
  r2  = args[2];  /*k21*/
  k31 = args[5];
  r = sqrt(r2);
  return -1./6*He2(x)*k31/(r2*r) - He0(x)*k12/r;
}
double q2k(double x, double *args)
{
  double k12, k22, k31, k41, r2;
  k12 = args[0];
  r2  = args[2];  /*k21*/  
  k22 = args[3];
  k31 = args[5];
  k41 = args[7];
  return -1./72*He5(x)*pow(k31, 2)/pow(r2, 3) - 1./24*(4*k12*k31 + k41)*He3(x)/pow(r2, 2) - 1./2*(pow(k12, 2) + k22)*He1(x)/r2;
}
double q3k(double x, double *args)
{
  double k12, k13, k22, k31, k32, k41, k51, r2, r;
  k12 = args[0];
  k13 = args[1];
  r2  = args[2];  /*k21*/
  k22 = args[3];
  k31 = args[5];
  k32 = args[6];
  k41 = args[7];
  k51 = args[9];
  r = sqrt(r2);
  return -1./1296*He8(x)*pow(k31, 3)/pow(r, 9) - 1./144*(2*k12*pow(k31, 2) + k31*k41)*He6(x)/pow(r, 7) - 1./120*(10*pow(k12, 2)*k31 + 10*k22*k31 + 5*k12*k41 + k51)*He4(x)/pow(r, 5) - 1./6*(pow(k12, 3) + 3*k12*k22 + k32)*He2(x)/(r2*r) - He0(x)*k13/r;
}
double q4k(double x, double *args) 
{
  double k12, k13, k22, k23, k31, k32, k41, k42, k51, k61, r2;
  k12 = args[0];
  k13 = args[1];
  r2  = args[2];  /*k21*/
  k22 = args[3];
  k23 = args[4];
  k31 = args[5];
  k32 = args[6];
  k41 = args[7];
  k42 = args[8];
  k51 = args[9];
  k61 = args[10];
  return -1./31104*He11(x)*pow(k31, 4)/pow(r2, 6) - 1./5184*(4*k12*pow(k31, 3) + 3*pow(k31, 2)*k41)*He9(x)/pow(r2, 5) - 1./5760*(40*pow(k12, 2)*pow(k31, 2) + 40*k22*pow(k31, 2) + 40*k12*k31*k41 + 5*pow(k41, 2) + 8*k31*k51)*He7(x)/pow(r2, 4) - 1./720*(20*pow(k12, 3)*k31 + 60*k12*k22*k31 + 15*pow(k12, 2)*k41 + 20*k31*k32 + 15*k22*k41 + 6*k12*k51 + k61)*He5(x)/pow(r2, 3) - 1./24*(pow(k12, 4) + 6*pow(k12, 2)*k22 + 3*pow(k22, 2) + 4*k13*k31 + 4*k12*k32 + k42)*He3(x)/pow(r2, 2) - 1./2*(2*k12*k13 + k23)*He1(x)/r2;
}

double q1s(double x, double *args)
{
  double lam3;
  lam3 = args[0];
  return 1./6*(2*pow(x, 2) + 1)*lam3;
}
double q2s(double x, double *args)
{
  double lam3, lam4;
  lam3 = args[0];
  lam4 = args[1];
  return -1./18*(pow(x, 5) + 2*pow(x, 3) - 3*x)*pow(lam3, 2) - 1./4*pow(x, 3) + 1./12*(pow(x, 3) - 3*x)*lam4 - 3./4*x;
}
double q3s(double x, double *args)
{
  double lam3, lam4, lam5;
  lam3 = args[0];
  lam4 = args[1];
  lam5 = args[2];
  return 1./1296*(8*pow(x, 8) + 28*pow(x, 6) - 210*pow(x, 4) - 525*pow(x, 2) - 105)*pow(lam3, 3) - 1./144*(4*pow(x, 6) - 30*pow(x, 4) - 90*pow(x, 2) - 15)*lam3*lam4 + 1./24*(2*pow(x, 6) - 3*pow(x, 4) - 6*pow(x, 2))*lam3 - 1./40*(2*pow(x, 4) + 8*pow(x, 2) + 1)*lam5;
}
double q4s(double x, double *args)
{
  double lam3, lam4, lam5, lam6;
  lam3 = args[0];
  lam4 = args[1];
  lam5 = args[2];
  lam6 = args[3];
  return -1./32*pow(x, 7) - 1./1944*(pow(x, 11) + 5*pow(x, 9) - 90*pow(x, 7) - 450*pow(x, 5) + 45*pow(x, 3) + 945*x)*pow(lam3, 4) - 5./96*pow(x, 5) - 1./72*(pow(x, 9) - 6*pow(x, 7) - 12*pow(x, 5) - 18*pow(x, 3) - 9*x)*pow(lam3, 2) - 1./288*(pow(x, 7) - 21*pow(x, 5) + 33*pow(x, 3) + 111*x)*pow(lam4, 2) + 1./60*(pow(x, 7) + 8*pow(x, 5) - 5*pow(x, 3) - 30*x)*lam3*lam5 - 7./96*pow(x, 3) + 1./432*(9*pow(x, 7) - 63*pow(x, 5) + 2*(pow(x, 9) - 12*pow(x, 7) - 90*pow(x, 5) + 36*pow(x, 3) + 261*x)*pow(lam3, 2) + 81*pow(x, 3) + 189*x)*lam4 - 1./90*(2*pow(x, 5) - 5*pow(x, 3) - 15*x)*lam6 - 7./32*x;
}

double He0(double x)
{
  return 1;
}
double He1(double x)
{
  return x;
}
double He2(double x)
{
  return pow(x, 2) - 1;
}
double He3(double x)
{
  return pow(x, 3) - 3*x;
}
double He4(double x)
{
  return pow(x, 4) - 6*pow(x, 2) + 3;
}
double He5(double x)
{
  return pow(x, 5) - 10*pow(x, 3) + 15*x;
}
double He6(double x)
{
  return pow(x, 6) - 15*pow(x, 4) + 45*pow(x, 2) - 15;
}
double He7(double x)
{
  return pow(x, 7) - 21*pow(x, 5) + 105*pow(x, 3) - 105*x;
}
double He8(double x)
{
  return pow(x, 8) - 28*pow(x, 6) + 210*pow(x, 4) - 420*pow(x, 2) + 105;
}
double He9(double x)
{
  return pow(x, 9) - 36*pow(x, 7) + 378*pow(x, 5) - 1260*pow(x, 3) + 945*x;
}
double He10(double x)
{
  return pow(x, 10) - 45*pow(x, 8) + 630*pow(x, 6) - 3150*pow(x, 4) + 4725*pow(x, 2) - 945;
}
double He11(double x)
{
  return pow(x, 11) - 55*pow(x, 9) + 990*pow(x, 7) - 6930*pow(x, 5) + 17325*pow(x, 3) - 10395*x;
}
