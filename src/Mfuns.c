#include <math.h>

double M3two(double m3, int n_x, int n_y)
{
  return m3*n_x*n_y*(n_x + n_y)/(pow(n_x, 2)*n_y + n_x*pow(n_y, 2) - 6*n_x*n_y + 2*n_x + 2*n_y);
}

double M4two(double m2, double m4, int n_x, int n_y)
{
  return n_x*n_y*(3*pow(m2, 2)*(pow(n_x, 2) + 2*n_x*n_y + pow(n_y, 2))*(4*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + 3*pow(n_x, 2) - 5*n_x*pow(n_y, 2) + 3*pow(n_y, 2)) - m4*n_x*n_y*(n_x + n_y)*(pow(n_x, 3)*n_y + 2*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + n_x*pow(n_y, 3) - 5*n_x*pow(n_y, 2) + 12*n_x*n_y - 3*n_x - 3*n_y))/(3*(pow(n_x, 2)*n_y + n_x*pow(n_y, 2) - 4*n_x*n_y + n_x + n_y)*(4*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + 3*pow(n_x, 2) - 5*n_x*pow(n_y, 2) + 3*pow(n_y, 2)) - (pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 6*pow(n_x, 2)*n_y - 3*pow(n_x, 2) + 6*n_x*pow(n_y, 2) - 3*pow(n_y, 2))*(pow(n_x, 3)*n_y + 2*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + n_x*pow(n_y, 3) - 5*n_x*pow(n_y, 2) + 12*n_x*n_y - 3*n_x - 3*n_y));
}

double M2pow2two(double m2, double m4, int n_x, int n_y)
{
  return n_x*n_y*(-pow(m2, 2)*(pow(n_x, 2) + 2*n_x*n_y + pow(n_y, 2))*(pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 6*pow(n_x, 2)*n_y - 3*pow(n_x, 2) + 6*n_x*pow(n_y, 2) - 3*pow(n_y, 2)) + m4*n_x*n_y*(n_x + n_y)*(pow(n_x, 2)*n_y + n_x*pow(n_y, 2) - 4*n_x*n_y + n_x + n_y))/(3*(pow(n_x, 2)*n_y + n_x*pow(n_y, 2) - 4*n_x*n_y + n_x + n_y)*(4*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + 3*pow(n_x, 2) - 5*n_x*pow(n_y, 2) + 3*pow(n_y, 2)) - (pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 6*pow(n_x, 2)*n_y - 3*pow(n_x, 2) + 6*n_x*pow(n_y, 2) - 3*pow(n_y, 2))*(pow(n_x, 3)*n_y + 2*pow(n_x, 2)*pow(n_y, 2) - 5*pow(n_x, 2)*n_y + n_x*pow(n_y, 3) - 5*n_x*pow(n_y, 2) + 12*n_x*n_y - 3*n_x - 3*n_y));
}

double M5two(double m2, double m3, double m5, int n_x, int n_y)
{
  return pow(n_x, 2)*pow(n_y, 2)*(10*m2*m3*(pow(n_x, 2) + 2*n_x*n_y + pow(n_y, 2))*(-2*pow(n_x, 3)*pow(n_y, 3) + 5*pow(n_x, 3)*pow(n_y, 2) - 8*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 5*pow(n_x, 2)*pow(n_y, 3) - 8*n_x*pow(n_y, 3) + 4*pow(n_y, 3)) + m5*n_x*n_y*(n_x + n_y)*(pow(n_x, 4)*pow(n_y, 2) + 2*pow(n_x, 3)*pow(n_y, 3) - 12*pow(n_x, 3)*pow(n_y, 2) + 2*pow(n_x, 3)*n_y + pow(n_x, 2)*pow(n_y, 4) - 12*pow(n_x, 2)*pow(n_y, 3) + 60*pow(n_x, 2)*pow(n_y, 2) - 42*pow(n_x, 2)*n_y + 20*pow(n_x, 2) + 2*n_x*pow(n_y, 3) - 42*n_x*pow(n_y, 2) + 20*pow(n_y, 2)))/(10*(pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 5*pow(n_x, 2)*n_y - 2*pow(n_x, 2) + 5*n_x*pow(n_y, 2) - 2*pow(n_y, 2))*(-2*pow(n_x, 3)*pow(n_y, 3) + 5*pow(n_x, 3)*pow(n_y, 2) - 8*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 5*pow(n_x, 2)*pow(n_y, 3) - 8*n_x*pow(n_y, 3) + 4*pow(n_y, 3)) + (pow(n_x, 4)*pow(n_y, 3) + pow(n_x, 3)*pow(n_y, 4) - 10*pow(n_x, 3)*pow(n_y, 3) + 10*pow(n_x, 3)*pow(n_y, 2) - 10*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 10*pow(n_x, 2)*pow(n_y, 3) - 10*n_x*pow(n_y, 3) + 4*pow(n_y, 3))*(pow(n_x, 4)*pow(n_y, 2) + 2*pow(n_x, 3)*pow(n_y, 3) - 12*pow(n_x, 3)*pow(n_y, 2) + 2*pow(n_x, 3)*n_y + pow(n_x, 2)*pow(n_y, 4) - 12*pow(n_x, 2)*pow(n_y, 3) + 60*pow(n_x, 2)*pow(n_y, 2) - 42*pow(n_x, 2)*n_y + 20*pow(n_x, 2) + 2*n_x*pow(n_y, 3) - 42*n_x*pow(n_y, 2) + 20*pow(n_y, 2)));
}

double M2M3two(double m2, double m3, double m5, int n_x, int n_y)
{
  return pow(n_x, 2)*pow(n_y, 2)*(m2*m3*(pow(n_x, 2) + 2*n_x*n_y + pow(n_y, 2))*(pow(n_x, 4)*pow(n_y, 3) + pow(n_x, 3)*pow(n_y, 4) - 10*pow(n_x, 3)*pow(n_y, 3) + 10*pow(n_x, 3)*pow(n_y, 2) - 10*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 10*pow(n_x, 2)*pow(n_y, 3) - 10*n_x*pow(n_y, 3) + 4*pow(n_y, 3)) - m5*n_x*n_y*(n_x + n_y)*(pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 5*pow(n_x, 2)*n_y - 2*pow(n_x, 2) + 5*n_x*pow(n_y, 2) - 2*pow(n_y, 2)))/(10*(pow(n_x, 3)*pow(n_y, 2) + pow(n_x, 2)*pow(n_y, 3) - 8*pow(n_x, 2)*pow(n_y, 2) + 5*pow(n_x, 2)*n_y - 2*pow(n_x, 2) + 5*n_x*pow(n_y, 2) - 2*pow(n_y, 2))*(-2*pow(n_x, 3)*pow(n_y, 3) + 5*pow(n_x, 3)*pow(n_y, 2) - 8*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 5*pow(n_x, 2)*pow(n_y, 3) - 8*n_x*pow(n_y, 3) + 4*pow(n_y, 3)) + (pow(n_x, 4)*pow(n_y, 3) + pow(n_x, 3)*pow(n_y, 4) - 10*pow(n_x, 3)*pow(n_y, 3) + 10*pow(n_x, 3)*pow(n_y, 2) - 10*pow(n_x, 3)*n_y + 4*pow(n_x, 3) + 10*pow(n_x, 2)*pow(n_y, 3) - 10*n_x*pow(n_y, 3) + 4*pow(n_y, 3))*(pow(n_x, 4)*pow(n_y, 2) + 2*pow(n_x, 3)*pow(n_y, 3) - 12*pow(n_x, 3)*pow(n_y, 2) + 2*pow(n_x, 3)*n_y + pow(n_x, 2)*pow(n_y, 4) - 12*pow(n_x, 2)*pow(n_y, 3) + 60*pow(n_x, 2)*pow(n_y, 2) - 42*pow(n_x, 2)*n_y + 20*pow(n_x, 2) + 2*n_x*pow(n_y, 3) - 42*n_x*pow(n_y, 2) + 20*pow(n_y, 2)));
}

