#ifndef CUBIC_SPLINE_HPP
#define CUBIC_SPLINE_HPP

/**
 * Nonuniform cubic splines with n intervals
 */

typedef struct CubicSpline
{
  int     n;
  double dx;
  double* x;
  double* a;
  double* b;
  double* c;
  double* d;
} CubicSpline;

void init_spline( CubicSpline* cs, int n, double dx, double *y );
void destroy_spline( CubicSpline* cs );

double cubic_spline( CubicSpline* cs, double x );
double d_cubic_spline( CubicSpline* cs, double x );
double d2_cubic_spline( CubicSpline* cs, double x );

double mf_cubic_spline( CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b);
void mf_d_cubic_spline(  double mf[3], double *mfwi, double *mfwj, CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b);
void mf_d_cubic_spline(  double mf[3], CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b);

#endif
