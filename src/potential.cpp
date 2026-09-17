/*
- EAM potential for palladium hydride
- data obtained from the appendix of X.W. Zhou et al. 2007
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "input.h"
#include "potential.h"
#include "cubic_spline.h"
#include "eam_setfl.h"

using namespace std;

Input input;
const double pi = 3.141592653589793;
//double rs = input.file.rc*0.8; // cutoff distance for smoothing pair functions

double EAM::F_M(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // embedding energy for M 
{
	return cubic_spline(&CSembed[0], rho);
}

double EAM::F_M_deriv(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
	return d_cubic_spline(&CSembed[0], rho);
}

double EAM::F_H(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // embedding energy for H
{
	return cubic_spline(&CSembed[1], rho);
}

double EAM::F_H_deriv(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
	return d_cubic_spline(&CSembed[1], rho);
}

double EAM::f_M(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // electron density for M 
{
  double rs = 5.0; // This rs only works for PdH.
  if (r<rs)
	  return cubic_spline(&CSrho[0], r);
  else { 
    double funRc = cubic_spline(&CSrho[0], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSrho[0], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth f_M. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funRc*exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::f_M_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
  double rs = 5.0;
  if (r<rs)
    return d_cubic_spline(&CSrho[0], r);
  else {
    double funRc = cubic_spline(&CSrho[0], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSrho[0], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth f_M_deriv. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funDerivRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::f_H(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // electron density for H 
{
  double rs = 5.2;
  if (r<rs)
    return cubic_spline(&CSrho[1], r);
  else {
    double funRc = cubic_spline(&CSrho[1], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSrho[1], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth f_H. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::f_H_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
  double rs = 5.2;
  if (r<rs)
    return d_cubic_spline(&CSrho[1], r);
  else {
    double funRc = cubic_spline(&CSrho[1], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSrho[1], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth f_H_deriv. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funDerivRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::phi_M(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // pair energy for M
{
  double rs = 5.05;
  if (r<rs)
    return cubic_spline(&CSpair[0], r);
  else {
    double funRc = cubic_spline(&CSpair[0], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[0], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_M. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::phi_M_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
  double rs = 5.05;
  if (r<rs)
    return d_cubic_spline(&CSpair[0], r);
  else {
    double funRc = cubic_spline(&CSpair[0], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[0], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_M_deriv. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funDerivRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::phi_MH(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // pair energy for MH
{
  double rs = 4.8;
  if (r<rs)
    return cubic_spline(&CSpair[1], r);
  else {
    double funRc = cubic_spline(&CSpair[1], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[1], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_MH. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::phi_MH_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
  double rs = 4.8;
  if (r<rs)
    return d_cubic_spline(&CSpair[1], r);
  else {
    double funRc = cubic_spline(&CSpair[1], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[1], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_MH_deriv. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funDerivRc * exp(funDerivRc/funRc*(r-rs));
  }
}

double EAM::phi_H(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair) // pair energy for H
{
  double rs = 5.0;
  if (r<rs)
    return cubic_spline(&CSpair[2], r);
  else {
    double funRc = cubic_spline(&CSpair[2], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[2], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_H. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funRc * exp(funDerivRc/funRc*(r-rs));
  } 
}

double EAM::phi_H_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
  double rs = 5.0;
  if (r<rs)
    return d_cubic_spline(&CSpair[2], r);
  else {
    double funRc = cubic_spline(&CSpair[2], rs); // function value at r = rc-rs
    double funDerivRc = d_cubic_spline(&CSpair[2], rs); // gradient value at r = rc-rs 
    if (!(funRc*funDerivRc<0.0)) {
      printf("Warning: fail to smooth phi_H_driv. Function value: %lf. Gradient value: %lf.\n", funRc, funDerivRc);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    }
    return funDerivRc * exp(funDerivRc/funRc*(r-rs));
  }
}

