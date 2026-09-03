#ifndef _POTENTIAL_H_
#define _POTENTIAL_H_
#include <vector>
#include "cubic_spline.h"
using namespace std;

class EAM
{
public:
    double F_M(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double F_M_deriv(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double F_H(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double F_H_deriv(double rho, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double f_M(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double f_M_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double f_H(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double f_H_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_M(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_M_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_MH(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_MH_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_H(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
    double phi_H_deriv(double r, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);
};
#endif
