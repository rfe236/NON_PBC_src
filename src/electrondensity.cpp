/*
- calculate averaged electron density of one M or H atom
*/
#include <iostream>
#include <vector>
#include <math.h>
#include "Vector3D.h"
#include "potential.h"
#include "cubic_spline.h"
using namespace std;

EAM eam;
inline Vec3D minimumImage(Vec3D dr, double Lx, double Ly, double Lz)
{
    dr[0] -= Lx * round(dr[0] / Lx);
    dr[1] -= Ly * round(dr[1] / Ly);
    dr[2] -= Lz * round(dr[2] / Lz);

    return dr;
}

double electronDensity(char atype, int n, double rc, double Lx, double Ly, double Lz, vector<vector<int> > &MM, 
                    vector<vector<int> > &MH, vector<vector<int> > &HM, 
                    vector<vector<int> > &HH, vector<Vec3D> &q_M, vector<Vec3D> &q_H, 
                    vector<double> &sigma_M, vector<double> &sigma_H, vector<double> &x, 
                    vector<vector<Vec3D> > &QPp, double &QWp, 
										vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
    double rho = 0.0;
    if (atype=='P') {
        for(int i=0; i<MM[n].size(); i++) {
            for(int k=0; k<QPp.size(); k++) {
                Vec3D r = q_M[n]+sqrt(2.0)*sigma_M[n]*QPp[k][0]-q_M[MM[n][i]]-sqrt(2.0)*sigma_M[MM[n][i]]*QPp[k][1];
r = minimumImage(r, Lx, Ly, Lz);
//rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;
                rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;    
						}
        }

        for(int i=0; i<MH[n].size(); i++) {
            for(int k=0; k<QPp.size(); k++) {
                Vec3D r = q_M[n]+sqrt(2.0)*sigma_M[n]*QPp[k][0]-q_H[MH[n][i]]-sqrt(2.0)*sigma_H[MH[n][i]]*QPp[k][1];
r = minimumImage(r, Lx, Ly, Lz);
//rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;
                rho += x[MH[n][i]]*eam.f_H(r.norm(), CSembed, CSrho, CSpair)*QWp;
            }
        }
    }
    else if (atype=='H') {
        for(int i=0; i<HM[n].size(); i++) {
            for(int k=0; k<QPp.size(); k++) {
                Vec3D r = q_H[n]+sqrt(2.0)*sigma_H[n]*QPp[k][0]-q_M[HM[n][i]]-sqrt(2.0)*sigma_M[HM[n][i]]*QPp[k][1];
r = minimumImage(r, Lx, Ly, Lz);
//rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;
                rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;
            }
        }

        for(int i=0; i<HH[n].size(); i++) {
            for(int k=0; k<QPp.size(); k++) {
                Vec3D r = q_H[n]+sqrt(2.0)*sigma_H[n]*QPp[k][0]-q_H[HH[n][i]]-sqrt(2.0)*sigma_H[HH[n][i]]*QPp[k][1];
r = minimumImage(r, Lx, Ly, Lz);
//rho += eam.f_M(r.norm(), CSembed, CSrho, CSpair)*QWp;
                rho += x[HH[n][i]]*eam.f_H(r.norm(), CSembed, CSrho, CSpair)*QWp;
            }
        }
    }
    else {
        cerr << "Error: First input in electronDensity must be either 'P' or 'H'." << endl; cout.flush();
        exit(-1);
    }

    return rho;
}


