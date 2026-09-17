#ifndef _MINIMIZER_H_
#define _MINIMIZER_H_
#include <petsctao.h>
#include <stdio.h>
#include <vector>
#include "cubic_spline.h"
struct Vec3D;
struct Input;
using namespace std;

struct AppCtx
{
  Input *input;
  double mubd;
  vector<vector<int> > *MM;
  vector<vector<int> > *MH;
  vector<vector<int> > *HM;
  vector<vector<int> > *HH;
  int numDOFPerSite;
  int globDim;
  int locDim;
  int first_site, last_site_plus_one;
  double *alldof;
  vector<vector<Vec3D> > *QPp;
  double QWp;
	vector<CubicSpline> *CSembed;
 	vector<CubicSpline> *CSrho;
 	vector<CubicSpline> *CSpair;


    // by XS
    double *send_M, *recv_M;
    double *send_H, *recv_H;

  // the following vectors will be populated (repeatedly) when TAO calls "formFunctionAndGradient"
  vector<Vec3D> q_M0;
  vector<Vec3D> q_H0;
  vector<double> sigma_M0;
  vector<double> sigma_H0;
  vector<double> x0;
  vector<Vec3D> Fq_M0;
  vector<Vec3D> Fq_H0;
  vector<double> Fsigma_M0;
  vector<double> Fsigma_H0;
  vector<double> Fx0;
vector<int> rigidAtoms; // rigidAtoms[0]=A, rigidAtoms[1]=B, rigidAtoms[2]=C
  AppCtx(Input *input_, vector<vector<int> > *MM_, vector<vector<int> > *MH_,
         vector<vector<int> > *HM_, vector<vector<int> > *HH_, vector<Vec3D> *q_M0_, 
         vector<Vec3D> *q_H0_, vector<double> *sigma_M0_, vector<double> *sigma_H0_, 
         vector<double> *x0_, vector<int> *rigidAtoms_,  
         vector<CubicSpline> *CSembed_, vector<CubicSpline> *CSrho_, vector<CubicSpline> *CSpair_,
				 int numDOFPerSite_, int globDim_, int locDim_, int first_site_, 
         int last_site_plus_one_, double *alldof_, vector<vector<Vec3D> > *QPp_, double QWp_); 
  ~AppCtx();
}; 

class Minimizer
{
  MPI_Comm      comm;
  PetscMPIInt   size, rank;
  Tao           tao;
  AppCtx        *application_context; 
  Vec           yy;    /* all-inclusive solution vector */
  Vec           Ub, Lb; // variable bounds
 
  int numDOFPerSite;
  int globDim; // total # of dofs (not sites) to be solved
  int locDim; // # of dofs (not sites) to be handled by this CPU core
  int first_site, last_site_plus_one; //global index of first element, and last_plus_one
  double *alldof;
  
  void updateApplicationContext(double mubd_, vector<Vec3D> *q_M0_, vector<Vec3D> *q_H0_, 
                                vector<double> *sigma_M0_, vector<double> *sigma_H0_,
                                vector<double> *x0_);
  static PetscErrorCode applicationData2Petsc(void* ptr/*AppCtx*/,
                                              vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, 
                                              vector<double> &sigma_M0,
                                              vector<double> &sigma_H0, 
                                              vector<double> &x0, Vec &yy);
  static PetscErrorCode petscData2Application(void* ptr/*AppCtx*/,
                                              Vec &yy, vector<Vec3D> &q_M, vector<Vec3D> &q_H, 
                                              vector<double> &sigma_M, vector<double> &sigma_H, 
                                              vector<double> &x0);
  static PetscErrorCode formFunctionAndGradient(Tao, Vec/*input vector*/, PetscReal* /*func value*/,
                                         Vec/*gradient*/, void* /*AppCtx*/);
  static PetscErrorCode formVariableBounds(Tao, Vec Ub, Vec Lb, void* /*AppCtx*/);

public:
  Minimizer(int *argc, char ***argv);
  ~Minimizer();
  PetscErrorCode initialize(Input *input, vector<vector<int> > *MM, vector<vector<int> > *MH,
                  vector<vector<int> > *HM, vector<vector<int> > *HH, 
                  vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, vector<double> &sigma_M0,
                  vector<double> &sigma_H0, vector<double> &x0, vector<int> &rigidAtoms, vector<vector<Vec3D> > *QPp, double &QWp,
									vector<CubicSpline> *CSembed, vector<CubicSpline> *CSrho, vector<CubicSpline> *CSpair); 
  PetscErrorCode minimizeFreeEntropy(double mubd, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, 
                          vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0, 
                          vector<Vec3D> &q_M, vector<Vec3D> &q_H, vector<double> &sigma_M, 
                          vector<double> &sigma_H, vector<double> &x, double &obj_fun);
  PetscErrorCode calculateFormationEnergy(vector<double> &x0, double gammabd, vector<double> &gamma, vector<double> &f, vector<int> &full_H);
  PetscErrorCode calculateMeanPotentialEnergy(vector<double> &x0, vector<double> &V_M, vector<double> &V_H);
	PetscErrorCode calculateLocalStress(vector<double> &x0, vector<vector<double> > &pi_M, vector<vector<double> > &pi_H);
  PetscErrorCode findLocalNeighbors(vector<vector<int> > &MM_ext, vector<vector<int> > &MH_ext,
                  vector<vector<int> > &HM_ext, vector<vector<int> > &HH_ext,
                  vector<vector<int> > &MM, vector<vector<int> > &MH,
                  vector<vector<int> > &HM, vector<vector<int> > &HH,
                  int &maxNeib);
};
#endif
