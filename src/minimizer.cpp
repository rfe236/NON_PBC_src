#include <iostream>
#include <math.h>
#include <algorithm>
#include "Vector3D.h"
#include "minimizer.h"
#include "input.h"
#include "potential.h"
#include "KDTree.h"
#include "cubic_spline.h"
#include <fstream>
using namespace std;
inline Vec3D minimumImage(Vec3D dr, double Lx, double Ly, double Lz)
{
    dr[0] -= Lx * round(dr[0] / Lx);
    dr[1] -= Ly * round(dr[1] / Ly);
    dr[2] -= Lz * round(dr[2] / Lz);

    return dr;
}

inline Vec3D nearestImage(Vec3D p, Vec3D center, double Lx, double Ly, double Lz)
{
    Vec3D dr(p[0]-center[0], p[1]-center[1], p[2]-center[2]);
    dr = minimumImage(dr, Lx, Ly, Lz);
    return Vec3D(center[0]+dr[0], center[1]+dr[1], center[2]+dr[2]);
}
static char help[] = "minimization of grand-canonical free entropy\n";

//// Slide atom "p" to its nearest periodic copy, as seen from "center".


double electronDensity(char atype, int n, double rc, double Lx, double Ly, double Lz, vector<vector<int> > &MM, 
                    vector<vector<int> > &MH, vector<vector<int> > &HM, 
                    vector<vector<int> > &HH, vector<Vec3D> &q_M, vector<Vec3D> &q_H, 
                    vector<double> &sigma_M, vector<double> &sigma_H, vector<double> &x, 
                    vector<vector<Vec3D> > &QPp, double &QWp, 
										vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair);

AppCtx::AppCtx(Input *input_, vector<vector<int> > *MM_, vector<vector<int> > *MH_,
               vector<vector<int> > *HM_, vector<vector<int> > *HH_, vector<Vec3D> *q_M0_, 
               vector<Vec3D> *q_H0_, vector<double> *sigma_M0_, vector<double> *sigma_H0_, 
               vector<double> *x0_, vector<int> *rigidAtoms_,  
							 vector<CubicSpline> *CSembed_, vector<CubicSpline> *CSrho_, vector<CubicSpline> *CSpair_,
							 int numDOFPerSite_, int globDim_, int locDim_, 
               int first_site_, int last_site_plus_one_, double *alldof_, 
               vector<vector<Vec3D> > *QPp_, double QWp_) : 
  input(input_), MM(MM_), MH(MH_), HM(HM_), HH(HH_), QPp(QPp_), CSembed(CSembed_), CSrho(CSrho_), CSpair(CSpair_), send_M(NULL),recv_M(NULL), send_H(NULL), recv_H(NULL)

{

  numDOFPerSite = numDOFPerSite_;
  globDim = globDim_;
  locDim = locDim_;
  first_site = first_site_;
  last_site_plus_one = last_site_plus_one_;
  alldof = alldof_;
  QWp = QWp_; 
  // allocate memory
  q_M0      = *q_M0_;
  q_H0       = *q_H0_;
  sigma_M0  = *sigma_M0_;
  sigma_H0   = *sigma_H0_;
  x0 = *x0_;
rigidAtoms = *rigidAtoms_;
  Fq_M0     = *q_M0_;
  Fq_H0      = *q_H0_;
  Fsigma_M0 = *sigma_M0_;
  Fsigma_H0  = *sigma_H0_;
  Fx0 = *x0_;
  
    // by XS
    // allocate memory for  electron density of each site.
    send_M = new double[q_M0.size()];
    recv_M = new double[q_M0.size()];
    send_H = new double[q_H0.size()];
    recv_H = new double[q_H0.size()];
}

AppCtx::~AppCtx() {
    delete[] send_M;
    delete[] recv_M;
    delete[] send_H;
    delete[] recv_H;
}

Minimizer::Minimizer(int *argc, char ***argv)
{
  PetscErrorCode ierr;
  PetscInitialize(argc, argv, *argc>=3 ? (*argv)[2] : (char*)0, help); 
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  comm = PETSC_COMM_WORLD; //same as MPI_COMM_WORLD
  application_context = NULL;
  alldof = NULL;
}


Minimizer::~Minimizer()
{
  if(alldof) delete[] alldof;
  PetscErrorCode ierr;
  ierr = TaoDestroy(&tao);
  ierr = VecDestroy(&yy);
  ierr = VecDestroy(&Ub);
  ierr = VecDestroy(&Lb);
}


PetscErrorCode Minimizer::formFunctionAndGradient(Tao tao, Vec X, PetscReal* fun, Vec Grad, void* ptr)
{
  /* TODO: calculate f and Grad at X */
  // THIS IS THE MOST IMPORTANT FUNCTION
  // THIS FUNCTION MUST BE PARALLELIZED, AND EFFICIENT!
  
  PetscErrorCode ierr = 0;
//PetscPrintf(PETSC_COMM_WORLD, "Inside H loop\n");
  AppCtx *app = (AppCtx *) ptr;
 int myrank;
MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  EAM eam; // by XS
  const double kB = 8.6173324e-5;
  const double hbar = 6.582119596e-16; // [eV*s]
  const double mH = 1.6735328115271e-27; // [kg] atomic mass of H
  const double mM = 1.7671455e-25; // [kg] atomic mass of M
  const double eVJ = 1.602176634e-19; // [J/eV] eV to J
  const double Am = 1.0e-10; // [m/A] Angstrom to meter
  double &T = app->input->file.T0;
  double &mubd = app->mubd;
double a_M = app->input->file.a_M;
double Fxx = app->input->file.Fxx_current;
double Lx = Fxx * app->input->file.N * a_M;   // compressed
//double Lx = app->input->file.N * a_M;
double Ly = app->input->file.N * a_M;
double Lz = app->input->file.N * a_M;
  int nM = app->q_M0.size();
  int nH  = app->q_H0.size();
  vector<vector<int> > &MM = *(app->MM); //create reference
  vector<vector<int> > &MH = *(app->MH); 
  vector<vector<int> > &HM = *(app->HM); 
  vector<vector<int> > &HH = *(app->HH);
  vector<vector<Vec3D> > &QPp = *(app->QPp);
  double &QWp = app->QWp;
	vector<CubicSpline> &CSembed = *(app->CSembed);
 	vector<CubicSpline> &CSrho = *(app->CSrho);
 	vector<CubicSpline> &CSpair = *(app->CSpair);
  double fun_local = 0.0;
  double sum_V = 0.0; // potential energy
  double sum_sigma = 0.0; // vibrational energy
  double sum_conf = 0.0; // configurational entropy 
  double sum_x = 0.0; // sum of atomic fractions

  // populate q_M0, q_H0, ..., in the application context
  ierr = petscData2Application(app, X, app->q_M0, app->q_H0, app->sigma_M0, app->sigma_H0, app->x0);
  CHKERRQ(ierr);

  // clear force vectors
    for(int i=0; i<nM; i++) {
        app->Fq_M0[i] = 0.0;
        app->Fsigma_M0[i] = 0.0;
    }
    for(int i=0; i<nH; i++) {
        app->Fq_H0[i] = 0.0;
        app->Fsigma_H0[i] = 0.0;
        app->Fx0[i] = 0.0;
    } 

    // initiallization for electron density calculation
    for(int i=0; i<nM; i++)
        app->send_M[i] = 0.0;
    for(int i=0; i<nH; i++)
        app->send_H[i] = 0.0;

/*
fun_local and app->F should be revised according to EAM potential.
fun_local is the total free entropy and app->Fq is its gradient with respect to positions.
app->Fsigma is the gradient with respect to sigma.
*/

  // ------------------------------------------------------------
  // Objective function = EAM potential, revised by XS.
  // ------------------------------------------------------------i
    // ------------------------------------------------------------
    // Parallelization I: calculate electron density of each site.
    // ------------------------------------------------------------
    if(app->first_site<nM) { // the first local site is a M
        for(int i=app->first_site; i<min(nM, app->last_site_plus_one); i++) {
            app->send_M[i] = electronDensity('P', i, app->input->file.rc, Lx, Ly, Lz, MM, MH, HM, HH, app->q_M0, app->q_H0, app->sigma_M0, app->sigma_H0, app->x0, QPp, QWp, CSembed, CSrho, CSpair);
				}
        for(int ii=nM; ii<app->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nM; // index in the H vectors
            app->send_H[i] = electronDensity('H', i, app->input->file.rc, Lx, Ly, Lz, MM, MH, HM, HH, app->q_M0, app->q_H0, app->sigma_M0, app->sigma_H0, app->x0, QPp, QWp, CSembed, CSrho, CSpair);
       }
    }
    else { //all the local sites are H
        for(int ii=app->first_site; ii<app->last_site_plus_one; ii++) {
            int i = ii - nM; // index in the H vectors
            app->send_H[i] = electronDensity('H', i, app->input->file.rc, Lx, Ly, Lz,MM, MH, HM, HH, app->q_M0, app->q_H0, app->sigma_M0, app->sigma_H0, app->x0, QPp, QWp, CSembed, CSrho, CSpair);
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(app->send_M, app->recv_M, nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(app->send_H, app->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // ------------------------------------------------------------
    // Parallelization II: calculate objective value and gradients.
    // ------------------------------------------------------------
    if(app->first_site<nM) { // the first local site is a M
        for(int i=app->first_site; i<min(nM, app->last_site_plus_one); i++) {
        // calc. force for Site i
            Vec3D &qi = app->q_M0[i];
            double &sigmai = app->sigma_M0[i];

            // Gaussian quadrature for embedding energy
            sum_V += eam.F_M(app->recv_M[i], CSembed, CSrho, CSpair);

            // M-M pairs
            for(int j=0; j<MM[i].size(); j++) {
                int &jj = MM[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_M0[jj];
                double &sigmaj = app->sigma_M0[jj];
 
                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_M0[i] += QWp*eam.F_M_deriv(app->recv_M[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += QWp*eam.F_M_deriv(app->recv_M[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_M0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    
                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*QWp*eam.phi_M(r.norm(), CSembed, CSrho, CSpair);
                    
                    // calculate gradient for pair energy
                    app->Fq_M0[i] += QWp*eam.phi_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += QWp*eam.phi_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // M-H pairs
            for(int j=0; j<MH[i].size(); j++) {
                int &jj = MH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_M0[i] += app->x0[jj]*QWp*eam.F_M_deriv(app->recv_M[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += app->x0[jj]*QWp*eam.F_M_deriv(app->recv_M[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_M0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*app->x0[jj]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);

                    // calculate gradient for pair energy
                    app->Fq_M0[i] += app->x0[jj]*QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_M0[i] += app->x0[jj]*QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_M0[i] *= 1.0/T;
            app->Fsigma_M0[i] *= 1.0/T*sqrt(2.0);
            app->Fsigma_M0[i] -= 3.0*kB/sigmai;
            sum_sigma += 1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mM*sigmai*sigmai*Am*Am))-1.0);
        }
        for(int ii=nM; ii<app->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nM; // index in the H vectors
            Vec3D &qi = app->q_H0[i]; // by XS
            double &sigmai = app->sigma_H0[i];

            // Gaussian quadrature for embedding energy
            sum_V += app->x0[i]*eam.F_H(app->recv_H[i], CSembed, CSrho, CSpair);
            app->Fx0[i] += eam.F_H(app->recv_H[i], CSembed, CSrho, CSpair); 

            // H-M pairs
            for(int j=0; j<HM[i].size(); j++) {
                int &jj = HM[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_M0[jj];
                double &sigmaj = app->sigma_M0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                     r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fx0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*app->x0[i]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);
                    app->Fx0[i] += QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);
                    // calculate gradient for pair energy
                    app->Fq_H0[i] += QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // H-H pairs
            for(int j=0; j<HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fx0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*app->x0[i]*app->x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);
                    app->Fx0[i] += app->x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_H0[i] *= 1.0/T*app->x0[i];
            app->Fsigma_H0[i] *= 1.0/T*sqrt(2.0)*app->x0[i];
           app->Fsigma_H0[i] -= 3.0*kB/sigmai*app->x0[i];
// {{
    //           int myrank;
      //          MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
        //        double grad_potential_raw = app->Fx0[i];
          //      double grad_vib  = 1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am)) - 1.0);
            //    double grad_conf = kB*log(app->x0[i]/(1.0-app->x0[i]));
              //  double grad_chem = -mubd/T;
            //    double grad_analytical = grad_potential_raw*(1.0/T) + grad_vib + grad_conf + grad_chem;
       //         if (myrank == 1 && i < 5) {
        //            cerr << "==============================" << endl;
       //             cerr << "H site: " << i << endl;
       //             cerr << "x value            = " << app->x0[i] << endl;
       //             cerr << "Potential part     = " << grad_potential_raw << endl;
       //             cerr << "Vibrational part   = " << grad_vib << endl;
      //              cerr << "Configurational    = " << grad_conf << endl;
      //              cerr << "Chemical part      = " << grad_chem << endl;
      //              cerr << "Analytical gradient= " << grad_analytical << endl;
      //              cerr << "==============================" << endl;
        //        }
        //    }
        //

if (app->x0[i] >= 0.9999) {
    // Treat as fully occupied 
      app->x0[i] = 1.0;
      app->Fx0[i] = 0.0;
          }
else {
            app->Fx0[i] *= 1.0/T;
            app->Fx0[i] += (1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am))-1.0) + kB*log(app->x0[i]/(1.0-app->x0[i])) - mubd/T);
}
            sum_sigma += 1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am))-1.0)*app->x0[i];
            sum_conf += (app->x0[i]*log(app->x0[i])+(1.0-app->x0[i])*log(1.0-app->x0[i]));
            sum_x += app->x0[i];
        }
    } 
    else { //all the local sites are H
        for(int ii=app->first_site; ii<app->last_site_plus_one; ii++) {
            int i = ii - nM; // index in the H vectors
            Vec3D &qi = app->q_H0[i]; // by XS
            double &sigmai = app->sigma_H0[i];

            // Gaussian quadrature for embedding energy
            sum_V += app->x0[i]*eam.F_H(app->recv_H[i], CSembed, CSrho, CSpair);
            app->Fx0[i] += eam.F_H(app->recv_H[i], CSembed, CSrho, CSpair);

            // H-M pairs
            for(int j=0; j<HM[i].size(); j++) {
                int &jj = HM[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_M0[jj];
                double &sigmaj = app->sigma_M0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                     r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fx0[i] += QWp*eam.F_M_deriv(app->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*app->x0[i]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);
                    app->Fx0[i] += QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }

            // H-H pairs
            for(int j=0; j<HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = app->q_H0[jj];
                double &sigmaj = app->sigma_H0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // calculate gradient for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                     r = minimumImage(r, Lx, Ly, Lz);
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;                    app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
             app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                    app->Fx0[i] += app->x0[jj]*QWp*eam.F_H_deriv(app->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_V += 0.5*app->x0[i]*app->x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);
                    app->Fx0[i] += app->x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);

                    // calculate gradient for pair energy
                    app->Fq_H0[i] += app->x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
                    app->Fsigma_H0[i] += app->x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r*QPp[k][0];
                }
                //NOTE: although the pair force to site j is the same, we MUST NOT apply 
                //the force to site j because site j might be assigned to another processor.
            }
            app->Fq_H0[i] *= 1.0/T*app->x0[i];
            app->Fsigma_H0[i] *= 1.0/T*sqrt(2.0)*app->x0[i];
            app->Fsigma_H0[i] -= 3.0*kB/sigmai*app->x0[i];
  // {
          //      double grad_potential_raw = app->Fx0[i];
            //    double grad_vib  = 1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am)) - 1.0);
      //          double grad_conf = kB*log(app->x0[i]/(1.0-app->x0[i]));
       //         double grad_chem = -mubd/T;
       //         double grad_analytical = grad_potential_raw*(1.0/T) + grad_vib + grad_conf + grad_chem;
        //        if (myrank == 1 && i < 5) {
       //             cerr << "==============================" << endl;
      //              cerr << "H site: " << i << endl;
     //               cerr << "x value            = " << app->x0[i] << endl;
      //              cerr << "Potential part     = " << grad_potential_raw << endl;
       //             cerr << "Vibrational part   = " << grad_vib << endl;
       //             cerr << "Configurational    = " << grad_conf << endl;
       //             cerr << "Chemical part      = " << grad_chem << endl;
      //              cerr << "Analytical gradient= " << grad_analytical << endl;
      //              cerr << "==============================" << endl;
       //         }
       //     }
            app->Fx0[i] *= 1.0/T;
            app->Fx0[i] += (1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am))-1.0) + kB*log(app->x0[i]/(1.0-app->x0[i])) - mubd/T);            
            sum_sigma += 1.5*kB*(log(hbar*hbar*eVJ/(kB*T*mH*sigmai*sigmai*Am*Am))-1.0)*app->x0[i];
            sum_conf += (app->x0[i]*log(app->x0[i])+(1.0-app->x0[i])*log(1.0-app->x0[i]));
            sum_x += app->x0[i];
        }
    }
   
    fun_local = 1.0/T*sum_V + sum_sigma + kB*sum_conf - mubd/T*sum_x;

  // ------------------------------------------------------------
  // Revision end.
  // ------------------------------------------------------------

  // copy forces to Grad
  ierr = applicationData2Petsc(app, app->Fq_M0, app->Fq_H0, app->Fsigma_M0, app->Fsigma_H0, app->Fx0, Grad);
  CHKERRQ(ierr);

// Analytical objectgive and gradient norm
//static int iter_counter = 0;

//double grad_sum = 0.0;

//for (int i = 0; i < nM; i++)
//{
//    grad_sum += app->Fq_M0[i][0] * app->Fq_M0[i][0];
//    grad_sum += app->Fq_M0[i][1] * app->Fq_M0[i][1];
//    grad_sum += app->Fq_M0[i][2] * app->Fq_M0[i][2];

//    grad_sum += app->Fsigma_M0[i] * app->Fsigma_M0[i];
//}
//for (int i = 0; i < nH; i++)
//{
//    grad_sum += app->Fq_H0[i][0] * app->Fq_H0[i][0];
//    grad_sum += app->Fq_H0[i][1] * app->Fq_H0[i][1];
//    grad_sum += app->Fq_H0[i][2] * app->Fq_H0[i][2];

//    grad_sum += app->Fsigma_H0[i] * app->Fsigma_H0[i];

//    grad_sum += app->Fx0[i] * app->Fx0[i];
//}

//double grad_norm = sqrt(grad_sum);

//cout << "ITER " << iter_counter
  //   << " OBJECTIVE " << fun_local
 //    << " GRAD_NORM " << grad_norm
 //    << endl;

//iter_counter++;












//if (myrank == 0) {
//  for (int i = 0; i < (int)app->Fq_M0.size(); i++) {
//    cerr << " " << app->Fq_M0[i][0] << endl;
//     cerr << " " << app->Fq_M0[i][1] << endl;
//   cerr << "  " << app->Fq_M0[i][2] << endl;
//  cerr << " " << app->Fsigma_M0[i] << endl;
//   cerr << " " << 0.0 << endl;
// }
//}
//if (myrank == 0) {
 // for (int i = 0; i < (int)app->Fq_H0.size(); i++) {
// cerr << " " << app->Fq_H0[i][0] << endl;
//  cerr << " " << app->Fq_H0[i][1] << endl;
//cerr << " " << app->Fq_H0[i][2] << endl;
// cerr << " " << app->Fsigma_H0[i] << endl;   
//cerr << " " << app->Fx0[i] << endl;
//  }                                                          

//}

  //VecView(Grad,PETSC_VIEWER_STDOUT_WORLD);

  // sum up fun from all the processors
  ierr = (PetscErrorCode)MPI_Allreduce((void*)&fun_local, (void*)fun, 1, MPIU_REAL, MPIU_SUM, 
                                       MPI_COMM_WORLD);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


void Minimizer::updateApplicationContext(double mubd_, vector<Vec3D> *q_M0_, vector<Vec3D> *q_H0_,
                                         vector<double> *sigma_M0_, vector<double> *sigma_H0_, 
                                         vector<double> *x0_)

{
  application_context->mubd = mubd_;

  application_context->q_M0     = *q_M0_;
  application_context->q_H0     = *q_H0_;
  application_context->sigma_M0 = *sigma_M0_;
  application_context->sigma_H0 = *sigma_H0_;
  application_context->x0       = *x0_;
}



PetscErrorCode Minimizer::applicationData2Petsc(void* ptr, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, 
                                      vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0, Vec &yy)
{
  PetscErrorCode ierr; 
  AppCtx *app = (AppCtx *) ptr;

  int nM = q_M0.size();
  int nH  = q_H0.size();
  int counter = 0;

  // create locData and loc2glob
  PetscReal locData[app->locDim];
  PetscInt  loc2glob[app->locDim];
  for(int i=0; i<app->locDim; i++) 
    loc2glob[i] = app->first_site*app->numDOFPerSite + i;

  // fill locData
  if(app->first_site<nM) { // the first site is an M
    for(int i=app->first_site; i<min(nM, app->last_site_plus_one); i++) {
      locData[counter++] = q_M0[i][0];
      locData[counter++] = q_M0[i][1];
      locData[counter++] = q_M0[i][2];
      locData[counter++] = sigma_M0[i];
      locData[counter++] = 0.0; // The atomic fraction at M sites is not updated
    }
    for(int i=nM; i<app->last_site_plus_one; i++) {
      locData[counter++] = q_H0[i-nM][0];
      locData[counter++] = q_H0[i-nM][1];
      locData[counter++] = q_H0[i-nM][2];
      locData[counter++] = sigma_H0[i-nM];
      locData[counter++] = x0[i-nM];
    }
  } else { // all the sites belong to H
    for(int i=app->first_site; i<app->last_site_plus_one; i++) {
      locData[counter++] = q_H0[i-nM][0];
      locData[counter++] = q_H0[i-nM][1];
      locData[counter++] = q_H0[i-nM][2];
      locData[counter++] = sigma_H0[i-nM];
      locData[counter++] = x0[i-nM];
    }
  }

  if(counter!=app->locDim) {
    PetscPrintf(MPI_COMM_WORLD, "Error in applicationData2Petsc!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // copy locData to Petsc vector 
  ierr = VecSetValues(yy, app->locDim, loc2glob, locData, INSERT_VALUES);
  ierr = VecAssemblyBegin(yy);
  ierr = VecAssemblyEnd(yy);

  return ierr;
}


PetscErrorCode Minimizer::petscData2Application(void* ptr,
                                                Vec &yy, vector<Vec3D> &q_M, vector<Vec3D> &q_H, 
                                                vector<double> &sigma_M, vector<double> &sigma_H,
                                                vector<double> &x)
{
  PetscErrorCode ierr;
  AppCtx *app = (AppCtx *)ptr;

  PetscScalar *locCopy;
  ierr = VecGetArray(yy, &locCopy);

  for(int i=0; i<app->globDim; i++)
    (app->alldof)[i] = 0.0;

//  double recv_alldof[app->globDim];
    double *recv_alldof;
    recv_alldof = new double[app->globDim];

  for(int i=0; i<app->globDim; i++)
    recv_alldof[i] = 0.0;
 
  for(int i=0; i<app->locDim; i++)
    (app->alldof)[app->first_site*app->numDOFPerSite+i] = locCopy[i];

  // collect data
  MPI_Allreduce(app->alldof, recv_alldof, app->globDim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  ierr = VecRestoreArray(yy, &locCopy);

  // copy recv_alldof to application data 
  int counter = 0;
  for(int i=0; i<q_M.size(); i++) { //M sites
    for(int j=0; j<3; j++)
      q_M[i][j] = recv_alldof[counter++];
    sigma_M[i] = recv_alldof[counter++];
    counter++; // skip atomic fraction, as we don't need to deal with it at the M sites
  }
  for(int i=0; i<q_H.size(); i++) {
    for(int j=0; j<3; j++)
      q_H[i][j] = recv_alldof[counter++];
    sigma_H[i] = recv_alldof[counter++];
    x[i] = recv_alldof[counter++];
  }
  if(counter!=app->globDim) {
    PetscPrintf(MPI_COMM_WORLD, "Error in petscData2Application!\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  delete[] recv_alldof;

  return ierr;
}


PetscErrorCode Minimizer::initialize(Input *input, vector<vector<int> > *MM, 
                                     vector<vector<int> > *MH, vector<vector<int> > *HM, 
                                     vector<vector<int> > *HH, vector<Vec3D> &q_M0, 
                                     vector<Vec3D> &q_H0, vector<double> &sigma_M0, 
                                     vector<double> &sigma_H0, vector<double> &x0, vector<int> &rigidAtoms,
                                     vector<vector<Vec3D> > *QPp, double &QWp, 
																		 vector<CubicSpline> *CSembed, vector<CubicSpline> *CSrho, 
																		 vector<CubicSpline> *CSpair) 
{
  PetscErrorCode ierr;

  // split the solution vector over all the CPU cores
  numDOFPerSite = 5; // each site has 5 DOFs, i.e., three mean positions + one standard deviation + one atomic fraction

  int globSite = q_M0.size() + q_H0.size();
  int locSite = PETSC_DECIDE;
  PetscSplitOwnership(comm, &locSite, &globSite); //calculates locSite
  locDim = numDOFPerSite*locSite;
  globDim = PETSC_DECIDE;
  PetscSplitOwnership(comm, &locDim, &globDim); //calculates globDim
 
  // create the solution vector
  ierr = VecCreateMPI(comm, locDim, globDim, &yy); CHKERRQ(ierr);
  int first_dof, last_dof_plus_one;
  ierr = VecSet(yy, 0.0); CHKERRQ(ierr);//initialize 0.0

//  cout << globSite<< " " << locSite << " " << globDim << " " << locDim << endl;

  // create the bound vectors
  ierr = VecDuplicate(yy, &Lb); CHKERRQ(ierr);
  ierr = VecDuplicate(yy, &Ub); CHKERRQ(ierr);
  ierr = VecSet(Lb, 0.0); CHKERRQ(ierr); //initialize 0.0
  ierr = VecSet(Ub, 0.0); CHKERRQ(ierr); //initialize 0.0

  // create a "global" data array
  alldof = new double[globDim];

  // get local site range
  VecGetOwnershipRange(yy, &first_dof, &last_dof_plus_one);
  first_site = first_dof/numDOFPerSite;
  last_site_plus_one = last_dof_plus_one/numDOFPerSite;

  //TODO: the above calculation assumes there are NO ghost sites in the domain; may need to be 
  //      modified. Ghost sites should be included in q_M0, ...; but excluded from globDim
	//cout << (&CSembed[0])->dx << " " << CSrho[0].dx << " " << CSpair[0].dx << endl;
  // create application context
  application_context = new AppCtx(input, MM, MH, HM, HH, &q_M0, &q_H0, &sigma_M0, &sigma_H0, 
                                   &x0, &rigidAtoms, CSembed, CSrho, CSpair, 
																	 numDOFPerSite, globDim, locDim, first_site, last_site_plus_one,
                                   alldof, QPp, QWp);
  // create Tao 
  ierr = TaoCreate(PETSC_COMM_WORLD,&tao);CHKERRQ(ierr);
  ierr = TaoSetType(tao,TAOBLMVM);CHKERRQ(ierr); //TODO: let's start with quasi-Newton
  ierr = TaoSetObjectiveAndGradientRoutine(tao,Minimizer::formFunctionAndGradient,application_context);
  CHKERRQ(ierr);
  // set bounds
  ierr = TaoSetVariableBoundsRoutine(tao,Minimizer::formVariableBounds,application_context); CHKERRQ(ierr);

  ierr = TaoSetFromOptions(tao);CHKERRQ(ierr); //apply command-line options (if any)

  ierr = PetscPrintf(MPI_COMM_WORLD,"- Initialized Petsc/Tao for optimization.\n");

  return ierr;
}

PetscErrorCode Minimizer::minimizeFreeEntropy(double mubd, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0,
                          vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0,
                          vector<Vec3D> &q_M, vector<Vec3D> &q_H, vector<double> &sigma_M,
                          vector<double> &sigma_H, vector<double> &x, double &obj_fun)
{
  PetscErrorCode ierr;

  // TODO: I don't think we need to update/rebuild Tao -- but I could be wrong...


//{
//    const double x_lo = 1.0e-6;
//    const double x_hi = 1.0 - 1.0e-6;
//    const double s_lo = 1.0e-6;
//    for (size_t i=0; i<x0.size(); i++) {
//      if (x0[i] < x_lo) x0[i] = x_lo;
//      if (x0[i] > x_hi) x0[i] = x_hi;
//    }
//    for (size_t i=0; i<sigma_M0.size(); i++)
//      if (sigma_M0[i] < s_lo) sigma_M0[i] = s_lo;
//    for (size_t i=0; i<sigma_H0.size(); i++)
//      if (sigma_H0[i] < s_lo) sigma_H0[i] = s_lo;
 // }

  // set up yy 
  ierr = applicationData2Petsc(application_context, q_M0, q_H0, sigma_M0, sigma_H0, x0, yy);CHKERRQ(ierr);

//  VecView(yy,PETSC_VIEWER_STDOUT_WORLD);

  // update context
  updateApplicationContext(mubd, &q_M0, &q_H0, &sigma_M0, &sigma_H0, &x0);

  // set initial guess (yy)
  ierr = TaoSetInitialVector(tao, yy);CHKERRQ(ierr);

  // SOLVE THE OPTIMIZATION PROBLEM --> yy is updated to store the solution
  ierr = TaoSolve(tao);CHKERRQ(ierr);

  /* Get termination information */
  TaoConvergedReason reason;
  PetscReal fmin = 0.0;
  
  ierr = TaoGetSolutionStatus(tao, NULL, &fmin, NULL, NULL, NULL, &reason);
  obj_fun = (double)fmin;

  if (reason <= 0) {
    ierr = PetscPrintf(MPI_COMM_WORLD,"  WARNING: TAO did not converge!\n");CHKERRQ(ierr);
  }

  // update q_M, q_H, sigma_M, sigma_H, x 
  ierr = petscData2Application(application_context, yy, q_M, q_H, sigma_M, sigma_H, x);CHKERRQ(ierr);

  return ierr;
}


PetscErrorCode Minimizer::calculateFormationEnergy(vector<double> &x0, double gammabd, vector<double> &gamma, vector<double> &f, vector<int> &full_H)
{
    PetscErrorCode ierr = 0;
    EAM eam;

    //const double xUb = 1.0-1.0e-16; // The upper bound of hydrogen atomic fraction in the sample
    //const double xLb = application_context->input->file.xe; // The lower bound of hydrogen atomic fraction in the sample
    const double kB = 8.6173324e-5; // [eV/K]
    const double hbar = 6.582119596e-16; // [eV*s]
    const double mH = 1.6735328115271e-27; // [kg] atomic mass of H
    const double eVJ = 1.602176634e-19; // [J/eV] eV to J
    const double Am = 1.0e-10; // [m/A] Angstrom to meter
    double &T = application_context->input->file.T0; // [T]
double a_M = application_context->input->file.a_M;
//double Lx = application_context->input->file.N * a_M;
double Fxx = application_context->input->file.Fxx_current;
double Lx = Fxx * application_context->input->file.N * a_M;   // compressed
double Ly = application_context->input->file.N * a_M;
double Lz = application_context->input->file.N * a_M;



    int nM = application_context->q_M0.size();
    int nH  = application_context->q_H0.size();
    vector<vector<int> > &MM = *(application_context->MM); //create reference
    vector<vector<int> > &MH = *(application_context->MH);
    vector<vector<int> > &HM = *(application_context->HM);
    vector<vector<int> > &HH = *(application_context->HH);
    vector<vector<Vec3D> > &QPp = *(application_context->QPp);
    double &QWp = application_context->QWp;
    vector<CubicSpline> &CSembed = *(application_context->CSembed);
    vector<CubicSpline> &CSrho = *(application_context->CSrho);
    vector<CubicSpline> &CSpair = *(application_context->CSpair);

    // initiallization for electron density calculation
    for(int i=0; i<nM; i++)
        application_context->send_M[i] = 0.0;
    for(int i=0; i<nH; i++)
        application_context->send_H[i] = 0.0;

    // ------------------------------------------------------------
    // Parallelization I: calculate electron density of each site at each Gaussian point.
    // ------------------------------------------------------------
    if(application_context->first_site<nM) { // the first local site is a M
        for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
            application_context->send_M[i] = electronDensity('P', i, application_context->input->file.rc, Lx, Ly, Lz, MM, MH, HM, HH, application_context->q_M0, application_context->q_H0, application_context->sigma_M0, application_context->sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
        }
        for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
            int i = ii - nM; // index in the H vectors
            application_context->send_H[i] = electronDensity('H', i, application_context->input->file.rc, Lx, Ly, Lz, MM, MH, HM, HH, application_context->q_M0, application_context->q_H0, application_context->sigma_M0, application_context->sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
       }
    }
    else { //all the local sites are H
        for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
            int i = ii - nM; // index in the H vectors
            application_context->send_H[i] = electronDensity('H', i, application_context->input->file.rc, Lx, Ly, Lz, MM, MH, HM, HH, application_context->q_M0, application_context->q_H0, application_context->sigma_M0, application_context->sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_M, application_context->recv_M, nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    // Initialize chemical potential.
    for(int i=0; i<nH; i++)
        application_context->send_H[i] = 0.0;

    // ------------------------------------------------------------
    // Parallelization II: calculate chemical potential at H sites.
    // ------------------------------------------------------------ 
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);

    int locSize; // local size for each cores except the last one
    int lastSize; // local size for the last one
    if(nH%MPI_size == 0) {
        locSize = nH/MPI_size;
        lastSize = nH/MPI_size;
    } else {
        locSize = nH/(MPI_size-1);
        lastSize = nH%(MPI_size-1);
    }

    if((MPI_rank+1) != MPI_size) { // not the last core
        for(int ii=0; ii<locSize; ii++) {
            int i = MPI_rank*locSize + ii;
            // if(full_H[i]==1) continue; // skip the site which is fully occupied by H
            Vec3D &qi = application_context->q_H0[i]; // by XS
            double &sigmai = application_context->sigma_H0[i];
            double sum_Vx = 0.0; // the derivative of potential energy with respect to H fraction on site i

            // Gaussian quadrature for embedding energy
            sum_Vx += eam.F_H(application_context->recv_H[i], CSembed, CSrho, CSpair);

            // H-M pairs
            for(int j=0; j<HM[i].size(); j++) {
                int &jj = HM[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_M0[jj];
                double &sigmaj = application_context->sigma_M0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                     r = minimumImage(r, Lx, Ly, Lz);
                    sum_Vx += QWp*eam.F_M_deriv(application_context->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_Vx += QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);
                }
            }

            // H-H pairs
            for(int j=0; j<HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_H0[jj];
                double &sigmaj = application_context->sigma_H0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    sum_Vx += x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_Vx += x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);
                }
            }

            application_context->send_H[i] = sum_Vx;
        }
    } else { // the last core
        for(int ii=0; ii<lastSize; ii++) {
            int i = MPI_rank*locSize + ii;
            // if(full_H[i]==1) continue; // skip the site which is fully occupied by H
            Vec3D &qi = application_context->q_H0[i]; // by XS
            double &sigmai = application_context->sigma_H0[i];
            double sum_Vx = 0.0; // the derivative of potential energy with respect to H fraction on site i

            // Gaussian quadrature for embedding energy
            sum_Vx += eam.F_H(application_context->recv_H[i], CSembed, CSrho, CSpair);

            // H-M pairs
            for(int j=0; j<HM[i].size(); j++) {
                int &jj = HM[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_M0[jj];
                double &sigmaj = application_context->sigma_M0[jj];

                for(int k=0; k<QPp.size(); k++){
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    sum_Vx += QWp*eam.F_M_deriv(application_context->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_Vx += QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair);
                }
            }

            // H-H pairs
            for(int j=0; j<HH[i].size(); j++) {
                int &jj = HH[i][j]; // index of j-th neighbor of site i
                Vec3D &qj = application_context->q_H0[jj];
                double &sigmaj = application_context->sigma_H0[jj];

                for(int k=0; k<QPp.size(); k++) {
                    // Gaussian quadrature for embedding energy
                    Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                    r = minimumImage(r, Lx, Ly, Lz);
                    sum_Vx += x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H(r.norm(), CSembed, CSrho, CSpair);

                    // Gaussian quadrature for pair energy
                    sum_Vx += x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair);
                }
            }

            application_context->send_H[i] = sum_Vx;
        }
    }

    ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    CHKERRQ(ierr);
    for (int i=0; i<nH; i++) {
// if full occupancy =1, F(xi)=0, to skip full occupancy issue
 if (x0[i] >= 0.9999) {
full_H[i] = 1;
f[i] = 0.0;
 gamma[i] = gammabd;

        continue;
    }
full_H[i] = 0;
        //if(full_H[i]==1) continue; // skip the site which is fully occupied by H
        // gamma[i] = application_context->recv_H[i] + 1.5;
        f[i] = application_context->recv_H[i] + 1.5*kB*T*(log(hbar*hbar*eVJ/(kB*T*mH*application_context->sigma_H0[i]*application_context->sigma_H0[i]*Am*Am))-1.0);
        gamma[i] = f[i]/(kB*T)+log(x0[i]/(1.0-x0[i]));
        //if(full_H[i]==1) continue; // skip the site which is fully occupied by H
        //gamma[i] = min(f[i]/(kB*T)+log(x0[i]/(1.0-x0[i])), gammabd);
        /*
        if(x0[i]<xUb)
            gamma[i] = min((f[i]/(kB*T)+log(x0[i]/(1.0-x0[i]))), gammabd); // Relation between chemical potential, formation energy, and configurational entropy
        else
            gamma[i] = gammabd;
        */
    }
    return ierr;
}






PetscErrorCode Minimizer::calculateMeanPotentialEnergy(vector<double> &x0, vector<double> &V_M, vector<double> &V_H)
{
	PetscErrorCode ierr = 0;
	EAM eam;

	int nM = application_context->q_M0.size();
	int nH  = application_context->q_H0.size();
	
	// create reference
	double &rc = application_context->input->file.rc;
double a_M = application_context->input->file.a_M;
double Fxx = application_context->input->file.Fxx_current;

double Lx = Fxx*application_context->input->file.N * a_M;   // compressed
//double Lx = application_context->input->file.N * a_M;
double Ly =application_context->input->file.N * a_M;
double Lz = application_context->input->file.N * a_M;
	vector<vector<int> > &MM = *(application_context->MM);
	vector<vector<int> > &MH = *(application_context->MH);
	vector<vector<int> > &HM = *(application_context->HM);
	vector<vector<int> > &HH = *(application_context->HH);
	vector<Vec3D> &q_M0 = application_context->q_M0;
	vector<Vec3D> &q_H0 = application_context->q_H0;
	vector<double> &sigma_M0 = application_context->sigma_M0;
	vector<double> &sigma_H0 = application_context->sigma_H0;
	vector<CubicSpline> &CSembed = *(application_context->CSembed);
  vector<CubicSpline> &CSrho = *(application_context->CSrho);
	vector<CubicSpline> &CSpair = *(application_context->CSpair);	

	vector<vector<Vec3D> > &QPp = *(application_context->QPp);
	double &QWp = application_context->QWp;

	// initiallization for electron density calculation
	for(int i=0; i<nM; i++)
		application_context->send_M[i] = 0.0;
	for(int i=0; i<nH; i++)
		application_context->send_H[i] = 0.0;
	
	// ------------------------------------------------------------
	// Parallelization I: calculate electron density of each site at each Gaussian point.
	// ------------------------------------------------------------
	if(application_context->first_site<nM) { // the first local site is a M
		for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
			application_context->send_M[i] = electronDensity('P', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
		}
		for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
			int i = ii - nM; // index in the H vectors
			application_context->send_H[i] = electronDensity('H', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
		}
	}
	else { //all the local sites are H
		for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
			int i = ii - nM; // index in the H vectors
			application_context->send_H[i] = electronDensity('H', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
		}
	}	

	ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_M, application_context->recv_M, nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Initialize mean potential energy.
	for(int i=0; i<nM; i++)
		application_context->send_M[i] = 0.0;
	for(int i=0; i<nH; i++)
		application_context->send_H[i] = 0.0;

	// ------------------------------------------------------------
	// Parallelization II: calculate mean potential energy at each site.
	// ------------------------------------------------------------
	if(application_context->first_site<nM) { // the first local site is a M
		for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
		// calc. mean energy for Site i
			Vec3D &qi = q_M0[i];
			double &sigmai = sigma_M0[i];
	
			// Gaussian quadrature for embedding energy
			application_context->send_M[i] += eam.F_M(application_context->recv_M[i], CSembed, CSrho, CSpair);
			// M-M pairs
			for(int j=0; j<MM[i].size(); j++) {
				int &jj = MM[i][j]; // index of j-th neighbor of site i
				Vec3D &qj = q_M0[jj];
				double &sigmaj = sigma_M0[jj];
				
				for(int k=0; k<QPp.size(); k++) {
					Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                        r = minimumImage(r, Lx, Ly, Lz);
					// Gaussian quadrature for pair energy
					application_context->send_M[i] += 0.5*QWp*eam.phi_M(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
        }
			}

			// M-H pairs
			for(int j=0; j<MH[i].size(); j++) {
				int &jj = MH[i][j]; // index of j-th neighbor of site i
				Vec3D &qj = q_H0[jj];
				double &sigmaj = sigma_H0[jj];

				for(int k=0; k<QPp.size(); k++) {
					Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                         r = minimumImage(r, Lx, Ly, Lz);
					// Gaussian quadrature for pair energy
					application_context->send_M[i] += 0.5*x0[jj]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
        }
			}
		}

		for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
			int i = ii - nM; // index in the H vectors
			Vec3D &qi = q_H0[i]; // by XS
			double &sigmai = sigma_H0[i];

			// Gaussian quadrature for embedding energy
			application_context->send_H[i] += x0[i]*eam.F_H(application_context->recv_H[i], CSembed, CSrho, CSpair);

			// H-M pairs
			for(int j=0; j<HM[i].size(); j++) {
				int &jj = HM[i][j]; // index of j-th neighbor of site i
				Vec3D &qj = q_M0[jj];
				double &sigmaj = sigma_M0[jj];

				for(int k=0; k<QPp.size(); k++) {
					Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                         r = minimumImage(r, Lx, Ly, Lz);
					// Gaussian quadrature for pair energy
					application_context->send_H[i] += 0.5*x0[i]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
        }
			}
			// H-H pairs
			for(int j=0; j<HH[i].size(); j++) {
				int &jj = HH[i][j]; // index of j-th neighbor of site i
			  Vec3D &qj = q_H0[jj];
			  double &sigmaj = sigma_H0[jj];

				for(int k=0; k<QPp.size(); k++) {
					Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                        r = minimumImage(r, Lx, Ly, Lz);
					// Gaussian quadrature for pair energy
					application_context->send_H[i] += 0.5*x0[i]*x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
        }
			}
		}
	}
	else { //all the local sites are H
		for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
			int i = ii - nM; // index in the H vectors
			Vec3D &qi = q_H0[i]; // by XS
			double &sigmai = sigma_H0[i];

			// Gaussian quadrature for embedding energy
			application_context->send_H[i] += x0[i]*eam.F_H(application_context->recv_H[i], CSembed, CSrho, CSpair);
			
			// H-M pairs
			for(int j=0; j<HM[i].size(); j++) {
				int &jj = HM[i][j]; // index of j-th neighbor of site i
				Vec3D &qj = q_M0[jj];
				double &sigmaj = sigma_M0[jj];
			
				for(int k=0; k<QPp.size(); k++) {
			  	Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                r = minimumImage(r, Lx, Ly, Lz);
				  // Gaussian quadrature for pair energy
				  application_context->send_H[i] += 0.5*x0[i]*QWp*eam.phi_MH(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
				}
			}

			// H-H pairs
			for(int j=0; j<HH[i].size(); j++) {
				int &jj = HH[i][j]; // index of j-th neighbor of site i
			  Vec3D &qj = q_H0[jj];
			  double &sigmaj = sigma_H0[jj];
			
			  for(int k=0; k<QPp.size(); k++) {
				  Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
                                   r = minimumImage(r, Lx, Ly, Lz);
			    // Gaussian quadrature for pair energy
			    application_context->send_H[i] += 0.5*x0[i]*x0[jj]*QWp*eam.phi_H(r.norm(), CSembed, CSrho, CSpair); //Note: half of the pair energy contributes to one site
				}
			}
		}
	}

	ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_M, application_context->recv_M, nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	CHKERRQ(ierr);
	ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	CHKERRQ(ierr);

	for (int i=0; i<nM; i++) V_M[i] = application_context->recv_M[i];
	for (int i=0; i<nH; i++) V_H[i] = application_context->recv_H[i];

	return ierr;
}


PetscErrorCode Minimizer::calculateLocalStress(vector<double> &x0, vector<vector<double> > &pi_M, vector<vector<double> > &pi_H)
{
  PetscErrorCode ierr = 0;
  EAM eam;

  int nM = application_context->q_M0.size();
  int nH  = application_context->q_H0.size();
  // create reference
  double &rc = application_context->input->file.rc;
double a_M = application_context->input->file.a_M;
double Fxx = application_context->input->file.Fxx_current;
double Lx = Fxx * application_context->input->file.N * a_M;   // compressed
//double Lx = application_context->input->file.N * a_M;
double Ly = application_context->input->file.N * a_M;
double Lz = application_context->input->file.N * a_M;
  vector<vector<int> > &MM = *(application_context->MM);
  vector<vector<int> > &MH = *(application_context->MH);
  vector<vector<int> > &HM = *(application_context->HM);
  vector<vector<int> > &HH = *(application_context->HH);
  vector<Vec3D> &q_M0 = application_context->q_M0;
  vector<Vec3D> &q_H0 = application_context->q_H0;
  vector<double> &sigma_M0 = application_context->sigma_M0;
  vector<double> &sigma_H0 = application_context->sigma_H0;
  vector<CubicSpline> &CSembed = *(application_context->CSembed);
  vector<CubicSpline> &CSrho = *(application_context->CSrho);
  vector<CubicSpline> &CSpair = *(application_context->CSpair);

  const double rc_pi = rc; // cut-off distance for stress average

  const double eVJ = 1.602176634e-19; // [J/eV] eV to J
  const double Am = 1.0e-10; // [m/A] Angstrom to meter

  vector<vector<Vec3D> > &QPp = *(application_context->QPp);
  double &QWp = application_context->QWp;

  // initiallization for electron density calculation
  for(int i=0; i<nM; i++)
    application_context->send_M[i] = 0.0;
  for(int i=0; i<nH; i++)
    application_context->send_H[i] = 0.0;

  // temporary vectors for MPI
  double *send_M1, *recv_M1;
  double *send_H1, *recv_H1;
  
  send_M1 = new double[9*nM];
  recv_M1 = new double[9*nM];
  send_H1 = new double[9*nH];
  recv_H1 = new double[9*nH];

  // Initialization
  for(int i=0; i<9*nM; i++)
    send_M1[i] = 0.0;
  for(int i=0; i<9*nH; i++)
    send_H1[i] = 0.0;

  // ------------------------------------------------------------
  // Parallelization I: calculate electron density of each site at each Gaussian point.
  // ------------------------------------------------------------
  if(application_context->first_site<nM) { // the first local site is a M
    for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
      application_context->send_M[i] = electronDensity('P', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
    }
    for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
      int i = ii - nM; // index in the H vectors
      application_context->send_H[i] = electronDensity('H', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
    }
  }
  else { //all the local sites are H
    for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
      int i = ii - nM; // index in the H vectors
      application_context->send_H[i] = electronDensity('H', i, rc, Lx, Ly, Lz, MM, MH, HM, HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, QPp, QWp, CSembed, CSrho, CSpair);
    }
  }
  
  ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_M, application_context->recv_M, nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ierr = (PetscErrorCode)MPI_Allreduce(application_context->send_H, application_context->recv_H, nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // ------------------------------------------------------------
  // Parallelization II: calculate atomic virial stress at each site.
  // ------------------------------------------------------------
  if(application_context->first_site<nM) { // the first local site is a M
    for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
      // calc. force for Site i
      Vec3D &qi = application_context->q_M0[i];
      double &sigmai = application_context->sigma_M0[i];
  
      // M-M pairs
      for(int j=0; j<MM[i].size(); j++) {
        int &jj = MM[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_M0[jj];
        double &sigmaj = application_context->sigma_M0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j      

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
          r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij += QWp*eam.F_M_deriv(application_context->recv_M[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij += QWp*eam.F_M_deriv(application_context->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij += QWp*eam.phi_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.
      
        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_M1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }

      // M-H pairs
      for(int j=0; j<MH[i].size(); j++) {
        int &jj = MH[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_H0[jj];
        double &sigmaj = application_context->sigma_H0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
          r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij += x0[jj]*QWp*eam.F_M_deriv(application_context->recv_M[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij += x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij += x0[jj]*QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.
      
        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_M1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }
    }
    for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
      int i = ii - nM; // index in the H vectors
      Vec3D &qi = application_context->q_H0[i]; // by XS
      double &sigmai = application_context->sigma_H0[i];

      // H-M pairs
      for(int j=0; j<HM[i].size(); j++) {
        int &jj = HM[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_M0[jj];
        double &sigmaj = application_context->sigma_M0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
           r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij +=x0[i]*QWp*eam.F_H_deriv(application_context->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij +=QWp*eam.F_M_deriv(application_context->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij +=x0[i]* QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.
      
        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_H1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }

      // H-H pairs
      for(int j=0; j<HH[i].size(); j++) {
        int &jj = HH[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_H0[jj];
        double &sigmaj = application_context->sigma_H0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j      

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
           r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij += x0[i]*x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij += x0[i]*x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij += x0[i]*x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.
    
        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_H1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }
    }
  }
  else { //all the local sites are H
    for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
      int i = ii - nM; // index in the H vectors
      Vec3D &qi = application_context->q_H0[i]; // by XS
      double &sigmai = application_context->sigma_H0[i];

      // H-M pairs
      for(int j=0; j<HM[i].size(); j++) {
        int &jj = HM[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_M0[jj];
        double &sigmaj = application_context->sigma_M0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
          r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij +=x0[i]* QWp*eam.F_H_deriv(application_context->recv_H[i], CSembed, CSrho, CSpair)*eam.f_M_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij += x0[i]*QWp*eam.F_M_deriv(application_context->recv_M[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij +=x0[i]* QWp*eam.phi_MH_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.

        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_H1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }

      // H-H pairs
      for(int j=0; j<HH[i].size(); j++) {
        int &jj = HH[i][j]; // index of j-th neighbor of site i
        Vec3D &qj = application_context->q_H0[jj];
        double &sigmaj = application_context->sigma_H0[jj];
        Vec3D r_ij(0.0,0.0,0.0); // mean distance between i and j
        Vec3D F_ij(0.0,0.0,0.0); // force between i and j      

        for(int k=0; k<QPp.size(); k++) {
          // calculate gradient for embedding energy
          Vec3D r = qi+sqrt(2.0)*sigmai*QPp[k][0]-qj-sqrt(2.0)*sigmaj*QPp[k][1];
           r = minimumImage(r, Lx, Ly, Lz);
          r_ij += QWp*r;
          F_ij += x0[i]*x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[i], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
          F_ij += x0[i]*x0[jj]*QWp*eam.F_H_deriv(application_context->recv_H[jj], CSembed, CSrho, CSpair)*eam.f_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;

          // calculate gradient for pair energy
          F_ij += x0[i]*x0[jj]*QWp*eam.phi_H_deriv(r.norm(), CSembed, CSrho, CSpair)/r.norm()*r;
        }
        //NOTE: although the pair force to site j is the same, we MUST NOT apply 
        //the force to site j because site j might be assigned to another processor.

        for(int k=0; k<3; k++) {
          for(int l=0; l<3; l++) {
            send_H1[9*i+3*k+l] += F_ij[k]*r_ij[l];
          }
        }
      }
    }
  } 

  ierr = (PetscErrorCode)MPI_Allreduce(send_M1, recv_M1, 9*nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  CHKERRQ(ierr);
  ierr = (PetscErrorCode)MPI_Allreduce(send_H1, recv_H1, 9*nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  CHKERRQ(ierr);

  // Initialization
  for(int i=0; i<9*nM; i++)
    send_M1[i] = 0.0;
  for(int i=0; i<9*nH; i++)
    send_H1[i] = 0.0;

/*
  int tmp = 100;
  cout << "Here 0: " << recv_M1[tmp] << endl;
  cout << "Here 0 M_M: " << MM[tmp].size() << endl;
  cout << "Here 0 M_H: " << MH[tmp].size() << endl;
  cout << "Here 1: " << recv_H1[tmp] << endl;
  cout << "Here 1 H_M: " << HM[tmp].size() << endl;
  cout << "Here 1 H_H: " << HH[tmp].size() << endl;
*/

  // ------------------------------------------------------------
  // Parallelization III: calculate averaged local stresses at each site.
  // ------------------------------------------------------------
  if(application_context->first_site<nM) { // the first local site is a M
    for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
      Vec3D &qi = application_context->q_M0[i];
/*
      // loop over all M and H sites, very time-consuming
      for (int j=0; j<nM; j++) {
        Vec3D &qj = application_context->q_M0[j];
        if ((qi-qj).norm()>rc_pi) continue;
        double r_n = ((qi-qj).norm()/rc_pi);
        double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n) 
                   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
        for (int k=0; k<9; k++) { // local stress at other M site j
          send_M1[9*i+k] += recv_M1[9*j+k]*psi; 
        }
      }
      for (int j=0; j<nH; j++) {
        Vec3D &qj = application_context->q_H0[j];
        if ((qi-qj).norm()>rc_pi) continue;
        double r_n = ((qi-qj).norm()/rc_pi);
        double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n) 
                   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
        for (int k=0; k<9; k++) { // local stress at other H site j
          send_M1[9*i+k] += recv_H1[9*j+k]*psi; 
        }
      }
*/
      // loop over neighbor list, computationally efficient
      double r_n = 0.0;
      double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
                   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
      for (int k=0; k<9; k++) { // local stress at central site i
        send_M1[9*i+k] += recv_M1[9*i+k]*psi;
      }
      for (int jj=0; jj<MM[i].size(); jj++) {
        int &j = MM[i][jj];
        Vec3D &qj = application_context->q_M0[j];
	Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = dr.norm() / rc_pi;
 //       double r_n = ((qi-qj).norm()/rc_pi);
        double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
                   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
        for (int k=0; k<9; k++) { // local stress at other M site j
          send_M1[9*i+k] += recv_M1[9*j+k]*psi; 
}
}
for (int jj=0; jj<MH[i].size(); jj++) {
int &j = MH[i][jj];
Vec3D &qj = application_context->q_H0[j];
Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = (dr.norm() /rc_pi);
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
* (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
for (int k=0; k<9; k++) { // local stress at other from neighboring H site j
send_M1[9*i+k] += recv_H1[9*j+k]*psi;
}
}
}
for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // same thing, for H sites
int i = ii - nM; // index in the H vectors
Vec3D &qi = application_context->q_H0[i];
/*
// loop over all M and H sites, very time-consuming
// for (int j=0; j<nM; j++) {
// Vec3D &qj = application_context->q_M0[j];
// if ((qi-qj).norm()>rc_pi) continue;
// double r_n = ((qi-qj).norm()/rc_pi);
// double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
// * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
// for (int k=0; k<9; k++) { // local stress at other M site j
// send_H1[9*i+k] += recv_M1[9*j+k]*psi;
//  }
//   }
//   for (int j=0; j<nH; j++) {
//   Vec3D &qj = application_context->q_H0[j];
//   if ((qi-qj).norm()>rc_pi) continue;
//   double r_n = ((qi-qj).norm()/rc_pi);
//   double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
//   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
//   for (int k=0; k<9; k++) { // local stress at other H site j
//   send_H1[9*i+k] += recv_H1[9*j+k]*psi;
//    }
//    }
//    */
//    // loop over neighbor list, computationally efficient
double r_n = 0.0;
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
* (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
for (int k=0; k<9; k++) { // local stress at central site i
send_H1[9*i+k] += recv_H1[9*i+k]*psi;
}
for (int jj=0; jj<HM[i].size(); jj++) {
int &j = HM[i][jj];
Vec3D &qj = application_context->q_M0[j];
Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = (dr.norm()/rc_pi);
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
                   * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
for (int k=0; k<9; k++) { // local stress at other M site j
send_H1[9*i+k] += recv_M1[9*j+k]*psi;
}
}
for (int jj=0; jj<HH[i].size(); jj++) {
int &j = HH[i][jj];
Vec3D &qj = application_context->q_H0[j];
Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = (dr.norm()/rc_pi);
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
* (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
for (int k=0; k<9; k++) { // local stress at other H site j
send_H1[9*i+k] += recv_H1[9*j+k]*psi;
}
}
}
}
else { //all the local sites are H
for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
int i = ii - nM; // index in the H vectors
Vec3D &qi = application_context->q_H0[i];
// loop over neighbor list, computationally efficient
double r_n = 0.0;
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
* (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
 for (int k=0; k<9; k++) { // local stress at central site i
send_H1[9*i+k] += recv_H1[9*i+k]*psi;
}
for (int jj=0; jj<HM[i].size(); jj++) {
int &j = HM[i][jj];
Vec3D &qj = application_context->q_M0[j];
Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = (dr.norm()/rc_pi);
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
 * (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
 for (int k=0; k<9; k++) { // local stress at other M site j
send_H1[9*i+k] += recv_M1[9*j+k]*psi;
}
}
for (int jj=0; jj<HH[i].size(); jj++) {
int &j = HH[i][jj];
Vec3D &qj = application_context->q_H0[j];
Vec3D dr = qi - qj;
dr = minimumImage(dr, Lx, Ly, Lz);
double r_n = (dr.norm()/rc_pi);
double psi = 2.3873241 / (rc_pi*rc_pi*rc_pi) * (1.0-r_n*r_n) * (1.0-r_n*r_n)
* (0.5-1.5*(r_n-0.5)+2.0*(r_n-0.5)*(r_n-0.5)*(r_n-0.5)); // normalized localization function: 7th order polynomial
for (int k=0; k<9; k++) { // local stress at other H site j
send_H1[9*i+k] += recv_H1[9*j+k]*psi;
}
}
}
}
ierr = (PetscErrorCode)MPI_Allreduce(send_M1, recv_M1, 9*nM, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
CHKERRQ(ierr);
ierr = (PetscErrorCode)MPI_Allreduce(send_H1, recv_H1, 9*nH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 CHKERRQ(ierr);
double omega = 4.0/3.0*3.1415926*rc_pi*rc_pi*rc_pi; // total volume
for (int i=0; i<nM; i++) { // averaged local stress at M site i
for (int k=0; k<9; k++) { // local stress at central site i
pi_M[i][k] = recv_M1[9*i+k]*0.5*eVJ/(Am*Am*Am)*1.0e-9; // [GPa]
}
}
for (int i=0; i<nH; i++) { // averaged local stress at H site i
for (int k=0; k<9; k++) { // local stress at central site i
pi_H[i][k] = recv_H1[9*i+k]*0.5*eVJ/(Am*Am*Am)*1.0e-9; // [GPa]
}
}
delete[] send_M1;
  delete[] recv_M1;
  delete[] send_H1;
  delete[] recv_H1;
  return ierr;
}
PetscErrorCode Minimizer::formVariableBounds(Tao tao, Vec Lb, Vec Ub, void* ptr)
{
PetscErrorCode ierr = 0;
 AppCtx *app = (AppCtx *) ptr;
vector<Vec3D> Ub_q_M(app->q_M0), Lb_q_M(app->q_M0), Ub_q_H(app->q_H0), Lb_q_H(app->q_H0);
    vector<double> Ub_sigma_M(app->sigma_M0), Lb_sigma_M(app->sigma_M0);
    vector<double> Ub_sigma_H(app->sigma_H0), Lb_sigma_H(app->sigma_H0);
    vector<double> Ub_x_H(app->x0), Lb_x_H(app->x0);

//for (int i=0; i<Ub_q_M.size(); i++) {
//Lb_q_M[i] = Vec3D(PETSC_NINFINITY, PETSC_NINFINITY, PETSC_NINFINITY);
// Ub_q_M[i] = Vec3D(PETSC_INFINITY,  PETSC_INFINITY,  PETSC_INFINITY);
//Lb_sigma_M[i] = 1.0e-6;
//Ub_sigma_M[i] = 0.5;
//} 
// double  maxDisp = 0.5; // Angstrom, temporary safety bound
for (int i=0; i<Ub_q_M.size(); i++) {
        Lb_q_M[i] = app->q_M0[i];
        Ub_q_M[i] = app->q_M0[i];
    //    Lb_q_M[i][0] = app->q_M0[i][0] - maxDisp;
    //    Ub_q_M[i][0] = app->q_M0[i][0] + maxDisp;

      //  Lb_q_M[i][1] = app->q_M0[i][1] - maxDisp;
    //   Ub_q_M[i][1] = app->q_M0[i][1] + maxDisp;

    //    Lb_q_M[i][2] = app->q_M0[i][2] - maxDisp;
    //    Ub_q_M[i][2] = app->q_M0[i][2] + maxDisp;

        Lb_sigma_M[i] = 1.0e-6;
        Ub_sigma_M[i] = PETSC_INFINITY;
}
//for (int i = 0; i < Ub_q_H.size(); i++) {
//
//    Lb_q_H[i][0] = app->q_H0[i][0] - maxDisp;
//    Ub_q_H[i][0] = app->q_H0[i][0] + maxDisp;

//    Lb_q_H[i][1] = app->q_H0[i][1] - maxDisp;
//    Ub_q_H[i][1] = app->q_H0[i][1] + maxDisp;

//    Lb_q_H[i][2] = app->q_H0[i][2] - maxDisp;
//    Ub_q_H[i][2] = app->q_H0[i][2] + maxDisp;

//    Lb_sigma_H[i] = 1.0e-6;
//    Ub_sigma_H[i] = PETSC_INFINITY;

//    Lb_x_H[i] = 1.0e-16;
//    Ub_x_H[i] = 1.0 - 1.0e-16;
//}
for (int i=0; i<Ub_q_H.size(); i++) {
        Lb_q_H[i] = app->q_H0[i];
        Ub_q_H[i] = app->q_H0[i];
        //Lb_q_H[i] = Vec3D(PETSC_NINFINITY, PETSC_NINFINITY, PETSC_NINFINITY);
        //Ub_q_H[i] = Vec3D(PETSC_INFINITY,  PETSC_INFINITY,  PETSC_INFINITY);
        Lb_sigma_H[i] = 1.0e-6;
        Ub_sigma_H[i] = 0.5;
        Lb_x_H[i] = 1.0e-16;
        Ub_x_H[i] = 1.0 - 1.0e-16;
}
//for (int i=0; i<Ub_q_M.size(); i++) {
 //       Lb_q_M[i] = app->q_M0[i];      // freeze at deformed position
//        Ub_q_M[i] = app->q_M0[i];
 
//        Lb_sigma_M[i] = 1.0e-6;
//        Ub_sigma_M[i] = PETSC_INFINITY;           // was PETSC_NINFINITY — that's a bug (see below)
//}
//for (int i=0; i<Ub_q_H.size(); i++) {
//        Lb_q_H[i] = app->q_H0[i];      // freeze at deformed position
//        Ub_q_H[i] = app->q_H0[i];
//
//        Lb_sigma_H[i] = 1.0e-6;
//        Ub_sigma_H[i] = PETSC_INFINITY;
//        Lb_x_H[i] = 1.0e-16;            // raise from 1e-16 to avoid log(0)
//        Ub_x_H[i] = 1.0 - 1.0e-16;
//}

//if (app->rigidAtoms.size() != 3) {
//    PetscPrintf(MPI_COMM_WORLD,
//                "ERROR: rigidAtoms must contain exactly 3 atoms.\n");
//    return PETSC_ERR_ARG_WRONG;
//}

//int atomA = app->rigidAtoms[0];
//int atomB = app->rigidAtoms[1];
//int atomC = app->rigidAtoms[2];

//int nRigidM = app->q_M0.size();

//if (atomA < 0 || atomA >= nRigidM ||
//    atomB < 0 || atomB >= nRigidM ||
//    atomC < 0 || atomC >= nRigidM) {
//    PetscPrintf(MPI_COMM_WORLD,
//                "ERROR: invalid rigid atom indices: A=%d B=%d C=%d\n",
//                atomA, atomB, atomC);
//    return PETSC_ERR_ARG_OUTOFRANGE;

/*if (app->rigidAtoms.size() != 4) {
    PetscPrintf(MPI_COMM_WORLD,
                "ERROR: rigidAtoms must contain exactly 3 atoms.\n");
    return PETSC_ERR_ARG_WRONG;
}
int atomA = app->rigidAtoms[0];
int atomB = app->rigidAtoms[1];
int atomC = app->rigidAtoms[2];
int atomD = app->rigidAtoms[3];
int nRigidM = static_cast<int>(app->q_M0.size());
//int nRigidM = app->q_M0.size();
//if (atomA < 0 || atomA >= nRigidM) {
//    PetscPrintf(MPI_COMM_WORLD,
//                "ERROR: invalid rigid atom A=%d\n", atomA);
//    return PETSC_ERR_ARG_OUTOFRANGE;
//}

if (atomA < 0 || atomA >= nRigidM ||
    atomB < 0 || atomB >= nRigidM ||
    atomC < 0 || atomC >= nRigidM ||
    atomD < 0 || atomD >= nRigidM) {
    PetscPrintf(MPI_COMM_WORLD,
                "ERROR: invalid rigid atom indices: A=%d B=%d C=%d D=%d\n",
                atomA, atomB, atomC, atomD);
    return PETSC_ERR_ARG_OUTOFRANGE;
}

//if (atomA < 0 || atomA >= nRigidM ||
//    atomB < 0 || atomB >= nRigidM ||
//    atomC < 0 || atomC >= nRigidM) {
//    PetscPrintf(MPI_COMM_WORLD,
//                "ERROR: invalid rigid atom indices: A=%d B=%d C=%d\n",
//                atomA, atomB, atomC);
//    return PETSC_ERR_ARG_OUTOFRANGE;

//}

Lb_q_M[atomA] = app->q_M0[atomA];
Ub_q_M[atomA] = app->q_M0[atomA];
//Atom B: pin y,z  -- only if B is active (>=0). Under PBC B=-1 and is skipped
if (atomB >= 0 && atomB < nRigidM) {
    Lb_q_M[atomB][1] = app->q_M0[atomB][1];
    Ub_q_M[atomB][1] = app->q_M0[atomB][1];
    Lb_q_M[atomB][2] = app->q_M0[atomB][2];
    Ub_q_M[atomB][2] = app->q_M0[atomB][2];
}
// Atom C: pin x,z  -- only if active
if (atomC >= 0 && atomC < nRigidM) {
    Lb_q_M[atomC][0] = app->q_M0[atomC][0];
    Ub_q_M[atomC][0] = app->q_M0[atomC][0];
    Lb_q_M[atomC][2] = app->q_M0[atomC][2];
    Ub_q_M[atomC][2] = app->q_M0[atomC][2];
}
// Atom D: pin x,y  -- only if active
if (atomD >= 0 && atomD < nRigidM) {
    Lb_q_M[atomD][0] = app->q_M0[atomD][0];
    Ub_q_M[atomD][0] = app->q_M0[atomD][0];
    Lb_q_M[atomD][1] = app->q_M0[atomD][1];
    Ub_q_M[atomD][1] = app->q_M0[atomD][1];
}

*/
//Lb_q_M[atomB][1] = app->q_M0[atomB][1];
//Ub_q_M[atomB][1] = app->q_M0[atomB][1];

//Lb_q_M[atomB][2] = app->q_M0[atomB][2];
//Ub_q_M[atomB][2] = app->q_M0[atomB][2];
//Lb_q_M[atomC][0] = app->q_M0[atomC][0];
//Ub_q_M[atomC][0] = app->q_M0[atomC][0];

//Lb_q_M[atomC][2] = app->q_M0[atomC][2];
//Ub_q_M[atomC][2] = app->q_M0[atomC][2];
//Lb_q_M[atomD][0] = app->q_M0[atomD][0];
//Ub_q_M[atomD][0] = app->q_M0[atomD][0];

//Lb_q_M[atomD][1] = app->q_M0[atomD][1];
//Ub_q_M[atomD][1] = app->q_M0[atomD][1];
//
//
//
//Lb_q_M[atomA] = app->q_M0[atomA];
//Ub_q_M[atomA] = app->q_M0[atomA];
//Lb_q_M[atomB][1] = app->q_M0[atomB][1];
//Ub_q_M[atomB][1] = app->q_M0[atomB][1];

//Lb_q_M[atomB][2] = app->q_M0[atomB][2];
//Ub_q_M[atomB][2] = app->q_M0[atomB][2];
//Lb_q_M[atomC][0] = app->q_M0[atomC][0];
//Ub_q_M[atomC][0] = app->q_M0[atomC][0];

//Lb_q_M[atomC][2] = app->q_M0[atomC][2];
//Ub_q_M[atomC][2] = app->q_M0[atomC][2];
//PetscPrintf(MPI_COMM_WORLD,
//            "Rigid constraints applied: A=%d fixed xyz, B=%d, C=%d, D=%d (>=0 means pinned)\n",
//            atomA, atomB, atomC, atomD);
//PetscPrintf(MPI_COMM_WORLD,
//            "Rigid constraints applied: A=%d fixed xyz, B=%d fixed yz, C=%d fixed xz, D=%d fixed xy\n",
//            atomA, atomB, atomC, atomD);
//PetscPrintf(MPI_COMM_WORLD,
//            "Rigid constraints applied: A=%d fixed xyz, B=%d fixed yz, C=%d fixed z\n",
//            atomA, atomB, atomC);

      ierr = applicationData2Petsc(app, Ub_q_M, Ub_q_H, Ub_sigma_M, Ub_sigma_H, Ub_x_H, Ub); CHKERRQ(ierr);
             ierr = applicationData2Petsc(app, Lb_q_M, Lb_q_H, Lb_sigma_M, Lb_sigma_H, Lb_x_H, Lb); CHKERRQ(ierr);

                 return ierr;
   }
//if (nM < 3) {
//    PetscPrintf(PETSC_COMM_WORLD,
//        "ERROR: Need at least 3 metal atoms to fix rigid body motion.\n");
//    MPI_Abort(PETSC_COMM_WORLD, 1);
//}
//Choose Atom A
//int atomA = 0;
//double min_corner_value = app->q_M0[0][0] + app->q_M0[0][1] + app->q_M0[0][2];

//for (int i = 1; i < nM; i++) {
//    double value = app->q_M0[i][0] + app->q_M0[i][1] + app->q_M0[i][2];

//    if (value < min_corner_value) {
//        min_corner_value = value;
//        atomA = i;
//    }
//}
// 2. Choose Atom B
//int atomB = -1;
//double max_dx = -1.0;

//for (int i = 0; i < nM; i++) {
//    if (i == atomA) continue;

//    double dx = fabs(app->q_M0[i][0] - app->q_M0[atomA][0]);

  //  if (dx > max_dx) {
//        max_dx = dx;
//        atomB = i;
//    }
//}
// 3. Choose Atom C
//int atomC = -1;
//double max_dist_from_line = -1.0;

//Vec3D A = app->q_M0[atomA];
//Vec3D B = app->q_M0[atomB];
//Vec3D AB = B - A;

//double AB_norm = AB.norm();

//if (AB_norm < 1.0e-12) {
//    PetscPrintf(PETSC_COMM_WORLD,
//        "ERROR: Atom A and Atom B are too close.\n");
//    MPI_Abort(PETSC_COMM_WORLD, 1);
//}

//for (int i = 0; i < nM; i++) {
//    if (i == atomA || i == atomB) continue;

//    Vec3D C = app->q_M0[i];
//    Vec3D AC = C - A;
// distance from point C to line AB
 
//double cx = AC[1]*AB[2] - AC[2]*AB[1];
//double cy = AC[2]*AB[0] - AC[0]*AB[2];
//double cz = AC[0]*AB[1] - AC[1]*AB[0];

//double cross_norm = sqrt(cx*cx + cy*cy + cz*cz);
//double dist_from_line = cross_norm / AB_norm;

//    if (dist_from_line > max_dist_from_line) {
//        max_dist_from_line = dist_from_line;
//        atomC = i;
//    }
//}

//if (atomC < 0 || max_dist_from_line < 1.0e-8) {
//    PetscPrintf(PETSC_COMM_WORLD,
//        "ERROR: Could not find a good third atom for rigid-body constraint.\n");
//    MPI_Abort(PETSC_COMM_WORLD, 1);
//}

//PetscPrintf(PETSC_COMM_WORLD,
//    "Rigid body fixed atoms: A = %d, B = %d, C = %d\n",
//    atomA, atomB, atomC);


// Atom A: fix x, y, z
//Lb_q_M[atomA][0] = app->q_M0[atomA][0];
//Ub_q_M[atomA][0] = app->q_M0[atomA][0];

//Lb_q_M[atomA][1] = app->q_M0[atomA][1];
//Ub_q_M[atomA][1] = app->q_M0[atomA][1];

//Lb_q_M[atomA][2] = app->q_M0[atomA][2];
//Ub_q_M[atomA][2] = app->q_M0[atomA][2];

// Atom B: fix y, z
 //Lb_q_M[atomB][1] = app->q_M0[atomB][1];
 //Ub_q_M[atomB][1] = app->q_M0[atomB][1];

 //Lb_q_M[atomB][2] = app->q_M0[atomB][2];
 //Ub_q_M[atomB][2] = app->q_M0[atomB][2];

 // Atom C: fix z
// Lb_q_M[atomC][2] = app->q_M0[atomC][2];
// Ub_q_M[atomC][2] = app->q_M0[atomC][2];


// Atom B: fix y, z
// Lb_q_M[atomB][1] = app->q_M0[atomB][1];
// Ub_q_M[atomB][1] = app->q_M0[atomB][1];

// Lb_q_M[atomB][2] = app->q_M0[atomB][2];
// Ub_q_M[atomB][2] = app->q_M0[atomB][2];

 // Atom C: fix z
// Lb_q_M[atomC][2] = app->q_M0[atomC][2];
// Ub_q_M[atomC][2] = app->q_M0[atomC][2];

 //int atomA = 3924;
 //int atomB = 3990;
 //int atomC = 3984;
 
 //Atom A: fix x, y, z
//  Lb_q_M[atomA][0] = app->q_M0[atomA][0];
//  Ub_q_M[atomA][0] = app->q_M0[atomA][0];
 
 // Lb_q_M[atomA][1] = app->q_M0[atomA][1];
 // Ub_q_M[atomA][1] = app->q_M0[atomA][1];
 
//  Lb_q_M[atomA][2] = app->q_M0[atomA][2];
//  Ub_q_M[atomA][2] = app->q_M0[atomA][2];
 
 // Atom B: fix y, z
// Lb_q_M[atomB][1] = app->q_M0[atomB][1];
// Ub_q_M[atomB][1] = app->q_M0[atomB][1];
 
// Lb_q_M[atomB][2] = app->q_M0[atomB][2];
// Ub_q_M[atomB][2] = app->q_M0[atomB][2];
 
 // Atom C: fix z
 //Lb_q_M[atomC][2] = app->q_M0[atomC][2];
 //Ub_q_M[atomC][2] = app->q_M0[atomC][2];



 




//    for (int i=0; i<Ub_q_M.size(); i++) {
//        
//        Lb_q_M[i] = Vec3D(PETSC_NINFINITY,PETSC_NINFINITY,PETSC_NINFINITY);
//        Ub_q_M[i] = Vec3D(PETSC_INFINITY,PETSC_INFINITY,PETSC_INFINITY);
//       Lb_q_M[i] = Vec3D(); 
        // for minimizing only with respect to sigma
        //Lb_q_M[i] = app->q_M0[i];
        //Ub_q_M[i] = app->q_M0[i];
        
//        Lb_sigma_M[i] = 1.0e-6;
//        Ub_sigma_M[i] = PETSC_INFINITY;
//    }

//const double q_window = 0.5;  // +/- 0.5 Angstrom leash around reference site
//
//    for (int i=0; i<Ub_q_M.size(); i++) {
//     Lb_q_M[i] = app->q_M0[i] - Vec3D(q_window, q_window, q_window);
//        Ub_q_M[i] = app->q_M0[i] + Vec3D(q_window, q_window, q_window);

//        Lb_sigma_M[i] = 1.0e-6;
//        Ub_sigma_M[i] = PETSC_INFINITY;
//    }

 //   for (int i=0; i<Ub_q_H.size(); i++) {
//
//        Lb_q_H[i] = app->q_H0[i] - Vec3D(q_window, q_window, q_window);
//        Ub_q_H[i] = app->q_H0[i] + Vec3D(q_window, q_window, q_window);
//
//        Lb_sigma_H[i] = 1.0e-6;
//        Ub_sigma_H[i] = PETSC_INFINITY;
//
//        Lb_x_H[i] = 1.0e-6;
 //       Ub_x_H[i] = 1.0-1.0e-6;
//    }

// Rigid-body constraint
// 3 translation + 3 rotation zero modes.
// freezing that atom in place.
// {
//    int nM = app->q_M0.size();
//    const int A = 0;
//    const int B = nM / 2;
//    const int C = nM - 1;
//if (app->first_site == 0) {
//        printf("nM = %d\n", nM);
//        printf("Anchor A (%d): %g %g %g\n", A, app->q_M0[A][0], app->q_M0[A][1], app->q_M0[A][2])
PetscErrorCode Minimizer::findLocalNeighbors(vector<vector<int> > &MM_ext, vector<vector<int> > &MH_ext,
                  vector<vector<int> > &HM_ext, vector<vector<int> > &HH_ext,
                  vector<vector<int> > &MM, vector<vector<int> > &MH,
                  vector<vector<int> > &HM, vector<vector<int> > &HH,
                  int &maxNeib)

{
    PetscErrorCode ierr = 0;

    int maxNeibTmp = 500;    
    int nM = application_context->q_M0.size();
    int nH  = application_context->q_H0.size();

    vector<Vec3D> &q_M0 = application_context->q_M0;
    vector<Vec3D> &q_H0 = application_context->q_H0;
    double rc = application_context->input->file.rc;

		//double rc_diff = application_context->input->file.rc_diff;
  
  double a_M = application_context->input->file.a_M;
double Fxx =
    application_context->input->file.Fxx_current; 
double Lx = Fxx * application_context->input->file.N * a_M;   // compressed
//double Lx = application_context->input->file.N * a_M;
double Ly = application_context->input->file.N * a_M;
double Lz = application_context->input->file.N * a_M;

    // temporary vectors for MPI
    int *send_M1, *recv_M1;
    int *send_M2, *recv_M2;
    int *send_H1, *recv_H1;
    int *send_H2, *recv_H2;
    int send_MH = 0;

    send_M1 = new int[nM*maxNeibTmp];
    recv_M1 = new int[nM*maxNeibTmp];
    send_M2 = new int[nM*maxNeibTmp];
    recv_M2 = new int[nM*maxNeibTmp];
    send_H1 = new int[nH*maxNeibTmp];
    recv_H1 = new int[nH*maxNeibTmp];
    send_H2 = new int[nH*maxNeibTmp];
    recv_H2 = new int[nH*maxNeibTmp];


    // initiallization
    for(int i=0; i<nM*maxNeibTmp; i++) {
        send_M1[i] = 0;
        send_M2[i] = 0;
    }

    for(int i=0; i<nH*maxNeibTmp; i++) {
        send_H1[i] = 0;
        send_H2[i] = 0;
    }

    // Find and store neighbors
    maxNeib = 0;
    int maxNei = 500; //max number of neighbors for each atom. 
    double qLocal[3];

    // ------------------------------------------------------------
    // Parallelization 
    // ------------------------------------------------------------
    if(application_context->first_site<nM) { // the first local site is a M
        for(int i=application_context->first_site; i<min(nM, application_context->last_site_plus_one); i++) {
        Vec3D qi = q_M0[i]; 
            // Store local sites in a KD-Tree with K = 3.
           PointIn3D *atoms = new PointIn3D[MM_ext[i].size()+MH_ext[i].size()];
            for(int j=0; j<MM_ext[i].size(); j++){
               Vec3D pj = nearestImage(q_M0[MM_ext[i][j]], qi, Lx, Ly, Lz);
                atoms[j] = PointIn3D((int)j, pj);
}
     
            for(int j=0; j<MH_ext[i].size(); j++){
Vec3D pj = nearestImage(q_H0[MH_ext[i][j]], qi, Lx, Ly, Lz);
                atoms[(int)(MM_ext[i].size() + j)] = PointIn3D((int)(MM_ext[i].size() + j), pj);
 }
            KDTree<PointIn3D> atomTree(MM_ext[i].size()+MH_ext[i].size(), atoms); //Note: atoms are re-ordered
            
            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++)
                qLocal[j] = q_M0[i][j]; // position of i
          // find neighbors of i within rc --> populates "candidates"
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
          // debug
            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                    << ". Increase maxNei! " << endl;
                exit(-1);
            }

            int iM1 = 0;
            int iM2 = 0;
          // populates MM, MH  
  for(int j=0; j<nFound; j++) {

    int nei = candidates[j].pid(); // index inside combined MM_ext + MH_ext pool

    Vec3D dr(candidates[j].val(0)-qLocal[0],
             candidates[j].val(1)-qLocal[1],
             candidates[j].val(2)-qLocal[2]);

    if(dr.norm()*dr.norm() > rc*rc)
        continue;

    if(nei < MM_ext[i].size()) { // this neighbor is M

        int realM = MM_ext[i][nei];

        bool duplicate = false;
        for(int k=0; k<iM1; k++) {
            if(send_M1[i*maxNeibTmp+k] == realM + 1) {
                duplicate = true;
                break;
            }
        }

        if(!duplicate) {
            if(iM1 >= maxNeibTmp) {
                cerr << "ERROR: M-M local neighbors exceed maxNeibTmp for M atom "
                     << i << endl;
                exit(-1);
            }

            send_M1[i*maxNeibTmp+iM1] = realM + 1;
            iM1++;
            nNeib++;
        }

    } else { // this neighbor is H

        int realH = MH_ext[i][nei - MM_ext[i].size()];

        bool duplicate = false;
        for(int k=0; k<iM2; k++) {
            if(send_M2[i*maxNeibTmp+k] == realH + 1) {
                duplicate = true;
                break;
            }
        }

        if(!duplicate) {
            if(iM2 >= maxNeibTmp) {
                cerr << "ERROR: M-H local neighbors exceed maxNeibTmp for M atom "
                     << i << endl;
                exit(-1);
            }

            send_M2[i*maxNeibTmp+iM2] = realH + 1;
            iM2++;
            nNeib++;
        }
    }
}
           if(nNeib>send_MH) send_MH = nNeib;
            delete[] atoms;
        }
     for(int ii=nM; ii<application_context->last_site_plus_one; ii++) { // H sites
            int i = ii - nM;
            Vec3D qi = q_H0[i];   // search center (this H site)
            PointIn3D *atoms = new PointIn3D[HM_ext[i].size()+HH_ext[i].size()];
            for(int j=0; j<HM_ext[i].size(); j++)
{
Vec3D pj = nearestImage(q_M0[HM_ext[i][j]], qi, Lx, Ly, Lz);
atoms[j] = PointIn3D(j, pj);
}
            for(int j=0; j<HH_ext[i].size(); j++){
Vec3D pj = nearestImage(q_H0[HH_ext[i][j]], qi, Lx, Ly, Lz);
                atoms[(int)(HM_ext[i].size() + j)] = PointIn3D((int)(HM_ext[i].size() + j), pj);
            }
            KDTree<PointIn3D> atomTree(HM_ext[i].size()+HH_ext[i].size(), atoms);

            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++)
 qLocal[j] = qi[j];
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);

            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                     << ". Increase maxNei! " << endl;
                exit(-1);
            }
int iH1 = 0, iH2 = 0;
            for(int j=0; j<nFound; j++) {
                int nei = candidates[j].pid();
                Vec3D dr(candidates[j].val(0)-qLocal[0],
                         candidates[j].val(1)-qLocal[1],
                         candidates[j].val(2)-qLocal[2]);
              //  dr = minimumImage(dr, Lx, Ly, Lz);
                if(dr.norm()*dr.norm() > rc*rc) continue;
if(nei < HM_ext[i].size()) { // this neighbor is M

    int realM = HM_ext[i][nei];

    bool duplicate = false;
    for(int k=0; k<iH2; k++) {
        if(send_H2[i*maxNeibTmp+k] == realM + 1) {
            duplicate = true;
            break;
        }
    }

    if(!duplicate) {
        if(iH2 >= maxNeibTmp) {
            cerr << "ERROR: H-M local neighbors exceed maxNeibTmp for H atom "
                 << i << endl;
            exit(-1);
        }

        send_H2[i*maxNeibTmp+iH2] = realM + 1;
        iH2++;
        nNeib++;
    }

} else { // this neighbor is H

    int realH = HH_ext[i][nei - HM_ext[i].size()];

    bool duplicate = false;
    for(int k=0; k<iH1; k++) {
        if(send_H1[i*maxNeibTmp+k] == realH + 1) {
            duplicate = true;
            break;
        }
    }

    if(!duplicate) {
        if(iH1 >= maxNeibTmp) {
            cerr << "ERROR: H-H local neighbors exceed maxNeibTmp for H atom "
                 << i << endl;
            exit(-1);
        }

        send_H1[i*maxNeibTmp+iH1] = realH + 1;
        iH1++;
        nNeib++;
    }

}
}
            if(nNeib>send_MH) send_MH = nNeib;
            delete[] atoms;
        }
    }
    else { // all local sites are H
        for(int ii=application_context->first_site; ii<application_context->last_site_plus_one; ii++) {
            int i = ii - nM;
            Vec3D qi = q_H0[i];   // search center (this H site)
            PointIn3D *atoms = new PointIn3D[HM_ext[i].size()+HH_ext[i].size()];
            for(int j=0; j<HM_ext[i].size(); j++)
{
Vec3D pj = nearestImage(q_M0[HM_ext[i][j]], qi, Lx, Ly, Lz);
atoms[j] = PointIn3D(j, pj);
}
for(int j=0; j<HH_ext[i].size(); j++){
                Vec3D pj = nearestImage(q_H0[HH_ext[i][j]], qi, Lx, Ly, Lz);
                atoms[(int)(HM_ext[i].size() + j)] = PointIn3D((int)(HM_ext[i].size() + j), pj);
            }
            KDTree<PointIn3D> atomTree(HM_ext[i].size()+HH_ext[i].size(), atoms);
            int nNeib = 0;
            PointIn3D candidates[maxNei];
            for(int j=0; j<3; j++) qLocal[j] = qi[j];
            int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
            if(nFound>maxNei) {
                cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei
                     << ". Increase maxNei! " << endl;
                exit(-1);
            }

            int iH1 = 0, iH2 = 0;
            for(int j=0; j<nFound; j++) {
                int nei = candidates[j].pid();
                Vec3D dr(candidates[j].val(0)-qLocal[0],
                         candidates[j].val(1)-qLocal[1],
                         candidates[j].val(2)-qLocal[2]);
                //dr = minimumImage(dr, Lx, Ly, Lz);
                if(dr.norm()*dr.norm() > rc*rc) continue;
   if(nei < HM_ext[i].size()) { // this neighbor is M

    int realM = HM_ext[i][nei];

    bool duplicate = false;
    for(int k=0; k<iH2; k++) {
        if(send_H2[i*maxNeibTmp+k] == realM + 1) {
            duplicate = true;
            break;
        }
    }

    if(!duplicate) {
        if(iH2 >= maxNeibTmp) {
            cerr << "ERROR: H-M local neighbors exceed maxNeibTmp for H atom "
                 << i << endl;
            exit(-1);
        }

        send_H2[i*maxNeibTmp+iH2] = realM + 1;
        iH2++;
        nNeib++;
    }

} else { // this neighbor is H

    int realH = HH_ext[i][nei - HM_ext[i].size()];

    bool duplicate = false;
    for(int k=0; k<iH1; k++) {
        if(send_H1[i*maxNeibTmp+k] == realH + 1) {
            duplicate = true;
            break;
        }
    }

    if(!duplicate) {
        if(iH1 >= maxNeibTmp) {
            cerr << "ERROR: H-H local neighbors exceed maxNeibTmp for H atom "
                 << i << endl;
            exit(-1);
        }

        send_H1[i*maxNeibTmp+iH1] = realH + 1;
        iH1++;
        nNeib++;
    }
}
}
            if(nNeib>send_MH) send_MH = nNeib;
            delete[] atoms;
        }
    }

    
    ierr = (PetscErrorCode)MPI_Allreduce(send_M1, recv_M1, nM*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_M2, recv_M2, nM*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_H1, recv_H1, nH*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    ierr = (PetscErrorCode)MPI_Allreduce(send_H2, recv_H2, nH*maxNeibTmp, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    ierr = (PetscErrorCode)MPI_Allreduce(&send_MH, &maxNeib, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    delete[] send_M1;
    delete[] send_M2;
    delete[] send_H1;
    delete[] send_H2;

    if(!(MM.empty() && MH.empty() && HM.empty() && HH.empty())) {
        MM.clear();
        MH.clear();
        HM.clear();
        HH.clear();
    }

    // Initialize neighbor lists
    vector<int> tmp; //an empty vector
    tmp.reserve(maxNeibTmp);
    MM.assign(nM, tmp); //M neighbors of M
    MH.assign(nM, tmp); //H neighbors of M
    HM.assign(nH, tmp); //M neighbors of H
    HH.assign(nH, tmp); //H neighbors of H

    vector<int> tmp1; //a zero vector
    tmp1.assign(3, 0);


    for(int i=0; i<nM; i++) {
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_M1[i*maxNeibTmp+j]==0) break;
            MM[i].push_back(recv_M1[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_M2[i*maxNeibTmp+j]==0) break;
            MH[i].push_back(recv_M2[i*maxNeibTmp+j]-1);
        }
    }

    for(int i=0; i<nH; i++) {
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_H1[i*maxNeibTmp+j]==0) break;
            HH[i].push_back(recv_H1[i*maxNeibTmp+j]-1);
        }
        for (int j=0; j<maxNeibTmp; j++) {
            if(recv_H2[i*maxNeibTmp+j]==0) break;
            HM[i].push_back(recv_H2[i*maxNeibTmp+j]-1);
        }
    }

    delete[] recv_M1;
    delete[] recv_M2;
    delete[] recv_H1;
    delete[] recv_H2;
int cmin=1e9, cmax=0;
for(int i=0;i<nM;i++){ int c=MM[i].size(); if(c<cmin)cmin=c; if(c>cmax)cmax=c; }
PetscPrintf(MPI_COMM_WORLD, "M-M coordination: min=%d max=%d\n", cmin, cmax);
    return ierr;
}

