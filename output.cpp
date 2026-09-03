#include <time.h>
#include "output.h"
#include "input.h"
#include "Vector3D.h"
#include <cmath>
using namespace std;
//static inline double wrapCoord(double x, double lo, double L)
//{
//    while(x < lo) x += L;
//    while(x >= lo + L) x -= L;
//    return x;
//}
static inline double wrapCoord(double x, double lo, double L)
{
    const double tol = 1.0e-8;

    double wrapped = x - L * floor((x - lo) / L);
if (wrapped >= lo + L - tol) {
        wrapped = lo;
    }

    // Also snap values numerically close to the lower boundary.
 if (fabs(wrapped - lo) < tol) {
        wrapped = lo;
    }

    return wrapped;
}
Output::Output(MPI_Comm *comm_, Input* input) : comm(comm_)
{
  sprintf(full_filename_base, "%s/%s", input->file.foldername, input->file.filename_base);
  sprintf(filename_base, "%s", input->file.filename_base);
  char f1[256];
  sprintf(f1, "%s_summary.txt", full_filename_base);

  summaryfile.open(f1, ios::out);
  summaryfile << "Computation started at:" << endl;
  summaryfile << "  " << getCurrentDateTime() << endl;
  summaryfile << "Using Code Revision:" << endl;
  summaryfile << endl;
  summaryfile << "Inputs" << endl;
  summaryfile << "  T0 = " << input->file.T0 << endl;
  summaryfile << "  rc = " << input->file.rc << endl;
// summaryfile << "  dt = " << input->file.dt << endl;
// summaryfile << "  t_final = " << input->file.t_final << endl;
  summaryfile << "  N = " << input->file.N << endl;
//  summaryfile << "  xe = " << input->file.xe << endl;
  summaryfile << "  a_M = " << input->file.a_M << endl;
  summaryfile << "  sigma_M = " << input->file.sigma_M << endl;
  summaryfile << "  sigma_H = " << input->file.sigma_H << endl;
  //summaryfile << "  minimize_frequency = " << input->file.minimize_frequency << endl;
  //summaryfile << "  output_frequency = " << input->file.output_frequency << endl;
  //summaryfile << "  Qm = " << input->file.Qm << endl;
  //summaryfile << "  v = " << input->file.v << endl;
  //summaryfile << "  mubd = " << input->file.mubd << endl;

}

Output::~Output()
{
  if(summaryfile.is_open()) summaryfile.close();
}


void Output::output_solution(int iFrame, int iTimeStep, vector<Vec3D> &q_M, vector<Vec3D> &q_H,
                     vector<double> &sigma_M, vector<double> &sigma_H, vector<double> &x, 
                     vector<double> &gamma, vector<double> &f, vector<double> &V_M, 
                     vector<double> &V_H, vector<vector<double> > &pi_M, 
                     vector<vector<double> > &pi_H, 
										 vector<int> &full_H, int &nH_oct, int &nH_tet, Input &input)
{
  int MPI_rank = 0;
  MPI_Comm_rank(*comm, &MPI_rank);
  if(MPI_rank)
    return;
  
  const double kB = 8.6173324e-5;

  // find bounds of simulation box
  //double xmax, xmin, ymax, ymin, zmax, zmin;
  //xmax = q_M[0][0]; xmin = q_M[0][0];
  //ymax = q_M[0][1]; ymin = q_M[0][1];
  //zmax = q_M[0][2]; zmin = q_M[0][2];

  //for(int i=0; i<q_M.size(); i++) {
  //  if(q_M[i][0]>xmax) xmax = q_M[i][0];
  //  if(q_M[i][0]<xmin) xmin = q_M[i][0];
  //  if(q_M[i][1]>ymax) ymax = q_M[i][1];
  //  if(q_M[i][1]<ymin) ymin = q_M[i][1];
  //  if(q_M[i][2]>zmax) zmax = q_M[i][2];
  //  if(q_M[i][2]<zmin) zmin = q_M[i][2];
 // }

  // create solution file
  char full_fname[128];

  sprintf(full_fname, "%s.%d.dat", full_filename_base, iFrame);
  ofstream file(full_fname, ios::out);

  // write time step
  file << "ITEM: TIMESTEP" << endl;
  file << iTimeStep  << endl;

  // write number of atomic sites
  file << "ITEM: NUMBER OF ATOMS" << endl;
  file <<  q_M.size()+q_H.size() << endl;
//double a_M  = input.file.a_M;
//double Lbox = input.file.N * a_M;
double a_M = input.file.a_M;
double Lx =input.file.Fxx_current *input.file.N *a_M;
//double Lx = input.file.N * a_M;   // compressed
double Ly = input.file.N * a_M;
double Lz = input.file.N * a_M;

double xlo=-0.5*Lx, xhi=0.5*Lx;
double ylo=-0.5*Ly, yhi=0.5*Ly;
double zlo=-0.5*Lz, zhi=0.5*Lz;
//double xlo = -0.5 * Lbox;
//double xhi =  0.5 * Lbox;

//double ylo = -0.5 * Lbox;
//
//double yhi =  0.5 * Lbox;

//double zlo = -0.5 * Lbox;
//double zhi =  0.5 * Lbox;

//double a_M  = input.file.a_M;
//  double Lbox = input.file.N * a_M;
// double xlo = xmin - 0.25*a_M, xhi = xlo + Lbox;
// double ylo = ymin - 0.25*a_M, yhi = ylo + Lbox;
// double zlo = zmin - 0.25*a_M, zhi = zlo + Lbox;

  // write bounds of simulation box
  file << "ITEM: BOX BOUNDS pp pp pp" << endl;
  file << scientific << xlo << " " << scientific << xhi  << endl;
  file << scientific << ylo << " " << scientific << yhi  << endl;
  file << scientific << zlo << " " << scientific << zhi  << endl;

  // write properties on each site
  file << "ITEM: ATOMS id type x y z sigma atomic_fraction chemical_potential formation_energy potential_energy stress_xx stress_xy stress_xz stress_yx stress_yy stress_yz stress_zx stress_zy stress_zz full_occupancy" << endl;

  double V_out = 0.0;
  vector<double > pi_out;
  pi_out.assign(9, 0.0); 

  // M sites
  for(int i=0; i<q_M.size(); i++) {
    if (input.file.output_pot == 1) {
      V_out = V_M[i];
    }    

    if (input.file.output_str == 1) {
      for(int j=0; j<9; j++) {
        pi_out[j] = pi_M[i][j];
      }
    }

 //   file << i+1 << " " << 1 << " " 
 //        << scientific << q_M[i][0] << " " 
 //        << scientific << q_M[i][1] << " " 
 //        << scientific << q_M[i][2] << " " 
// double xout = wrapCoord(q_M[i][0], xlo, Lbox);
//double yout = wrapCoord(q_M[i][1], ylo, Lbox);
//double zout = wrapCoord(q_M[i][2], zlo, Lbox);
double xout = wrapCoord(q_M[i][0], xlo, Lx);
double yout = wrapCoord(q_M[i][1], ylo, Ly);
double zout = wrapCoord(q_M[i][2], zlo, Lz);

file << i+1 << " " << 1 << " " 
     << scientific << xout << " " 
     << scientific << yout << " " 
     << scientific << zout << " "
         << scientific << sigma_M[i] << " " 
         << scientific << 1.0 << " " 
         << scientific << 0.0 << " " 
         << scientific << 0.0 << " " 
         << scientific << V_out << " " 
         << scientific << pi_out[0] << " "
         << scientific << pi_out[1] << " "
         << scientific << pi_out[2] << " "
         << scientific << pi_out[3] << " "
         << scientific << pi_out[4] << " "
         << scientific << pi_out[5] << " "
         << scientific << pi_out[6] << " "
         << scientific << pi_out[7] << " "
         << scientific << pi_out[8] << " "
         << 1 << endl;
  }
  // H octahedral sites
  for(int i=0; i<nH_oct; i++) {
    if (input.file.output_pot == 1) {
      V_out = V_H[i];
    }

    if (input.file.output_str == 1) {
      for(int j=0; j<9; j++) {
        pi_out[j] = pi_H[i][j];
      }
    }

//    file << i+1+q_M.size() << " " << 2 << " " 
//         << scientific << q_H[i][0] << " " 
//         << scientific << q_H[i][1] << " " 
//         << scientific << q_H[i][2] << " "
//double xout = wrapCoord(q_H[i][0], xlo, Lbox);
//double yout = wrapCoord(q_H[i][1], ylo, Lbox);
//double zout = wrapCoord(q_H[i][2], zlo, Lbox);
double xout = wrapCoord(q_H[i][0], xlo, Lx);
double yout = wrapCoord(q_H[i][1], ylo, Ly);
double zout = wrapCoord(q_H[i][2], zlo, Lz);


file << i+1+q_M.size() << " " << 2 << " " 
     << scientific << xout << " " 
     << scientific << yout << " " 
     << scientific << zout << " " 
         << scientific << sigma_H[i] << " " 
         << scientific << x[i] << " " 
         << scientific << gamma[i]*kB*input.file.T0 << " " 
         << scientific << f[i] << " " 
         << scientific << V_out << " "
         << scientific << pi_out[0] << " "
         << scientific << pi_out[1] << " "
         << scientific << pi_out[2] << " "
         << scientific << pi_out[3] << " "
         << scientific << pi_out[4] << " "
         << scientific << pi_out[5] << " "
         << scientific << pi_out[6] << " "
         << scientific << pi_out[7] << " "
         << scientific << pi_out[8] << " "
         << full_H[i] << endl;
  }

	// H tetrahedral sites
	for(int i=nH_oct; i<nH_oct+nH_tet; i++) {
    if (input.file.output_pot == 1) {
      V_out = V_H[i];
    }

    if (input.file.output_str == 1) {
      for(int j=0; j<9; j++) {
        pi_out[j] = pi_H[i][j];
      }
    }


	//	file << i+1+q_M.size() << " " << 3 << " " 
        // << scientific << q_H[i][0] << " " 
        // << scientific << q_H[i][1] << " " 
        // << scientific << q_H[i][2] << " " 
        // << scientific << sigma_H[i] << " "
        double xout = wrapCoord(q_H[i][0], xlo, Lx);
    double yout = wrapCoord(q_H[i][1], ylo, Ly);
    double zout = wrapCoord(q_H[i][2], zlo, Lz);

		file << i+1+q_M.size() << " " << 3 << " " 
         << scientific << xout << " " 
         << scientific << yout << " " 
         << scientific << zout << " " 
         << scientific << sigma_H[i] << " " 
         << scientific << x[i] << " " 
         << scientific << gamma[i]*kB*input.file.T0 << " " 
         << scientific << f[i] << " " 
         << scientific << V_out << " "
         << scientific << pi_out[0] << " "
         << scientific << pi_out[1] << " "
         << scientific << pi_out[2] << " "
         << scientific << pi_out[3] << " "
         << scientific << pi_out[4] << " "
         << scientific << pi_out[5] << " "
         << scientific << pi_out[6] << " "
         << scientific << pi_out[7] << " "
         << scientific << pi_out[8] << " "
         << full_H[i] << endl;
  }
  return;
}


// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string Output::getCurrentDateTime() 
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}
