#include <iostream>
#include <vector>
#include <math.h>
#include "Vector3D.h"
#include "input.h"
#include <limits>
#include <cmath>
#include <mpi.h>
using namespace std;
void findRigidAtoms(vector<Vec3D> &q_M0, vector<int> &rigidAtoms)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
    int nM = q_M0.size();
    rigidAtoms.assign(4, -1);

    if (nM < 3) {
        cerr << "ERROR: Need at least 3 host atoms to fix rigid-body motion." << endl;
        exit(-1);
    }
double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    for (int i = 0; i < nM; i++) {
        cx += q_M0[i][0];
        cy += q_M0[i][1];
        cz += q_M0[i][2];
    }

    cx /= nM;
    cy /= nM;
    cz /= nM;
int atomA = -1;
    double bestA = numeric_limits<double>::max();

    for (int i = 0; i < nM; i++) {
        double dx = q_M0[i][0] - cx;
        double dy = q_M0[i][1] - cy;
        double dz = q_M0[i][2] - cz;

        double dist2 = dx*dx + dy*dy + dz*dz;

        if (dist2 < bestA) {
            bestA = dist2;
            atomA = i;
        }
    }

    Vec3D A = q_M0[atomA];
int atomB = -1;
    double bestB = numeric_limits<double>::max();
    double bestB_sep = -1.0;

    for (int i = 0; i < nM; i++) {
        if (i == atomA) continue;

        double dx = q_M0[i][0] - A[0];
        double dy = q_M0[i][1] - A[1];
        double dz = q_M0[i][2] - A[2];

        double obj = dy*dy + dz*dz;   // distance from x-line
        double sep = fabs(dx);        // distance along x

        if (obj < bestB ||
            (fabs(obj - bestB) < 1.0e-10 && sep > bestB_sep)) {
            bestB = obj;
            bestB_sep = sep;
            atomB = i;
        }
    }
  int atomC = -1;
    double bestC = numeric_limits<double>::max();
    double bestC_sep = -1.0;

    for (int i = 0; i < nM; i++) {
        if (i == atomA || i == atomB) continue;

        double dx = q_M0[i][0] - A[0];
        double dy = q_M0[i][1] - A[1];
        double dz = q_M0[i][2] - A[2];

        double obj = dx*dx + dz*dz;   // distance from y-line
        double sep = fabs(dy);        // distance along y

        if (obj < bestC ||
            (fabs(obj - bestC) < 1.0e-10 && sep > bestC_sep)) {
            bestC = obj;
            bestC_sep = sep;
            atomC = i;
        }
    }
int atomD = -1;
double bestD = numeric_limits<double>::max();
double bestD_sep = -1.0;

for (int i = 0; i < nM; i++) {
    if (i == atomA || i == atomB || i == atomC) continue;

    double dx = q_M0[i][0] - A[0];
    double dy = q_M0[i][1] - A[1];
    double dz = q_M0[i][2] - A[2];

    double obj = dx*dx + dy*dy;   // distance from z-line
    double sep = fabs(dz);        // distance along z

    if (obj < bestD ||
        (fabs(obj - bestD) < 1.0e-10 && sep > bestD_sep)) {
        bestD = obj;
        bestD_sep = sep;
        atomD = i;
    }
}
if (atomA < 0 || atomB < 0 || atomC < 0 || atomD < 0) {
  if (rank == 0) {
        cerr << "ERROR: Could not find rigid atoms A, B, C, D." << endl;
    }

    MPI_Abort(MPI_COMM_WORLD, 1);
}
rigidAtoms.clear();
rigidAtoms.push_back(atomA);
rigidAtoms.push_back(atomB);
rigidAtoms.push_back(atomC);
rigidAtoms.push_back(atomD);

// if (atomA < 0 || atomB < 0 || atomC < 0) {
//        cerr << "ERROR: Could not find rigid atoms A, B, C." << endl;
//        exit(-1);
//    }

 

if (rank == 0) {
    cout << "DEBUG rigid atoms selected:" << endl;
    cout << "  center = " << cx << " " << cy << " " << cz << endl;

    cout << "  A = " << atomA << " closest to center, fixed xyz at "
         << q_M0[atomA][0] << " "
         << q_M0[atomA][1] << " "
         << q_M0[atomA][2] << endl;

    cout << "  B = " << atomB << " closest to x-line through A, fixed yz at "
         << q_M0[atomB][0] << " "
         << q_M0[atomB][1] << " "
         << q_M0[atomB][2] << endl;

    cout << "      B line distance squared = " << bestB
         << ", x separation = " << bestB_sep << endl;

    cout << "  C = " << atomC << " closest to y-line through A, fixed z at "
         << q_M0[atomC][0] << " "
         << q_M0[atomC][1] << " "
         << q_M0[atomC][2] << endl;

    cout << "      C line distance squared = " << bestC
         << ", y separation = " << bestC_sep << endl;
cout << "  D = " << atomD << " closest to z-line through A, fixed xy at "
     << q_M0[atomD][0] << " "
     << q_M0[atomD][1] << " "
     << q_M0[atomD][2] << endl;

cout << "      D line distance squared = " << bestD
     << ", z separation = " << bestD_sep << endl;
 cout.flush();
}
}



//{
//    int nM = q_M0.size();
//    rigidAtoms.assign(3, -1);

//    if (nM < 3) {
 //       cerr << "ERROR: Need at least 3 host atoms to fix rigid body motion." << endl;
//        exit(-1);
//    }
//  double cx = 0.0;
//    double cy = 0.0;
//    double cz = 0.0;

//    for (int i = 0; i < nM; i++) {
//        cx += q_M0[i][0];
//        cy += q_M0[i][1];
//        cz += q_M0[i][2];
//    }

//    cx /= nM;
//    cy /= nM;
 //   cz /= nM;
 //int atomB = -1;
 //   double bestB = numeric_limits<double>::max();
//   double bestB_sep = -1.0;

//    for (int i = 0; i < nM; i++) {
//        double dx = q_M0[i][0] - cx;
//        double dy = q_M0[i][1] - cy;
//        double dz = q_M0[i][2] - cz;

//        double obj = dy*dy + dz*dz;
//        double sep = fabs(dx);

//        if (obj < bestB || (fabs(obj - bestB) < 1.0e-12 && sep > bestB_sep)) {
//            bestB = obj;
//            bestB_sep = sep;
//            atomB = i;
//        }
//    }
//int atomC = -1;
//    double bestC = numeric_limits<double>::max();
//    double bestC_sep = -1.0;

 //   for (int i = 0; i < nM; i++) {
 //       if (i == atomB) continue;

//        double dx = q_M0[i][0] - cx;
//        double dy = q_M0[i][1] - cy;
//        double dz = q_M0[i][2] - cz;

//        double obj = dx*dx + dz*dz;
//        double sep = fabs(dy);

//        if (obj < bestC || (fabs(obj - bestC) < 1.0e-12 && sep > bestC_sep)) {
//            bestC = obj;
//            bestC_sep = sep;
//            atomC = i;
//        }
//    }
// int atomD = -1;
//    double bestD = numeric_limits<double>::max();
//    double bestD_sep = -1.0;

//    for (int i = 0; i < nM; i++) {
//        if (i == atomB || i == atomC) continue;

//        double dx = q_M0[i][0] - cx;
//        double dy = q_M0[i][1] - cy;
//        double dz = q_M0[i][2] - cz;

//        double obj = dx*dx + dy*dy;
//        double sep = fabs(dz);

//        if (obj < bestD || (fabs(obj - bestD) < 1.0e-12 && sep > bestD_sep)) {
//            bestD = obj;
//            bestD_sep = sep;
//            atomD = i;
//        }
//    }

//    if (atomB < 0 || atomC < 0 || atomD < 0) {
//        cerr << "ERROR: Could not find rigid fixing atoms B, C, D." << endl;
//        exit(-1);
//    }

//    rigidAtoms[0] = atomB;
//    rigidAtoms[1] = atomC;
//    rigidAtoms[2] = atomD;

//    cout << "DEBUG rigid fixing atoms using virtual center:" << endl;
//    cout << "  center = " << cx << " " << cy << " " << cz << endl;

//    cout << "  B = " << atomB << " along x, fixed y,z at "
//         << q_M0[atomB][0] << " "
//         << q_M0[atomB][1] << " "
//         << q_M0[atomB][2] << endl;

//    cout << "  C = " << atomC << " along y, fixed x,z at "
//         << q_M0[atomC][0] << " "
//         << q_M0[atomC][1] << " "
 //        << q_M0[atomC][2] << endl;

//     cout << "  D = " << atomD << " along z, fixed x,y at "
//         << q_M0[atomD][0] << " "
//         << q_M0[atomD][1] << " "
//         << q_M0[atomD][2] << endl;
//}
                                                   
void initializeStateVariables(Input &input, vector<vector<Int3> > &SS1, vector<vector<Int3> > &SS2,
                              vector<Vec3D> &SS3,
															vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, vector<double> &sigma_M0,
                              vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma,
														 	vector<double> &f, vector<double> &V_M, vector<double> &V_H,	
                              vector<vector<double> > &pi_M, vector<vector<double> > &pi_H,
                              vector<int> &HSubsurf, 
                              double &gammabd, vector<int> &full_H, vector<int> &HType,
	
  														int &nH_oct, int &nH_tet, vector<int> &rigidAtoms) 
{
    // initialize boundary conditions
    double &T = input.file.T0;
    const double kB = 8.6173324e-5; // (eV/K)
    const double cutoff = 1.0e-20; // The minima of hydrogen atomic fraction
    const double xUb = 1.0-1.0e-16;
    const double xLb = 1.0e-16;

    //gammabd = input.file.mubd/(kB*T);

  if(!(q_M0.empty() && q_H0.empty() && sigma_M0.empty() && sigma_H0.empty() 
       && x0.empty() && gamma.empty() && f.empty() && V_M.empty() && V_H.empty() 
       && pi_M.empty() && pi_H.empty() && HSubsurf.empty())) {
    cerr << "WARNING: some vectors are not empty before initialization!" << endl;
    q_M0.clear();
    q_H0.clear();
    sigma_M0.clear();
    sigma_H0.clear();
    x0.clear();
    gamma.clear();
    f.clear();
		V_M.clear();
		V_H.clear();
    pi_M.clear();
    pi_H.clear();
    HSubsurf.clear();
  }

  // INITIALIZE q_M0
  q_M0.push_back(Vec3D(0.0,0.0,0.0)); //first atomic site

  int p = 1; // store index 
  // loop through shells to add sites
  for(int i=0; i<SS1.size(); i++)
    for(int j=0; j<SS1[i].size(); j++) {
      Int3 &site = SS1[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_M0.push_back(q_site);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //rhombic dodecahedron
        if(fabs(site[0])+fabs(site[1])<=N && fabs(site[1])+fabs(site[2])<=N && fabs(site[0])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_M0.push_back(q_site);
          p++;
        }
      }
      else if(input.file.sample_shape == 3) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_M0.push_back(q_site);
          p++;
        }
      }
      else if(input.file.sample_shape == 4) { // sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_M0.push_back(q_site);
          p++;
        }
      }
    }
  
  // INITIALIZE q_H0 for octahetral interstitial sites
  p = 0; // store index 
  for(int i=0; i<SS2.size(); i++)
    for(int j=0; j<SS2[i].size(); j++) {
      Int3 &site = SS2[i][j];
      int N = input.file.N;
      if(input.file.sample_shape == 1) { //cube
        if(site[0]>-N && site[0]<=N && site[1]>-N && site[1]<=N && site[2]>-N && site[2]<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_H0.push_back(q_site);
					HType.push_back(0);
          p++;
        }
      }
      else if(input.file.sample_shape == 2) { //rhombic dodecahedron
        if(fabs(site[0])+fabs(site[1])<=N && fabs(site[1])+fabs(site[2])<=N && fabs(site[0])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_H0.push_back(q_site);
          HType.push_back(0);
          p++;
        }
      }
      else if(input.file.sample_shape == 3) { //octahedron
        if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_H0.push_back(q_site);
					HType.push_back(0);
          p++;
        }
      }
      else if(input.file.sample_shape == 4) { //sphere
        if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
          Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
          q_site *= 0.5*input.file.a_M;
          q_H0.push_back(q_site);
					HType.push_back(0);
          p++;
        }
      }
    }

	nH_oct = q_H0.size();


	// INITIALIZE q_H0 for tetrahedral interstitial sites
	// At this time, p should be equal to nH_oct

	/*for(int i=0; i<SS3.size(); i++) {
		Vec3D &site = SS3[i];
		int N = input.file.N;
		if(input.file.sample_shape == 1) { //cube
			if(site[0]>-(N-0.75) && site[0]<=N && site[1]>-(N-0.75) && site[1]<=N && site[2]>-(N-0.75) && site[2]<=N) {
				Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
				q_site *= 0.5*input.file.a_M;
				q_H0.push_back(q_site);
				HType.push_back(1);
				p++;
			}
		}
*/
// Initialize tetrahedral interstitial sites
if (input.file.sample_shape == 1) {

    const int N = input.file.N;
    const double a = input.file.a_M;
    const double L = N * a;
    const double lo = -0.5 * L;
const double tetOffset[2] = {0.25, 0.75};

  for (int ix = 0; ix < N; ix++) {
  for (int iy = 0; iy < N; iy++) {
  for (int iz = 0; iz < N; iz++) {

  for (int ox = 0; ox < 2; ox++) {
  for (int oy = 0; oy < 2; oy++) {
  for (int oz = 0; oz < 2; oz++) {

 const double x = lo + (ix + tetOffset[ox]) * a;

 const double y = lo + (iy + tetOffset[oy]) * a;

   const double z = lo + (iz + tetOffset[oz]) * a;

    q_H0.push_back(Vec3D(x, y, z));
 HType.push_back(1);

 p++;
}
 }
 }
 }
}
}
}
else {
 for (int i = 0; i < static_cast<int>(SS3.size()); i++) {

        Vec3D &site = SS3[i];
        const int N = input.file.N;

        if (input.file.sample_shape == 2) {
//rhombic
if (fabs(site[0]) + fabs(site[1]) <= N &&
                fabs(site[1]) + fabs(site[2]) <= N &&
                fabs(site[0]) + fabs(site[2]) <= N)
            {
                Vec3D q_site(
                    static_cast<double>(site[0]),
                    static_cast<double>(site[1]),
                    static_cast<double>(site[2])
                );

                q_site *= 0.5 * input.file.a_M;

                q_H0.push_back(q_site);
                HType.push_back(1);
                p++;
            }
        }
        else if (input.file.sample_shape == 3) {
//sphere
  if (site[0] * site[0] +
                site[1] * site[1] +
                site[2] * site[2] <= N * N)
            {
                Vec3D q_site(
                    static_cast<double>(site[0]),
                    static_cast<double>(site[1]),
                    static_cast<double>(site[2])
                );

                q_site *= 0.5 * input.file.a_M;

                q_H0.push_back(q_site);
                HType.push_back(1);
                p++;
            }
        }
    }
}   
/*    else if(input.file.sample_shape == 2) { //rhombic dodecahedron
      if(fabs(site[0])+fabs(site[1])<=N && fabs(site[1])+fabs(site[2])<=N && fabs(site[0])+fabs(site[2])<=N) {
        Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
        q_site *= 0.5*input.file.a_M;
        q_H0.push_back(q_site);
        HType.push_back(1);
        p++;
      }
    }
		else if(input.file.sample_shape == 3) { //octahedron
			if(fabs(site[0])+fabs(site[1])+fabs(site[2])<=N) {
				Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
				q_site *= 0.5*input.file.a_M;
				q_H0.push_back(q_site);
				HType.push_back(1);
				p++;
			}
		}
		else if(input.file.sample_shape == 4) { //sphere
			if(site[0]*site[0] + site[1]*site[1] + site[2]*site[2] <= N*N) {
				Vec3D q_site((double)site[0], (double)site[1], (double)site[2]);
				q_site *= 0.5*input.file.a_M;
				q_H0.push_back(q_site);
				HType.push_back(1);
				p++;
			}
		}
	}

*/

  int nM = q_M0.size();
  int nH  = q_H0.size();
nH_tet = nH - nH_oct;
// Automatically select atoms A, B, and C
findRigidAtoms(q_M0, rigidAtoms);

//rigidAtoms[1] = -1;
//    rigidAtoms[2] = -1;
//    rigidAtoms[3] = -1;
// ===== UNIAXIAL COMPRESSION ALONG X (PBC: scale about origin) =====
//for (int i=0; i<nM; i++) q_M0[i][0] = FXX_STRAIN * q_M0[i][0];   // y,z untouched
//for (int i=0; i<nH; i++) q_H0[i][0] = FXX_STRAIN * q_H0[i][0];

  // Initialize sigma, x, gamma, f
  sigma_M0.assign(nM, input.file.sigma_M);
  sigma_H0.assign(nH, input.file.sigma_H);
  x0.assign(nH, min(max(input.file.x_ini, xLb), xUb));
	
  /*
  for(int i=0; i<nH_oct; i++) {
		sigma_H0.push_back(input.file.sigma_H);		
    if(input.file.abs==0)	x0.push_back(xUb);
    else x0.push_back(xLb);
	}
	for(int i=nH_oct; i<nH; i++) {
		sigma_H0.push_back(input.file.sigma_H);
		if(input.file.abs==0) x0.push_back(xUb);
    else x0.push_back(xLb);
	}
  */
  gamma.assign(nH, 0.0);
  f.assign(nH, 0.0);
  full_H.assign(nH, 0);

  if (input.file.output_pot == 1) {
	  V_M.assign(nM, 0.0);
	  V_H.assign(nH, 0.0);
  }
  
  if (input.file.output_str == 1) {
    pi_M.assign(nM, vector<double>(9, 0.0));
    pi_H.assign(nH, vector<double>(9, 0.0));
  }

/*
  // debug
  for(int i=0; i<nM; i++) 
    sigma_M0[i] = (double)i;
  for(int i=0; i<nH; i++) 
    sigma_H0[i] = (double)(i+nM);
*/

return;
}




