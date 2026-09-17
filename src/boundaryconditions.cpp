#include <iostream>
#include <vector>
#include <math.h>
#include <mpi.h>
#include "Vector3D.h"
#include "input.h"
using namespace std;

int findAdsorptionSites(Input &input, vector<double> &gamma, vector<int> &HSubsurf,
                        vector<int> &HAdsorp)
{
  /*
  int &N_ads = input.file.N_ads;
  HAdsorp.clear();

  vector<pair<double,int> > vect;
  for (int i=0; i<HSubsurf.size(); i++) {
      vect.push_back(make_pair(gamma[HSubsurf[i]],HSubsurf[i]));
  }
  // Using simple sort() function to sort based on the first column
  sort(vect.begin(), vect.end());  

  for (int i=0; i<N_ads; i++) {
      HAdsorp.push_back(vect[i].second);
  }

  vect.clear();
*/
  return 0;
}

int enforceBoundaryConditions(Input &input, double dt, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0,
                              vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<double> &f, vector<int> &HAdsorp, 
                              double &gammabd, vector<int> &full_H, int nH_oct, vector<vector<int> > &HH)
{
    const double xUb = 1.0-1.0e-16; // The upper bound of hydrogen atomic fraction in the sample
    const double xLb = 1.0e-16; // The lower bound of hydrogen atomic fraction in the sample
		const double kB = 8.6173324e-5;
		const double r_bd = 3.5;
    double &T = input.file.T0;
		//double x_boundary = min(input.file.xbd, xUb);
     
    for (int i=0; i<HAdsorp.size(); i++) {
      x0[HAdsorp[i]] = 1.0;
      gamma[HAdsorp[i]] = gammabd;
      //f[HAdsorp[i]] = kB*T*(gammabd - log(xUb/(1.0-xUb)));
      full_H[HAdsorp[i]] = 1;
/*
      for(int k=0; k<HH[HAdsorp[i]].size(); k++) {
        int j = HH[HAdsorp[i]][k];
        Vec3D r = q_H0[HAdsorp[i]]-q_H0[j];
        if(r.norm()>r_bd) continue;
        x0[j] = x_boundary;
        gamma[j] = gammabd;
        f[j] = kB*T*(gammabd - log(xUb/(1.0-xUb)));
        full_H[j] = 1;
      }
*/
    }

/*
    int tmp = 3855;
    
    double J_ij = 0.0;
    
    //J_ij = v*exp(-(Qm+0.5*(f[tmp]+0.0))/(kB*T))
    //                    *(1.0-x0[tmp])*(1.0-xUb)
    //                    *(exp(gammabd)-exp(gamma[tmp]));
    //x0[tmp] += dt*J_ij;
    //if (x0[tmp]>xUb) {
      x0[tmp] = xUb;
      gamma[tmp] = gammabd;
      f[tmp] = kB*T*(gammabd - log(xUb/(1.0-xUb)));
      full_H[tmp] = 1;
    //}
  */    
    //cout << "Adsorption rate: " << J_ij << endl;
    //cout << "Adsorption fraction: " << x0[tmp] << endl;
    /*
		for (int i=0; i<HSubsurf.size(); i++) {
			if (HSubsurf[i] < nH_oct) {
				gamma[HSubsurf[i]] = gammabd;
				full_H[HSubsurf[i]] = 1;
					x0[HSubsurf[i]] = x_boundary;
			}
		}
    */

    return 0;
}

/*
int enforceBoundaryConditions(Input &input, double dt, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0,
                              vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<int> &HSubsurf, double &gammabd, vector<int> &full_H)
{
  // TODO: Appropriate boundary conditions need to be specified in "Input" and enforced here.

    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    int MPI_rank, MPI_size;
    MPI_Comm_rank(comm, &MPI_rank);
    MPI_Comm_size(comm, &MPI_size);

    const double kB = 8.6173324e-5;
    const double xUb = 1.0-1.0e-16; // The upper bound of hydrogen atomic fraction in the sample
    const double xLb = input.file.xe; // The lower bound of hydrogen atomic fraction in the sample
    
    if(input.file.subsurf_equil == 1) {
        double &Bad = input.file.Bad; // Bondwise adsorption coefficient

        int nHSubsurf  = HSubsurf.size();
        int locSize; // local size for each core except the last one
        int lastSize; // local size for the last one

        double *send_x, *recv_x;
        send_x = new double[nHSubsurf];
        recv_x = new double[nHSubsurf];

        for(int i=0; i<nHSubsurf; i++)
            send_x[i] = 0.0;

        if(nHSubsurf%MPI_size == 0) {
            locSize = nHSubsurf/MPI_size;
            lastSize = nHSubsurf/MPI_size;
        } else {
            locSize = nHSubsurf/(MPI_size-1);
            lastSize = nHSubsurf%(MPI_size-1);
        }

        if((MPI_rank+1) != MPI_size) { // not the last core
            for(int ii=0; ii<locSize; ii++) {
                int i = MPI_rank*locSize + ii;
                if(full_H[HSubsurf[i]]==1) continue; // skip the site which is fully occupied by H
                double x_deriv = gammabd-gamma[HSubsurf[i]];
                send_x[i] = x0[HSubsurf[i]] + dt*x_deriv*kB*Bad/2.0;
            }
        } else { // the last core
            for(int ii=0; ii<lastSize; ii++) {
                int i = MPI_rank*locSize + ii;
                if(full_H[HSubsurf[i]]==1) continue; // skip the site which is fully occupied by H
                double x_deriv = gammabd-gamma[HSubsurf[i]];
                send_x[i] = x0[HSubsurf[i]] + dt*x_deriv*kB*Bad/2.0;
            }
        }

        MPI_Barrier(comm);
        MPI_Allreduce(send_x, recv_x, nHSubsurf, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Update the boundary values in x
        for(int i=0; i<nHSubsurf; i++) {
            if(full_H[HSubsurf[i]]==1) continue; // skip the site which is fully occupied by H 
            x0[HSubsurf[i]] = recv_x[i];
            if(x0[HSubsurf[i]]<xLb) x0[HSubsurf[i]] = xLb;
            else if(x0[HSubsurf[i]]>xUb) x0[HSubsurf[i]] = xUb;
        }

        delete[] send_x;
        delete[] recv_x;
    }
    else {
       for(int i=0; i<HSubsurf.size(); i++) {
            gamma[HSubsurf[i]] = gammabd;
            full_H[HSubsurf[i]] = 1;
        }
    }
 
    return 0;
}
*/
