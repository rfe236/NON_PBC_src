/* 
* File: read_setfl.cpp
* Description: read data of tabulated EAM potential from a setfl format file
*              details about setfl format can be found at
*              https://lammps.sandia.gov/doc/pair_eam.html
*              https://spasmmini.weebly.com/format-of-tabulated-eam-potential.html
* Author: Xingsheng Sun
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "cubic_spline.h"
#include "eam_setfl.h"

using namespace std;

void EAM_setfl_cubic_spline(string eam_setfl_file, vector<CubicSpline> &CSembed, vector<CubicSpline> &CSrho, vector<CubicSpline> &CSpair)
{
	int Nelements, Nrho, Nr;
	double drho, dr, cutoff;
  int i, j, k;
	FILE* f_setfl;
	double *F, *rho, *phi;

	CubicSpline CS_tmp;

	f_setfl = fopen(eam_setfl_file.c_str(), "r" );

	if (f_setfl == NULL) {
		char msg[160];
		sprintf(msg, "No potential file %s", (char *)eam_setfl_file.c_str());
		perror(msg);
		return;
	}

	//if(myrank==0)
	//printf("Reading eam potential file (in setfl format): %s\n", eam_setfl_file.c_str());

	// skip the first three lines
	fscanf(f_setfl, "%*[^\n]\n");
	fscanf(f_setfl, "%*[^\n]\n");
	fscanf(f_setfl, "%*[^\n]\n");
	
	// read data on the fourth line	
	fscanf(f_setfl, "%d", &Nelements); // read number of elements
	//printf("%d\n", Nelements);
	string Etypes[Nelements]; // used to store element types
	for(i = 0; i < Nelements; i++ ) {
		char str[16];
		fscanf(f_setfl, "%s", str); // read element types
		Etypes[i] = str;
	}
	/*for(i = 0; i < Nelements; i++ ) {
		std::cout << Etypes[i] << std::endl;
	}*/

	// read data on the fifth line
	fscanf(f_setfl, "%d %lf %d %lf %lf", &Nrho, &drho, &Nr, &dr, &cutoff);
	// printf("%d %e %d %e %e\n", Nrho, drho, Nr, dr, cutoff);

	int Anum[Nelements]; // used to store atomic number
	int Amass[Nelements]; // used to store atomic mass
	double Lconst[Nelements]; // used to store lattice constant
	string Ltype[Nelements]; // used to store lattice type

	F = new double[Nrho];
	rho = new double[Nr];
	phi = new double[Nr];

	// read embedding energy and electron density
	for(i = 0; i < Nelements; i++) {
		// read basic information for each element
		char str[16];
		fscanf(f_setfl, "%d %lf %lf %s", &Anum[i], &Amass[i], &Lconst[i], str);
		Ltype[i] = str;

		// read embedding energy for each element
		for(j = 0; j < Nrho; j++) {
			F[j] = 0.0;
			fscanf(f_setfl, "%lf", &F[j]); // read embedding energy for i-th element
		}
		//TODO: calculate spline coefficients
		init_spline(&CS_tmp, Nrho-1, drho, F);
		CSembed.push_back(CS_tmp);
		// read electron density for each element
		for(j = 0; j < Nr; j++) {
			rho[j] = 0.0;
			fscanf(f_setfl, "%lf", &rho[j]); // read embedding energy for i-th element
		}
		//TODO: calculate spline coefficients
		init_spline(&CS_tmp, Nr-1, dr, rho);
		CSrho.push_back(CS_tmp);
		//int ii = 250;
		//printf("%e %e\n", F[ii], rho[ii]);
		//printf("%d %e %e %e %e %e %e\n", CSembed.n, CSembed.dx, CSembed.x[ii], CSembed.d[ii],      
						//CSembed.c[ii], CSembed.b[ii], CSembed.a[ii]);
		//printf("%d %e %e %e %e %e %e\n", CSrho.n, CSrho.dx, CSrho.x[ii], CSrho.d[ii], 
						//CSrho.c[ii], CSrho.b[ii], CSrho.a[ii]);
	}
	
	// read pair energy
	for(i = 0; i < Nelements; i++) {
		for(j = 0; j < i+1; j++) {
			// read pair energy for i-j pair
			// setfl only stores pair energy for i>=j in the following sequence
			// phi_00 phi_10 phi_11 phi_20 phi_21 phi_22 phi_30 phi_31 phi_32 phi_33 ...
			for(k = 0; k < Nr; k++) {
				phi[k] = 0.0;
				fscanf(f_setfl, "%lf", &phi[k]);
        phi[k] /= max(1.0e-12, k*dr); // in the setfl format, pair potential is in the format of r*phi(r), where, r is inter-atomic distance, 
			}
			//TODO: calculate spline coefficients
			init_spline(&CS_tmp, Nr-1, dr, phi);
			CSpair.push_back(CS_tmp);
			//int ii = 250;
			//printf("%e\n", phi[ii]);
			//printf("%d %e %e %e %e %e %e\n", CSpair.n, CSpair.dx, CSpair.x[ii], CSpair.d[ii],
			        //CSpair.c[ii], CSpair.b[ii], CSpair.a[ii]);
		}
	}

	fclose(f_setfl);
	delete F;
  delete rho;
  delete phi;
  return;
}


