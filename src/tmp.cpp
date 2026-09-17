/**********************************************************************************
 * Copyright © Kevin G. Wang, Xingsheng Sun, 2015
 * (1) Redistribution and use in source and binary forms, with or without modification,
 *     are permitted, provided that this copyright notice is retained.
 * (2) Use at your own risk.
 **********************************************************************************/
#include <stdio.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include "Vector3D.h"
#include "input.h"
#include "output.h"
#include "minimizer.h"
#include "cubic_spline.h"
#include "eam_setfl.h"
using namespace std;

//--------------------------------------------------------------
void generateShells(int NC, vector<vector<Int3> > *shells);
void generateOctahedralInterstitialShells(int NC, vector<vector<Int3> > *shells);
void generateTetrahedralInterstitialSites(int NC, vector<Vec3D> *sites);
void initializeStateVariables(Input &input, vector<vector<Int3> > &SS1, vector<vector<Int3> > &SS2, 
                              vector<Vec3D> &SS3,
															vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, vector<double> &sigma_M0, 
                              vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma,
                              vector<double> &f,
														 	vector<double> &V_M, vector<double> &V_H,	
                              vector<vector<double> > &pi_M, vector<vector<double> > &pi_H,
                              vector<int> &HSubsurf, 
                              double &gammabd, vector<int> &full_H, vector<int> &HType,
															int &nH_oct, int &nH_tet); 
int findNeighbors(double a_M, double rc, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, 
                  vector<vector<int> > &MM, vector<vector<int> > &MH, 
                  vector<vector<int> > &HM, vector<vector<int> > &HH);
int findAdsorptionSites(Input &input, vector<double> &gamma, vector<int> &HSubsurf,
                        vector<int> &HAdsorp);
int enforceBoundaryConditions(Input &input, double dt, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0,
                              vector<double> &sigma_M0, vector<double> &sigma_H0, vector<double> &x0, 
                              vector<double> &gamma, vector<double> &f, vector<int> &HAdsorp, 
                              double &gammabd, vector<int> &full_H, int nH_oct, vector<vector<int> > &HH);
int integrateDiffusion(Input &input, double gammabd, vector<vector<int> > &HH, vector<double> &gamma, 
                       vector<double> &f, vector<double> &x0, 
                       vector<double> &x, vector<int> &full_H, vector<Vec3D> &q_H0, vector<double> &sigma_H0);
int calculateMacroResults(vector<double> &x, double &xH, 
                          double &obj_fun, double &free_energy, double &T, int nM);
void readResults(Input &input, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, vector<double> &sigma_M0,
                 vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, 
                 vector<double> &f, vector<int> &full_H, double &t);
void gaussianQuadrature(int dim, vector<vector<Vec3D> > &QP);
void printLogo(MPI_Comm comm);
//--------------------------------------------------------------


//--------------------------------------------------------------
// Main Function
//--------------------------------------------------------------
int main(int argc, char* argv[])
{
  clock_t start_time = clock(); //for timing purpose only
  //--------------------------------------------------------------
  // Setup MPI 
  //--------------------------------------------------------------
  MPI_Init(&argc, &argv);
  MPI_Comm comm;
  comm = MPI_COMM_WORLD;
  int MPI_rank, MPI_size;
  MPI_Comm_rank(comm, &MPI_rank);
  MPI_Comm_size(comm, &MPI_size);
  MPI_Barrier(comm);
  printLogo(comm);

  //--------------------------------------------------------------
  // Inputs
  //--------------------------------------------------------------
  Input input;
  input.readCmdLine(argc, argv);
  input.readCmdFile();

//  input.initialize(300.0/*T0*/, 1.0e-5/*eps*/, 5.35/*rc*/, 6/*N*/, 0.3/*xe*/, 4.3/*a_M*/,
//                   0.1/*sigma_M*/, 0.1/*sigma_H*/, 1.0e-6/*dt*/, 1.0e-5/*t_final*/,
//                   1/*output_frequency*/, "results"/*folder*/, "sol"/*filename_base*/);

  Output output(&comm, &input);

  //--------------------------------------------------------------
  // Create Shells
  //--------------------------------------------------------------
  int NC = 1.5*input.file.N*input.file.N + 4; //number of shells to be computed. TODO: enough?
  vector<vector<Int3> > SS1; // M
  vector<vector<Int3> > SS2; // H octahedral
	vector<Vec3D> SS3; // H tetrahedral
  generateShells(NC,&SS1);
  generateOctahedralInterstitialShells(NC,&SS2);
	//generateTetrahedralInterstitialSites(input.file.N,&SS3);
 
  if(!MPI_rank) {
    if(input.file.sample_shape == 1) {
        cout << "- Shape of nanoparticle: cube." << endl; cout.flush();}
    else if(input.file.sample_shape == 2) {
        cout << "- Shape of nanoparticle: rhombic dodecahedron." << endl; cout.flush();}
    else if(input.file.sample_shape == 3) {
        cout << "- Shape of nanoparticle: octahedron." << endl; cout.flush();}
    else if(input.file.sample_shape == 4) {
        cout << "- Shape of nanoparticle: sphere." << endl; cout.flush();}
    else {
        cerr << "Error: SampleShapeType in input.st must be 1 (cube), 2 (octahedron) or 3 (sphere)." << endl; cout.flush();
        exit(-1);}
    cout << "- Generated neighbor shells and sites. " << endl; cout.flush();}

  //--------------------------------------------------------------
  // Initialize state: atomic positions, vibration frequencies, H molar fractions
  //--------------------------------------------------------------
  vector<Vec3D> q_M0, q_H0;
  vector<double> sigma_M0, sigma_H0, x0, gamma, f, V_M, V_H;
  vector<vector<double> > pi_M, pi_H;
  vector<int> HSubsurf, full_H, HType;
	int nH_oct, nH_tet;
	// HType defines the type of interstitial sites: 0 for octahetral and 1 for tetrahedral
  double gammabd;
  double T = input.file.T0;

  initializeStateVariables(input, SS1, SS2, SS3, q_M0, q_H0, sigma_M0, sigma_H0, x0, gamma, f, V_M, V_H, pi_M, pi_H, HSubsurf, gammabd, full_H, HType, nH_oct, nH_tet); 

	SS1.clear();
	SS2.clear();
	SS3.clear();

  if(!MPI_rank) {
    cout << "- Initialized atomic sites and state variables. " << endl;
    cout << "    Number of host sites: " << q_M0.size() << ". " << endl;
    cout << "    Number of intersitial sites: " << q_H0.size() << ". " << endl;
		cout << "    Number of octahedral interstitial sites: " << nH_oct << ". " << endl;
		cout << "    Number of tetrahedral interstitial sites: " << nH_tet << ". " << endl;
    cout << "    Number of subsurface interstitial sites: " << HSubsurf.size() << ". " << endl;
    //cout << "- Boundary conditions:" << endl;
    //cout << "    Subsurface thickness: " << input.file.t_subsurf << " angstrom." << endl;
    //cout << "    Nondimensional chemical potential: " << gammabd << "." << endl;
    //cout << "    Chemical potential: " << input.file.mubd << " eV." << endl;
    //cout << "    Atomic fraction: " << input.file.xbd << "." << endl;
		//cout << "    Number of adsorption sites: " << input.file.N_ads << ". " << endl;
    cout.flush();
  }

  //--------------------------------------------------------------
  // Find neighbors of each atom within cut-off distance rc
  // Note: This algorithm does NOT rely on shells or the (i,j,k) indices. Hence allows amorphous layers
  //--------------------------------------------------------------
  double ext_dist = 1.5;
  double rc_ext = input.file.rc + ext_dist;
  vector<vector<int> > MM_ext;
  vector<vector<int> > MH_ext; 
  vector<vector<int> > HM_ext; 
  vector<vector<int> > HH_ext; 
 
  int maxNeighbors_ext = findNeighbors(input.file.a_M, rc_ext, q_M0, q_H0, MM_ext, MH_ext, HM_ext, HH_ext);

  if(!MPI_rank) {
    cout << "- Done with extended neighbor search using a KDTree. " << endl;
    cout << "    Original cutoff distance: " << input.file.rc << " angstrom. " << endl;
    cout << "    Extended cutoff distance: " << rc_ext << " angstrom. " << endl; 
		cout << "    Maximium number of extended neighbors: " << maxNeighbors_ext << ". " << endl;
    cout.flush();}

  //--------------------------------------------------------------
  // Calculate Gaussian Quadrature points and weights with 6 variables
  //--------------------------------------------------------------
  vector<vector<Vec3D> > QPp; // for function with 6 variables.
  gaussianQuadrature(3*2, QPp);
  double QWp = 0.5/(3.0*2.0);

	//--------------------------------------------------------------
	//  Read eam potential file and calculate cubic spline coefficients
	//--------------------------------------------------------------
	vector<CubicSpline> CSembed, CSrho, CSpair;
	EAM_setfl_cubic_spline(input.file.eam_filename, CSembed, CSrho, CSpair);
	if(!MPI_rank) {
		cout << "- Done with calcualtion of cubic spline coefficients for interatomic potential. " << endl;
		cout << "    EAM potential file: " << input.file.eam_filename << "." << endl;	
		cout.flush();}

/*
	for (int i=0; i<2; i++)
	cout <<  CSembed[i].dx << " " << CSembed[i].x[CSembed[i].n-1] << " " << CSembed[i].a[CSembed[i].n-1]
			 << " " << CSembed[i].b[CSembed[i].n-1] << " " << CSembed[i].c[CSembed[i].n-1] 
			 << " " << CSembed[i].d[CSembed[i].n-1] << endl;
	
	for (int i=0; i<2; i++)
  cout <<  CSrho[i].dx << " " << CSrho[i].x[CSrho[i].n-1] << " " << CSrho[i].a[CSrho[i].n-1]
       << " " << CSrho[i].b[CSrho[i].n-1] << " " << CSrho[i].c[CSrho[i].n-1]
       << " " << CSrho[i].d[CSrho[i].n-1] << endl;
	
	for (int i=0; i<3; i++)
  cout <<  CSpair[i].dx << " " << CSpair[i].x[CSpair[i].n-1] << " " << CSpair[i].a[CSpair[i].n-1]
       << " " << CSpair[i].b[CSpair[i].n-1] << " " << CSpair[i].c[CSpair[i].n-1]
       << " " << CSpair[i].d[CSrho[i].n-1] << endl;
*/

  //--------------------------------------------------------------
  // Initialize optimization solver
  //--------------------------------------------------------------
  vector<vector<int> > MM;//M neighbors of M;
  vector<vector<int> > MH; //H neighbors of M;
  vector<vector<int> > HM; //M neighbors of H;
  vector<vector<int> > HH;  // H neighbors of H;
  int err = 0; //error code
  
  Minimizer minimizer(&argc, &argv);
  
  int maxNeighbors = 0;

  //--------------------------------------------------------------
  // Initialization
  //--------------------------------------------------------------
  int iFrame = 0;
  int iTimeStep = 0;
  int startTimeStep = 0;
  double t = 0.0; //current simulation time
  double obj_fun = 0.0;
  double free_energy = 0.0;
  double dt = 0.0; // by XS
  vector<int> HAdsorp;

  const double kB = 8.6173324e-5;
  double mubd_lb = input.file.mubd_lb;
  double mubd_ub = input.file.mubd_ub;
  double mubd_step = input.file.mubd_step;
  double mubd = mubd_lb;
  gammabd = mubd / (kB*input.file.T0);
  int N_mubd = (int) (mubd_ub-mubd_lb) / mubd_step + 1;

  vector<Vec3D> q_M(q_M0), q_H(q_H0); //to be used for storing new states
  vector<double> sigma_M(sigma_M0), sigma_H(sigma_H0), x(x0); //same as above
  
  vector<Vec3D> q_M_ini(q_M0), q_H_ini(q_H0);
  vector<double> sigma_M_ini(sigma_M0), sigma_H_ini(sigma_H0), x_ini(x0), gamma_ini(gamma);  

  double xH = 0.0;
  double xHinterior = 0.0;

  double xH0 = 10.0;

  minimizer.initialize(&input, &MM, &MH, &HM, &HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, &QPp, QWp, &CSembed, &CSrho, &CSpair);

  err = minimizer.findLocalNeighbors(MM_ext, MH_ext, HM_ext, HH_ext, MM, MH, HM, HH, maxNeighbors);
  if(!MPI_rank) {
    cout << "- Done with local neighbor search using a KDTree." << endl;
    cout << "    Maximum number of local neighbors: " << maxNeighbors << ". " << endl;
    cout.flush();}
  MPI_Barrier(comm);
/*
  // Optimize the problem with initial structure
  if(!MPI_rank) {
    cout << "- Optimization for finding adsorption sites: " << endl; cout.flush();}

  err = minimizer.minimizeFreeEntropy(0.0, q_M0, q_H0, sigma_M0, sigma_H0, x0, q_M, q_H, sigma_M, sigma_H, obj_fun);
  if(!MPI_rank){
    cout << "- Done with Max-Ent." << endl;
    cout.flush();
  }
  MPI_Barrier(comm);

  q_M0 = q_M;
  q_H0 = q_H;
  sigma_M0 = sigma_M;
  sigma_H0 = sigma_H;
*/
/*
  err = minimizer.calculateFormationEnergy(x0, gammabd, gamma, f, full_H);
  if(!MPI_rank) {
    cout << "- Done with Formation Energy and Chemical Potential." << endl;
    cout.flush();
  }
  MPI_Barrier(comm);

  findAdsorptionSites(input, gamma, HSubsurf, HAdsorp);
  if(!MPI_rank) {
    cout << "- Done with adsorption site search." << endl;
  cout.flush();}
*/
  //-------------------------------------------------------------------------
  // open file to write average H/M ratio, free energy and potential energy
  // -------------------------------------------------------------------------
  char full_filename_base[128];
  sprintf(full_filename_base, "%s/%s", input.file.foldername, input.file.filename_base);

  char full_fname[128];
  sprintf(full_fname, "%s_history.txt", full_filename_base);
  
  ofstream result_file(full_fname, ios::app);
  if(!result_file){
    cerr << "WARNING! Cannot open the history file!" <<  endl;
  }

  //--------------------------------------------------------------
  // Check start or restart for initialization
  //--------------------------------------------------------------
  if (input.file.restart == 0) {
    if(!MPI_rank) {
      cout << "---------------------------------------" << endl;
      cout << "- Attention: start a new calculation! -" << endl;
      cout << "---------------------------------------" << endl;
      cout.flush();
    }

/*
		err = enforceBoundaryConditions(input, dt, q_M0, q_H0, sigma_M0, sigma_H0, x0, gamma, f, HAdsorp, gammabd, full_H, nH_oct, HH);

		if(!MPI_rank) {
			cout << "- Done with Boundary Condition." << endl;
			cout.flush();
		}
		MPI_Barrier(comm);
*/

    // Optimize the problem with initial condition
    if(!MPI_rank) {
      cout << "- Optimization for initialization: " << endl; cout.flush();}

    err = minimizer.minimizeFreeEntropy(mubd, q_M0, q_H0, sigma_M0, sigma_H0, x0, q_M, q_H, sigma_M, sigma_H, x, obj_fun);
    if(!MPI_rank){
      cout << "- Done with Max-Ent." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    q_M0 = q_M;
    q_H0 = q_H;
    sigma_M0 = sigma_M;
    sigma_H0 = sigma_H;
    x0 = x;

    q_M_ini = q_M;
    q_H_ini = q_H;
    sigma_M_ini = sigma_M;
    sigma_H_ini = sigma_H;
    x_ini = x0;

    if (input.file.output_pot == 1) {
		  err = minimizer.calculateMeanPotentialEnergy(x0, V_M, V_H);
		  if(!MPI_rank) {
			  cout << "- Done with Mean Potential Energy." << endl;
			  cout.flush();
		  }
		  MPI_Barrier(comm);
    }

    if (input.file.output_str == 1) {
      err = minimizer.calculateLocalStress(x0, pi_M, pi_H);
      if(!MPI_rank) {
        cout << "- Done with Local Stresses." << endl;
        cout.flush();
      }
      MPI_Barrier(comm);
    }

    err = calculateMacroResults(x0, xH, obj_fun, free_energy, T, q_M0.size());

    // write time history of macroscropic results to file
    if(!MPI_rank) {
      result_file << "  " << scientific << iFrame << "  " << scientific << 0.0 << "  " << scientific << xH << "  " << scientific << free_energy << endl;
      result_file.flush();
    }
    // write solution to file
    output.output_solution(iFrame++, iTimeStep, q_M0, q_H0, sigma_M0, sigma_H0, x0, gamma, f, V_M, V_H, pi_M, pi_H, full_H, nH_oct, nH_tet, input);
    MPI_Barrier(comm);
    mubd += mubd_step;
    gammabd = mubd / (kB*input.file.T0);
  }
  else if (input.file.restart == 1) {
    if(!MPI_rank) {
      cout << "------------------------------------------" << endl;
      cout << "- Attention: restart an old calculation! -" << endl;
      cout << "------------------------------------------" << endl;
      cout.flush();
    }
    iFrame = input.file.restart_file_num;
    //iTimeStep = input.file.restart_file_num*input.file.output_frequency; 
    //startTimeStep = iTimeStep;
    readResults(input, q_M0, q_H0, sigma_M0, sigma_H0, x0, gamma, f, full_H, t);
    MPI_Barrier(comm);
    
    q_M = q_M0;
    q_H = q_H0;
    sigma_M = sigma_M0;
    sigma_H = sigma_H0;
    x = x0;

    q_M_ini = q_M0;
    q_H_ini = q_H0;
    sigma_M_ini = sigma_M0;
    sigma_H_ini = sigma_H0;
    x_ini = x0;

    if(!MPI_rank) {
      cout << "- Read restart input file: " << input.file.restart_file_num << "." << endl; cout.flush();}

    minimizer.initialize(&input, &MM, &MH, &HM, &HH, q_M0, q_H0, sigma_M0, sigma_H0, x0, &QPp, QWp, &CSembed, &CSrho, &CSpair);
 
    err = minimizer.findLocalNeighbors(MM_ext, MH_ext, HM_ext, HH_ext, MM, MH, HM, HH, maxNeighbors);
    MPI_Barrier(comm);
    if(!MPI_rank) {
        cout << "- Done with local neighbor search using a KDTree:" << endl;
        cout << "    Maximum number of neighbors: " << maxNeighbors << ". " << endl;
        cout.flush();}
  }
  else {
    if(!MPI_rank) {
      cerr << "ERREOR! Wrong Restart in the input file!" <<  endl;
      exit(-1);
    }
  }


  //--------------------------------------------------------------
  // MAIN LOOP 
  //--------------------------------------------------------------
  for (int i=0; i<N_mubd; i++) { 
    // initialization
    q_M0 = q_M_ini;
    q_H0 = q_H_ini;
    sigma_M0 = sigma_M_ini;
    sigma_H0 = sigma_H_ini;
    x0 = x_ini;
    gamma = gamma_ini;

    for (int j=0; j<full_H.size(); j++)
      full_H[j] = 0;

    if(!MPI_rank) {
      cout << endl;
      cout << "------------------------------------" << endl;
      cout << "--           Main Loop            --" << endl;
      cout << "-- Chemical potential: " << mubd << " eV. --" << endl;
      cout << "------------------------------------" << endl;
      cout.flush();
    }

    //--------------------------------------------------------------
    // Minimize grand-canonical free energy (equiv. to free entropy at constant T)  
    // Update q_M, q_H, sigma_M, sigma_H, x
    //--------------------------------------------------------------
    obj_fun = 0.0;
    if(!MPI_rank){
      cout << "# Start with Max-Ent:" << endl;
      cout.flush();
    }
    err = minimizer.minimizeFreeEntropy(mubd, q_M0, q_H0, sigma_M0, sigma_H0, x0, q_M, q_H, sigma_M, sigma_H, x, obj_fun);
    if(!MPI_rank){
      cout << "# Done with Max-Ent." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);
    
    // update old states
    q_M0     = q_M;
    q_H0      = q_H;
    sigma_M0 = sigma_M;
    sigma_H0  = sigma_H;
    x0 = x;

    if(input.file.local_nei == 1) {
      err = minimizer.findLocalNeighbors(MM_ext, MH_ext, HM_ext, HH_ext, MM, MH, HM, HH, maxNeighbors);
          
      if(!MPI_rank) {
        cout << "# Done with local neighbor search using a KDTree. Max. number of neighbors: " << maxNeighbors << ". " << endl;
        cout.flush();}
      MPI_Barrier(comm);
    }

    //--------------------------------------------------------------
    // Calculate formation energy and chemical potential
    //--------------------------------------------------------------
    err = minimizer.calculateFormationEnergy(x0, gammabd, gamma, f, full_H);

    if(!MPI_rank) {
      cout << "# Done with Formation Energy and Chemical Potential." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    // Calculate H mole ratio 
    err = calculateMacroResults(x0, xH, obj_fun, free_energy, T, q_M0.size());
    if(!MPI_rank) {
      cout << "# Done with Macroscopic Parameters. " << endl;
      cout << " - H/M ratio: " << xH << "." << endl;
      cout << " - Chemical potential: " << mubd << "." << endl;
      cout.flush();
    }
    MPI_Barrier(comm);

    //--------------------------------------------------------------
    // Calculate mean potential energy
    //--------------------------------------------------------------
    if (input.file.output_pot == 1) {
      err = minimizer.calculateMeanPotentialEnergy(x0, V_M, V_H);
      if(!MPI_rank) {
        cout << "# Done with Mean Potential Energy." << endl;
        cout.flush();
      }
      MPI_Barrier(comm);
    }

    //--------------------------------------------------------------
    // Calculate local stress
    //--------------------------------------------------------------
    if (input.file.output_str == 1) {
      err = minimizer.calculateLocalStress(x0, pi_M, pi_H);
      if(!MPI_rank) {
        cout << "# Done with Local Stresses." << endl;
        cout.flush();
      }
      MPI_Barrier(comm); 
    }      

    // write time history of macroscropic results to file
    if(!MPI_rank) {
      result_file << "  " << scientific << iFrame << "  " << scientific << mubd << "  " << scientific << xH << "  " << scientific << free_energy << endl;
      result_file.flush();
    }
    // write solution to file
    output.output_solution(iFrame++, iTimeStep, q_M0, q_H0, sigma_M0, sigma_H0, x0, gamma, f, V_M, V_H, pi_M, pi_H, full_H, nH_oct, nH_tet, input);
    MPI_Barrier(comm);

    mubd += mubd_step;
    gammabd = mubd / (kB*input.file.T0);
 
    if(!MPI_rank) {
      cout << "------------------------------------" << endl;
      cout << "--       End of Main Loop         --" << endl;
      cout << "-- Chemical potential: " << mubd << " eV. --" << endl;
      cout << "------------------------------------" << endl;
      cout.flush();
    }
  }

  //--------------------------------------------------------------
  // close result file
  //-------------------------------------------------------------- 
  result_file.close();

  if(!MPI_rank) {
    cout << "==========================================" << endl;
    cout << "   NORMAL TERMINATION (t = " << t << "s) " << endl;
    cout << "==========================================" << endl;
  }
  MPI_Barrier(comm);
  if(!MPI_rank) cout << "Total Computation Time: " << ((double)(clock()-start_time))/CLOCKS_PER_SEC << " sec." << endl;
//  MPI_Finalize(); //will be called in the destructor of Tao

  return 0;
}

void printLogo(MPI_Comm comm)
{
  int MPI_rank = 0;
  MPI_Comm_rank(comm, &MPI_rank);
  if(!MPI_rank) {
    cout << endl;
    cout << " ---------------------------------------------------------------------- " << endl;
    cout << " - Maximum Entropy Calculation for Equilibrium Metal-Hydrogen Systems - " << endl;
    cout << " ---------------------------------------------------------------------- " << endl;
    cout << endl;
    cout.flush();
  }
}
