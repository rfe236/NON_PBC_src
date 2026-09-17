#ifndef _INPUT_H_
#define _INPUT_H_
#include <stdio.h>

struct InputFileData
{
  // constants
  double T0;  // temperature (K)
  double rc;  // cutoff distance (Ang)
double Fxx_current;
  // computational domain
  int sample_shape;
  int N; //N*N*N = number of unit cells;

  // initial condition
  double a_M; //M lattice constant (Ang)
  double sigma_M; //M frequency (Ang^-2)
  double sigma_H; //H frequency (Ang^-2)
  double x_ini; // H fraction

  // output
  int output_frequency;
  const char* foldername;
  const char* filename_base;

	// eam 
	const char* eam_filename;

  // boundary condition
  double mubd_lb;
  double mubd_ub;
  double mubd_step;

  // find local neighbors
  int local_nei;
  
  // output mean potential
  int output_pot;

  // output local stress
  int output_str;

  // continue setup
  int restart;
  int restart_file_num;

  // number of sites in restart file
  int nM_restart;
  int nH_restart;

  InputFileData();
  ~InputFileData();

  void setup(const char *);

  // the following function is obsolete. 
/*
  void initialize(double T0_, double eps_, double rc_, int N_, double xe_, double a_M_,
                  double sigma_M_, double sigma_H_, double dt_, double t_final_,
                  int output_frequency_, char* foldername_, char *filename_base_);
*/
};


class Input
{
  char *cmdFileName;
  FILE *cmdFilePtr;

public:

  InputFileData  file;
  
public:

  Input() {}
  ~Input() {}

  void readCmdLine(int, char**);
  void readCmdFile();
  void setupCmdFileVariables();
  
};
#endif
