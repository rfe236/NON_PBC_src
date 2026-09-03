#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
struct Vec3D;
struct Input;

class Output {
  MPI_Comm *comm;
  char filename_base[128];
  char full_filename_base[128];
  ofstream summaryfile;

public:
  Output(MPI_Comm *comm_, Input *input);
  ~Output();
  
  void output_solution(int iFrame, int iTimeStep, vector<Vec3D> &q_M, vector<Vec3D> &q_H,
                     vector<double> &sigma_M, vector<double> &sigma_H, vector<double> &x,
                     vector<double> &gamma, vector<double> &f, vector<double> &V_M, 
                     vector<double> &V_H, vector<vector<double> > &pi_M,
                     vector<vector<double> > &pi_H, 
										 vector<int> &full_H, int &nH_oct, int &nH_tet, Input &input); 

  const string getCurrentDateTime();
};
#endif
