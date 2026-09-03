#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "input.h"
#include "parser/Assigner.h"
#include "parser/Dictionary.h"

InputFileData::InputFileData() 
: Fxx_current(1.0), foldername("results/"), filename_base("sol")
{

  // set defaults
  //T0 = 300; 
  //rc = 5.35; 
  //rc_diff = 1.82;
	//N = 6; 
  //abs = 1; 
  //a_M = 3.8;
  //sigma_M = 1000.0; 
  //sigma_H = 1000.0; 
  
  //minimize_frequency = 1;
  //output_frequency = 1;

  //mubd = -2.0; // (Pa)
	//xbd = 1.0;

  //foldername = "";
  //filename_base = "";
 
  //restart = 0;
  //restart_file_num = 0;
}


InputFileData::~InputFileData()
{
//  delete[] foldername;
//  delete[] filename_base;
}


void InputFileData::setup(const char *name)
{
  ClassAssigner *ca = new ClassAssigner(name, 22, 0); //father);

  new ClassStr<InputFileData> (ca, "ResultFolder", this, &InputFileData::foldername);
//  sprintf(foldername,"%s%s", foldername, "/"); // append a slash
  new ClassStr<InputFileData> (ca, "ResultFilePrefix", this, &InputFileData::filename_base);

	new ClassStr<InputFileData> (ca, "EamFile", this, &InputFileData::eam_filename);

  new ClassDouble<InputFileData> (ca, "Temperature", this, &InputFileData::T0);
  new ClassDouble<InputFileData> (ca, "CutOffDistance", this, &InputFileData::rc);

  new ClassInt<InputFileData> (ca, "SampleShapeType", this, &InputFileData::sample_shape);
  new ClassInt<InputFileData> (ca, "NumberOfUnitCells", this, &InputFileData::N);
  
  new ClassDouble<InputFileData> (ca, "InitialHostLatticeConstant", this, &InputFileData::a_M);
  new ClassDouble<InputFileData> (ca, "InitialSigmaHost", this, &InputFileData::sigma_M);
  new ClassDouble<InputFileData> (ca, "InitialSigmaInterstitial", this, &InputFileData::sigma_H);
  new ClassDouble<InputFileData> (ca, "InitialFractionInterstitial", this, &InputFileData::x_ini);

  new ClassInt<InputFileData> (ca, "OutputFrequency", this, &InputFileData::output_frequency);

  new ClassDouble<InputFileData> (ca, "ChemicalPotentialLower", this, &InputFileData::mubd_lb);
  new ClassDouble<InputFileData> (ca, "ChemicalPotentialUpper", this, &InputFileData::mubd_ub);
  new ClassDouble<InputFileData> (ca, "ChemicalPotentialStep", this, &InputFileData::mubd_step);

  new ClassInt<InputFileData> (ca, "FindNeighborAfterOpt", this, &InputFileData::local_nei);
  new ClassInt<InputFileData> (ca, "OutputMeanPotential", this, &InputFileData::output_pot);
  new ClassInt<InputFileData> (ca, "OutputLocalStress", this, &InputFileData::output_str); 
 
  new ClassInt<InputFileData> (ca, "Restart", this, &InputFileData::restart);
  new ClassInt<InputFileData> (ca, "RestartFileNum", this, &InputFileData::restart_file_num); 
  new ClassInt<InputFileData> (ca, "HostSiteNum", this, &InputFileData::nM_restart);
  new ClassInt<InputFileData> (ca, "InterstitialSiteNum", this, &InputFileData::nH_restart); 
}


/*
void InputFileData::initialize(double T0_, double eps_, double rc_, int N_, double xe_, double a_M_,
                  double sigma_M_, double sigma_H_, double dt_, double t_final_,
                  int output_frequency_, char* foldername_, char *filename_base_)
{
  T0 = T0_;  eps = eps_;  rc = rc_;  N = N_;  xe = xe_;  a_M = a_M_;
  sigma_M = sigma_M_;  sigma_H = sigma_H_;  dt = dt_;  t_final = t_final_;
  output_frequency = output_frequency_;
  sprintf(foldername, "%s", foldername_);
  sprintf(filename_base, "%s", filename_base_);
}
*/

//-----------------------------------------------------

void Input::readCmdLine(int argc, char** argv)
{
  if(argc==1) {
    fprintf(stderr,"ERROR: Input file not provided!\n");
    exit(-1);
  }
  cmdFileName = argv[1];
}


void Input::readCmdFile()
{
  extern FILE *yyCmdfin;
  extern int yyCmdfparse();

  setupCmdFileVariables();
//  cmdFilePtr = freopen(cmdFileName, "r", stdin);
  yyCmdfin = cmdFilePtr = fopen(cmdFileName, "r");

  if (!cmdFilePtr) {
    fprintf(stderr,"*** Error: could not open \'%s\'\n", cmdFileName);
    exit(-1);
  }

  int error = yyCmdfparse();
  if (error) {
    fprintf(stderr,"*** Error: command file contained parsing errors\n");
    exit(error);
  }
  fclose(cmdFilePtr);
}


void Input::setupCmdFileVariables()
{
  file.setup("DMDInputs");
}
