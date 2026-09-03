#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <sstream>
#include "Vector3D.h"
#include "input.h"
using namespace std;

// read the results for restart
void readResults(Input &input, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0, vector<double> &sigma_M0,
                 vector<double> &sigma_H0, vector<double> &x0, vector<double> &gamma, 
                 vector<double> &f, vector<int> &full_H, double &t)
{
    double T = input.file.T0;
    const double kB = 8.6173324e-5;

    char full_filename_base[128];    
    sprintf(full_filename_base, "%s/%s", input.file.foldername, input.file.filename_base);

    char full_fname[128];
    
    sprintf(full_fname, "%s.%d.dat", full_filename_base, input.file.restart_file_num);
    ifstream file(full_fname, ios::in);

    if(!file){
        cerr << "ERREOR! The input solution file for restart does not exit!" <<  endl;
        exit(-1);
    }

    // read the solution file for restart line by line
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');
    file.ignore(10000,'\n');

    double tmp1;
    for(int i=0; i<q_M0.size(); i++)
        file >> tmp1 >> tmp1 
             >> q_M0[i][0] 
             >> q_M0[i][1] 
             >> q_M0[i][2] 
             >> sigma_M0[i] 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1 
             >> tmp1 
             >> tmp1 
             >> tmp1 
             >> tmp1;
    for(int i=0; i<q_H0.size(); i++)
        file >> tmp1 >> tmp1 
             >> q_H0[i][0] 
             >> q_H0[i][1] 
             >> q_H0[i][2] 
             >> sigma_H0[i] 
             >> x0[i] 
             >> gamma[i] 
             >> f[i] 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1
             >> tmp1 
             >> tmp1 
             >> tmp1 
             >> full_H[i];

    file.close();

		const double xUb = 1.0-1.0e-16; // The upper bound of hydrogen atomic fraction in the sample
		const double xLb = 1.0e-16; // The lower bound of hydrogen atomic fraction in the sample
		for(int i=0; i<q_H0.size(); i++) { // for numerical stability of log function
			gamma[i] /= (kB*T);
      if(x0[i]<xLb) x0[i] = xLb;
			else if(x0[i]>xUb) x0[i] = xUb;
		}

    // Read time step from the history file
    sprintf(full_fname, "%s_history.txt", full_filename_base);
    vector<string> hist;
    string line;

    ifstream history(full_fname, ios::in);

    if(!history){
      cerr << "WARNING! Cannot open the history file for restarting!" <<  endl;
    } else {
      while(getline(history, line)) {
        hist.push_back(line);
      }
    }

    history.close();
   
    stringstream ss(hist[hist.size()-1]);
    ss >> t;
    hist.clear();

		return; 
}

