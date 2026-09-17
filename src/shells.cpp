#include <iostream>
#include <vector>
#include "Vector3D.h"
#include <math.h>
using namespace std;

void generateShells(int NC, vector<vector<Int3> > *shells)
{
  int p = 1; //step size

  if(!shells->empty()) {
    cerr << "WARNING: shells not empty!" << endl;
    shells->clear();
  }
  for(int iShell=0; iShell<NC; iShell++) {
    vector<Int3> shell;
    int N = 0;
    while(shell.empty()) { //find the iShell
      int lmax = (int)floor(sqrt((double)p));
      for (int l1=-lmax; l1<=lmax; l1++) 
        for (int l2=-lmax; l2<=lmax; l2++) 
          for (int l3=-lmax; l3<=lmax; l3++)
            if((l1*l1+l2*l2+l3*l3 == p) && ((l1+l2+l3)%2 == 0))
              shell.push_back(Int3(l1,l2,l3));
      if(shell.empty())
        p++;
    }
    shells->push_back(shell);
    p++;
  }
  if(shells->size()!=NC) {
    cerr << "ERROR: size = " << shells->size() << ", NC = " << NC << endl;
    exit(-1);
  }
}

void generateOctahedralInterstitialShells(int NC, vector<vector<Int3> > *shells)
{
  int p = 1; //step size

  if(!shells->empty()) {
    cerr << "ERROR: shells not empty!" << endl;
    shells->clear();
  }
  for(int iShell=0; iShell<NC; iShell++) {
    vector<Int3> shell;
    int N = 0;
    while(shell.empty()) { //find the iShell
      int lmax = (int)floor(sqrt((double)p));
      for (int l1=-lmax; l1<=lmax; l1++)
        for (int l2=-lmax; l2<=lmax; l2++)
          for (int l3=-lmax; l3<=lmax; l3++) 
            if((l1*l1+l2*l2+l3*l3 == p) && ((l1+l2+l3)%2 != 0)) //Note: % behaves differently with 'mod' in MATLAB for negative numbers
              shell.push_back(Int3(l1,l2,l3));
      if(shell.empty())
        p++;
    }
    shells->push_back(shell);
    p++;
  }
  if(shells->size()!=NC) {
    cerr << "ERROR: size = " << shells->size() << ", NC = " << NC << endl;
    exit(-1);
  }
//  for(int i=0; i<shells->size(); i++) 
//    for(int j=0; j<(*shells)[i].size(); j++)
//      cout << "shell " << i << ", id " << j << ": " << (*shells)[i][j][0] << " " << (*shells)[i][j][1] << " " << (*shells)[i][j][2] << endl;
}

void generateTetrahedralInterstitialSites(int N, vector<Vec3D> *sites)
{
	int ext = 5; //number of extended lattice cells in one direction

	if(!sites->empty()) {
		cerr << "ERROR: sites not empty!" << endl;
		sites->clear();
	}
	for (int l1=-(N+ext); l1<=(N+ext); l1++)
		for (int l2=-(N+ext); l2<=(N+ext); l2++)
			for (int l3=-(N+ext); l3<=(N+ext); l3++)
				sites->push_back(Vec3D((double)l1+0.5,(double)l2+0.5,(double)l3+0.5));

//  for(int i=0; i<sites->size(); i++)
//      cout << "tetrahedral site: " << i << ": " << (*sites)[i][0] << " " << (*sites)[i][1] << " " << (*sites)[i][2] << endl;
}
