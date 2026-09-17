#include <iostream>
#include <vector>
#include <math.h>
#include "Vector3D.h"
using namespace std;

int calculateMacroResults(vector<double> &x, double &xH,  
                          double &obj_fun, double &free_energy, double &T, int nM)
{
    xH = 0.0;
    double x_sum0 = 0.0;
    double x_sum1 = 0.0;
    int nH = x.size();

    const double xUb = 1.0-1.0e-16;

    for(int i=0; i<nH; i++) {
      x_sum0 += x[i];
    }

    free_energy = T*obj_fun/(x_sum0+nM);
    xH = x_sum0/nM;
    
		return 0;
}
