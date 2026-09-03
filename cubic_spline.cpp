#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "cubic_spline.h"

using namespace std;

static double distance( double* q1, double* q2 )
{
  double dist = 0.0;

  for (int i=0; i<3; i++)
    dist += (q1[i] - q2[i]) * (q1[i] - q2[i]);

  return sqrt(dist);
}

void init_spline( CubicSpline* cs, int n, double dx, double *y )
{
  cs->n  = n;
  cs->dx = dx;
	cs->x  = new double[n];
	cs->a  = new double[n];
  cs->b  = new double[n];
  cs->c  = new double[n];
  cs->d  = new double[n];
	
	//printf("Test: %s\n", y[0]);

	// Generate the coefficients of a cubic Hermite spline.
	// n is the number of coefficients, so the number of data, i.e., the size of y, should be n+1.
	int i;
	double *y1;  // used to store first derivatives
	y1 = new double[n+1]; 

	// calculate first derivative at the first point by forward difference
	y1[0] = (y[1] - y[0]) / dx;
	// calculate first derivative at the last point by backward difference
	y1[n] = (y[n] - y[n-1]) / dx;
	// calculate first derivatives at other points by central difference
	for (i=1; i<n; i++) {
		y1[i] = (y[i+1] - y[i-1]) / (2.0*dx);
	}

	for (i=0; i<n; i++) {
		cs->x[i] = i * dx; // Note: the variable starts from zero
		cs->a[i] = (2.0*y[i] + dx*y1[i] - 2.0*y[i+1] + dx*y1[i+1]) / (dx*dx*dx);
		cs->b[i] = (-3.0*y[i] - 2.0*dx*y1[i] + 3.0*y[i+1] - dx*y1[i+1]) / (dx*dx);
		cs->c[i] = y1[i];
		cs->d[i] = y[i];
	}

	delete [] y1;
}

void destroy_spline( CubicSpline* cs )
{
  delete [] cs->x;
  delete [] cs->a;
  delete [] cs->b;
  delete [] cs->c;
  delete [] cs->d;
}

double cubic_spline( CubicSpline* cs, double x ) // XS: cubic spline function
{
//    int j;
//    int i = cs->n - 1;

  //if (isinf(x)) printf("cs 53 x %f\n",x);

  if (fabs(x) < cs->x[0]) {
    printf("Error: Atoms very close. Potential parameters are not defined for this interatomic distance: %lf", x);// %lf %lf %lf %lf %lf \n", x, cs->x[0], cs->d[0], cs->c[0], cs->b[0], cs->a[0]);
    exit(0);
  }

//    for(j=1; j<cs->n; j++ ){
//	
//        if( x < cs->x[j] ){
//	    i = j - 1;
//	    break;
//	}
//
//    }
    
  int i = static_cast<int>(x/cs->dx);
  i = min(i, cs->n-1);   
  //if (isinf(x)) printf("cs 71 i %i x %f csdx %f csn %i\n",i,x,cs->dx,cs->n-1);  
  double dx = x - cs->x[i];
  dx = min(dx,cs->dx);
  double dx2 = dx*dx;

  return cs->d[i] + ( cs->c[i] + cs->b[i]*dx + cs->a[i]*dx2 )*dx;
}

double d_cubic_spline( CubicSpline* cs, double x ) // XS: first derivative of cubic spline function
{
//    int j;
//    int i = cs->n - 1;
    
  if (fabs(x) < cs->x[0]) {
    printf("Error: Atoms very close. Potential parameters are not defined for this interatomic distance: %lf \n",x);
    exit(0);
  }

//    for( j = 1; j < cs->n; j++ ){
//	
//        if( x <cs->x[j] ){
//            i = j - 1;
//            break;
//	}
//    }    
        
  int i = static_cast<int>(x/cs->dx);
  i = min(i, cs->n-1);   
        
  double dx = x - cs->x[i];
  dx = min(dx,cs->dx);

  return cs->c[i] + (2.0*cs->b[i] + 3.0*cs->a[i]*dx)*dx;
}

double d2_cubic_spline( CubicSpline* cs, double x ) // XS: second derivative of cubic spline function
{
//    int j;
//    int i = cs->n - 1;
    
  if (fabs(x) < cs->x[0]) {
    printf("Error: Atoms very close. Potential parameters are not defined for this interatomic distance: %lf \n", x);
    exit(0);
  }

//    for( j = 1; j < cs->n; j++ ){
//	
//        if( x <cs->x[j] ){
//            i = j - 1;
//            break;
//	}
//    }
    
  int i = static_cast<int>(x/cs->dx);
  i = min(i, cs->n-1); 
    
  double dx = x - cs->x[i];
  dx = min(dx,cs->dx);

  return 2.0*cs->b[i] + 6.0*cs->a[i]*dx;
}

double mf_cubic_spline( CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b)
{
  int N_a = 2, k;
  double r= sqrt(3.0*(double)N_a)*98.227002603;//0.0101805; // (1/eV)^(1/2)*(gr/mol)^(1/2)*(1/ps)= 0.010180505/A 
  double W_k = 1.0/(6.0*(double)N_a);
  double desp, desp2, rijm, mf=0.0;
  double qm[3];  
           
  // for atom_i    
  desp  = r/( sqrt(b*mass_i)*w_i );
  desp2 = 2.0*desp; 
  
  for (k=0; k<3; k++) {
    qm[0] = qi[0];
    qm[1] = qi[1];
    qm[2] = qi[2];        

    // positive        
    qm[k]+= desp;
    rijm = distance (qm, qj);               
    mf += cubic_spline( cs, rijm );        

    // negative
    qm[k]-= desp2;
    rijm = distance (qm, qj);               
    mf += cubic_spline( cs, rijm );         
  }
    
  // for atom_j
  desp =r/( sqrt(b*mass_j)*w_j );
  desp2 = 2.0*desp;
  
  for (k=0; k<3; k++) {
    qm[0] = qj[0];
    qm[1] = qj[1];
    qm[2] = qj[2];        

    // positive        
    qm[k]+= desp;
    rijm = distance (qi, qm);        
    mf += cubic_spline( cs, rijm );          

    // negative
    qm[k]-= desp2;
    rijm  = distance (qi, qm);        
    mf   += cubic_spline( cs, rijm );      
  }   
    
  return mf*W_k;
}

void mf_d_cubic_spline( double mf[3], double *mfwi, double *mfwj, CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b)
{
  int N_a = 2, k;
  double r= sqrt(3.0*(double)N_a)*98.227002603;//0.0101805 = 1/98.227002603; // (1/eV)^(1/2)*(gr/mol)^(1/2)*(1/ps)= 0.010180505/A //  r = x * sqrt(2) = sqrt(3*Na/2) * sqrt(2) = sqrt(3*Na)
  double W_k = 1.0/(6.0*(double)N_a); // Wk * (1/pi^(1/2))^(3*Na) = ( (1/2/3/Na) * pi^(3*Na/2) ) * (1/pi^(1/2))^(3*Na) = 1/(6*Na)
  double X = sqrt(3*N_a/2); // sqrt(3*Na/2)
  double desp, desp2, rijm, value;
  double qm[3];
    
  //Initialization
  mf[0] = 0.0; 
  mf[1] = 0.0; 
  mf[2] = 0.0; 
  *mfwi = 0.0; 
  *mfwj = 0.0;
           
  // for atom_i    
  desp  = r/( sqrt(b*mass_i)*w_i );
  desp2 = 2.0*desp; 
    
  for (k=0; k<3; k++) {
    qm[0] = qi[0];
    qm[1] = qi[1];
    qm[2] = qi[2];        

    // positive        
    qm[k]+= desp;
    rijm  = distance (qm, qj);
    value = d_cubic_spline( cs, rijm )/rijm;
    mf[0]   += value*(qj[0]-qm[0]);
    mf[1]   += value*(qj[1]-qm[1]);
    mf[2]   += value*(qj[2]-qm[2]);
    (*mfwi) += value*(qj[k]-qm[k])*X;

    // negative
    qm[k]-= desp2;
    rijm  = distance (qm, qj);
    value = d_cubic_spline( cs, rijm )/rijm;
    mf[0]   += value*(qj[0]-qm[0]);
    mf[1]   += value*(qj[1]-qm[1]);
    mf[2]   += value*(qj[2]-qm[2]);
    (*mfwi) -= value*(qj[k]-qm[k])*X; //(qj-qm)*(-X)
  }
    
  // for atom_j
  desp  =r/( sqrt(b*mass_j)*w_j );
  desp2 = 2.0*desp;
  
  for (k=0; k<3; k++) {
    qm[0] = qj[0]; 
    qm[1] = qj[1]; 
    qm[2] = qj[2];        

    // positive        
    qm[k]+= desp;
    rijm  = distance (qi, qm);
    value = d_cubic_spline( cs, rijm )/rijm;       
    mf[0]   += value*(qm[0]-qi[0]);
    mf[1]   += value*(qm[1]-qi[1]);
    mf[2]   += value*(qm[2]-qi[2]);
    (*mfwj) -= value*(qm[k]-qi[k])*X;

    // negative
    qm[k]-= desp2;
    rijm  = distance (qi, qm);
    value = d_cubic_spline( cs, rijm )/rijm;        
    mf[0]   += value*(qm[0]-qi[0]);
    mf[1]   += value*(qm[1]-qi[1]);
    mf[2]   += value*(qm[2]-qi[2]);
    (*mfwj) += value*(qm[k]-qi[k])*X; //(qm-qi)*(-X)
  }   
    
  mf[0]   *= W_k;
  mf[1]   *= W_k;
  mf[2]   *= W_k;
  (*mfwi) *= W_k;
  (*mfwj) *= W_k;
}

void mf_d_cubic_spline( double mf[6], CubicSpline* cs, double qi[3], double qj[3], double mass_i, double mass_j, double w_i, double w_j, double b)
{
  int N_a = 2, i, k;
  double r= sqrt(3.0*(double)N_a)*98.227002603;//0.0101805 = 1/98.227002603; // (1/eV)^(1/2)*(gr/mol)^(1/2)*(1/ps)= 0.010180505/A //  r = x * sqrt(2) = sqrt(3*Na/2) * sqrt(2) = sqrt(3*Na)
  double W_k = 1.0/(6.0*(double)N_a); // Wk * (1/pi^(1/2))^(3*Na) = ( (1/2/3/Na) * pi^(3*Na/2) ) * (1/pi^(1/2))^(3*Na) = 1/(6*Na)
//  double X = sqrt(3*N_a/2); // sqrt(3*Na/2)
  double desp, desp2, rijm, value;
  double qm[3];
    
  //Initialization
  for (i=0; i<6; i++)
    mf[i] = 0.0;
           
  // for atom_i    
  desp  = r/( sqrt(b*mass_i)*w_i );
  desp2 = 2.0*desp; 
    
  for (k=0; k<3; k++) {
    qm[0] = qi[0];
    qm[1] = qi[1];
    qm[2] = qi[2];        

    // positive        
    qm[k]+= desp;
    rijm  = distance (qm, qj);
    value = d_cubic_spline( cs, rijm )/rijm;
    mf[0] += (qj[0]-qm[0])*value*(qj[0]-qm[0]);//xx
    mf[1] += (qj[1]-qm[1])*value*(qj[1]-qm[1]);//yy
    mf[2] += (qj[2]-qm[2])*value*(qj[2]-qm[2]);//zz
    mf[3] += (qj[0]-qm[0])*value*(qj[1]-qm[1]);//xy
    mf[4] += (qj[0]-qm[0])*value*(qj[2]-qm[2]);//xz
    mf[5] += (qj[1]-qm[1])*value*(qj[2]-qm[2]);//yz

    // negative
    qm[k]-= desp2;
    rijm  = distance (qm, qj);
    value = d_cubic_spline( cs, rijm )/rijm;
    mf[0] += (qj[0]-qm[0])*value*(qj[0]-qm[0]);//xx
    mf[1] += (qj[1]-qm[1])*value*(qj[1]-qm[1]);//yy
    mf[2] += (qj[2]-qm[2])*value*(qj[2]-qm[2]);//zz
    mf[3] += (qj[0]-qm[0])*value*(qj[1]-qm[1]);//xy
    mf[4] += (qj[0]-qm[0])*value*(qj[2]-qm[2]);//xz
    mf[5] += (qj[1]-qm[1])*value*(qj[2]-qm[2]);//yz            
  }
    
  // for atom_j
  desp  = r/( sqrt(b*mass_j)*w_j );
  desp2 = 2.0*desp;
  
  for (k=0; k<3; k++) {
    qm[0] = qj[0];
    qm[1] = qj[1];
    qm[2] = qj[2];        

    // positive        
    qm[k]+= desp;
    rijm  = distance (qi, qm);
    value = d_cubic_spline( cs, rijm )/rijm;       
    mf[0] += (qm[0]-qi[0])*value*(qm[0]-qi[0]);//xx
    mf[1] += (qm[1]-qi[1])*value*(qm[1]-qi[1]);//yy
    mf[2] += (qm[2]-qi[2])*value*(qm[2]-qi[2]);//zz
    mf[3] += (qm[0]-qi[0])*value*(qm[1]-qi[1]);//xy
    mf[4] += (qm[0]-qi[0])*value*(qm[2]-qi[2]);//xz
    mf[5] += (qm[1]-qi[1])*value*(qm[2]-qi[2]);//yz

    // negative
    qm[k]-= desp2;
    rijm  = distance (qi, qm);
    value = d_cubic_spline( cs, rijm )/rijm;        
    mf[0] += (qm[0]-qi[0])*value*(qm[0]-qi[0]);//xx
    mf[1] += (qm[1]-qi[1])*value*(qm[1]-qi[1]);//yy
    mf[2] += (qm[2]-qi[2])*value*(qm[2]-qi[2]);//zz
    mf[3] += (qm[0]-qi[0])*value*(qm[1]-qi[1]);//xy
    mf[4] += (qm[0]-qi[0])*value*(qm[2]-qi[2]);//xz
    mf[5] += (qm[1]-qi[1])*value*(qm[2]-qi[2]);//yz
  }   
    
  for (i=0; i<6; i++)
    mf[i] *= W_k;
}
