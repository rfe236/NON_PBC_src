#include <vector>
#include <iostream>
#include "Vector3D.h"
#include "KDTree.h"
#include "input.h"
#include <mpi.h>
using namespace std;

int findNeighbors(double a_M, double rc, double Lx, double Ly, double Lz, int pbc, vector<Vec3D> &q_M0, vector<Vec3D> &q_H0,
                  vector<vector<int> > &MM, vector<vector<int> > &MH,
                  vector<vector<int> > &HM, vector<vector<int> > &HH)

 {
    if(!(MM.empty() && MH.empty() && HM.empty() && HH.empty())) {
        cout << "WARNING: Some neighbor list(s) is/are not empty. Clearing all the lists." << endl;
        MM.clear();
        MH.clear();
        HM.clear();
        HH.clear();
    }

    int nM = q_M0.size();
    int nH = q_H0.size();
    int Natoms = nM + nH;
 vector<int> tmp;
    tmp.reserve(1000);

    MM.assign(nM, tmp); // M neighbors of M
    MH.assign(nM, tmp); // H neighbors of M
    HM.assign(nH, tmp); // M neighbors of H
    HH.assign(nH, tmp); // H neighbors of H

  //  double lo = -0.5 * L;
  //  double hi =  0.5 * L;
 double loX=-0.5*Lx, hiX=0.5*Lx;
double loY=-0.5*Ly, hiY=0.5*Ly;
double loZ=-0.5*Lz, hiZ=0.5*Lz;
 vector<PointIn3D> pts;
    pts.reserve(Natoms * 27);

    for(int a=0; a<Natoms; a++) {

        Vec3D p = (a < nM) ? q_M0[a] : q_H0[a-nM];
 pts.push_back(PointIn3D(a, p));

        if(pbc) {
            for(int sx=-1; sx<=1; sx++) {
                for(int sy=-1; sy<=1; sy++) {
                    for(int sz=-1; sz<=1; sz++) {

                        if(sx==0 && sy==0 && sz==0) continue;

 Vec3D g(p[0] + sx*Lx, p[1] + sy*Ly, p[2] + sz*Lz);
                  //      bool nearBox =
                  //          (g[0] > lo-rc && g[0] < hi+rc) &&
                  //          (g[1] > lo-rc && g[1] < hi+rc) &&
                  //          (g[2] > lo-rc && g[2] < hi+rc);

                  //      bool insideBox =
                  //          (g[0] > lo && g[0] < hi) &&
                 //           (g[1] > lo && g[1] < hi) &&
                 //           (g[2] > lo && g[2] < hi);
  bool nearBox   = (g[0]>loX-rc && g[0]<hiX+rc) && (g[1]>loY-rc && g[1]<hiY+rc) && (g[2]>loZ-rc && g[2]<hiZ+rc);
bool insideBox = (g[0]>loX && g[0]<hiX) && (g[1]>loY && g[1]<hiY) && (g[2]>loZ && g[2]<hiZ);
  if(!nearBox || insideBox) continue;
 pts.push_back(PointIn3D(a, g));
                    }
                }
            }
        }
    }
 int Npts = pts.size();

    PointIn3D *atoms = new PointIn3D[Npts];
    for(int i=0; i<Npts; i++) {
        atoms[i] = pts[i];
    }

    KDTree<PointIn3D> atomTree(Npts, atoms);

    int maxNeib = 0;
    int maxNei = 1000;
    PointIn3D candidates[maxNei];
    double qLocal[3];

    for(int i=0; i<Natoms; i++) {

        Vec3D qi = (i < nM) ? q_M0[i] : q_H0[i-nM];

        qLocal[0] = qi[0];
        qLocal[1] = qi[1];
        qLocal[2] = qi[2];

        int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);

        if(nFound > maxNei) {
            cerr << "ERROR: found " << nFound
                 << " neighbors, allocated space = " << maxNei
                 << ". Increase maxNei!" << endl;
            exit(-1);
        }
 for(int j=0; j<nFound; j++) {

            int nei = candidates[j].pid(); // real atom id

            if(i == nei) continue;

            double dx = candidates[j].val(0) - qLocal[0];
            double dy = candidates[j].val(1) - qLocal[1];
            double dz = candidates[j].val(2) - qLocal[2];

            if(dx*dx + dy*dy + dz*dz > rc*rc) continue;

            if(i < nM) { // center atom is M

                if(nei < nM) { // neighbor is M

                    bool duplicate = false;
                    for(int k=0; k<MM[i].size(); k++) {
                        if(MM[i][k] == nei) {
                            duplicate = true;
                            break;
                        }
                    }

                    if(!duplicate) MM[i].push_back(nei);

                } else { // neighbor is H

                    int hnei = nei - nM;

                    bool duplicate = false;
                    for(int k=0; k<MH[i].size(); k++) {
                        if(MH[i][k] == hnei) {
                            duplicate = true;
                            break;
                        }
                    }

                    if(!duplicate) MH[i].push_back(hnei);
                }

            } else { // center atom is H

                int ih = i - nM;

                if(nei < nM) { // neighbor is M

                    bool duplicate = false;
                    for(int k=0; k<HM[ih].size(); k++) {
                        if(HM[ih][k] == nei) {
                            duplicate = true;
                            break;
                        }
                    }

                    if(!duplicate) HM[ih].push_back(nei);

                } else { // neighbor is H

                    int hnei = nei - nM;

                    bool duplicate = false;
                    for(int k=0; k<HH[ih].size(); k++) {
                        if(HH[ih][k] == hnei) {
                            duplicate = true;
                            break;
                        }
                    }

                    if(!duplicate) HH[ih].push_back(hnei);
                }
            }
        }

        int nNeib = 0;

        if(i < nM) {
            nNeib = MM[i].size() + MH[i].size();
        } else {
            int ih = i - nM;
            nNeib = HM[ih].size() + HH[ih].size();
        }

        if(nNeib > maxNeib) maxNeib = nNeib;
    }

    delete[] atoms;

    int cmin = 1000000000;
    int cmax = 0;

    for(int i=0; i<nM; i++) {
        int c = MM[i].size();
        if(c < cmin) cmin = c;
        if(c > cmax) cmax = c;
    }

int rank = 0;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if (rank == 0) {
    cout << "EXT M-M pool coordination: min="
         << cmin << " max=" << cmax << endl;
} 
//  vector<int> tmp1; //a zero vector
//  tmp1.assign(3, 0);

  // Store all the sites in a KD-Tree with K = 3.
//  int Natoms = q_M0.size() + q_H0.size();
//  PointIn3D *atoms = new PointIn3D[Natoms];
//  for(int i=0; i<q_M0.size(); i++)
//    atoms[i] = PointIn3D(i, q_M0[i]);
//  for(int i=0; i<q_H0.size(); i++)
//    atoms[q_M0.size() + i] = PointIn3D(q_M0.size() + i, q_H0[i]);
//  KDTree<PointIn3D> atomTree(Natoms, atoms); //Note: atoms are re-ordered.

  // Find and store neighbors
//  int maxNeib = 0;
//  int maxNei = 1000; //max number of neighbors for each atom. 
//  PointIn3D candidates[maxNei];
//  double qLocal[3];
//  for(int i=0; i<Natoms; i++) {
//    int nNeib = 0; 
//    for(int j=0; j<3; j++)
 //     qLocal[j] = atoms[i].val(j); // position of i (could be either M or H)
    // find neighbors of i within rc --> populates "candidates"
//    int nFound = atomTree.findCandidatesWithin(qLocal, candidates, maxNei, rc);
//    // debug
//    if(nFound>maxNei) {
//      cerr << "ERROR: found " << nFound << " neighbors, allocated space = " << maxNei 
//           << ". Increase maxNei! " << endl;
//      exit(-1);
//    } 
/*
    for(int j=0; j<nFound; j++) {
      Vec3D thisone(qLocal[0],qLocal[1],qLocal[2]);
      Vec3D neib(candidates[j].val(0), candidates[j].val(1), candidates[j].val(2));
      Vec3D neib2;
      if (candidates[j].pid()<q_M0.size())
        neib2 = q_M0[candidates[j].pid()];
      else
        neib2 = q_H0[candidates[j].pid() - q_M0.size()];

      Vec3D dif = thisone - neib;
      Vec3D dif0 = neib - neib2;
      cout << "i = " << i << ", j = " << j << ", d = " << max(fabs(dif[0]), max(fabs(dif[1]), fabs(dif[2]))) << ", 0 = " << dif0.norm() << endl;
    }
*/

    // populates MM, MH, HM, HH  
 /*
    for(int j=0; j<nFound; j++) {
      int nei = candidates[j].pid(); //one neighbor
      //filter the same one and the outside ones -- KDTree works with "pseudo distance"
      if(atoms[i].pid()==nei) 
        continue;
      if((candidates[j].val(0)-qLocal[0])*(candidates[j].val(0)-qLocal[0]) +
         (candidates[j].val(1)-qLocal[1])*(candidates[j].val(1)-qLocal[1]) +
         (candidates[j].val(2)-qLocal[2])*(candidates[j].val(2)-qLocal[2]) > rc*rc)
        continue;
      nNeib++;
      if(atoms[i].pid()<q_M0.size()) {//working on a M site
        if(nei<q_M0.size()) {// this neighbor is M
          MM[atoms[i].pid()].push_back(nei);
        } else { // this neighbor is H
          MH[atoms[i].pid()].push_back(nei-q_M0.size());
        }
      } else {// working on an H site
        if(nei<q_M0.size()) {//this neighbor is M
          HM[atoms[i].pid()-q_M0.size()].push_back(nei);
        } else {//this neighbor is H
          HH[atoms[i].pid()-q_M0.size()].push_back(nei-q_M0.size());
        }
      }
    }

    if(nNeib>maxNeib) maxNeib = nNeib;
  }
/*
//  for(int i=0; i<q_M0.size(); i++) {
//    cout << "M[" << i << "]:" << q_M0[i][0] << " " << q_M0[i][1] << " " << q_M0[i][2] << endl;
//    for(int j=0; j<MM[i].size(); j++) {
//      Vec3D delta = q_M0[i] - q_M0[MM[i][j]];
//      double dist = delta.norm();
//      cout << "  Nei " << j << ": id = " << MM[i][j] << ", " << q_M0[MM[i][j]][0] << " " <<  q_M0[MM[i][j]][1] << " " << q_M0[MM[i][j]][2] << ", d = " << dist << endl;
//    }
//  }
*/

 

  return maxNeib;
}

