#ifndef _SAL_H
#define _SAL_H

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#define HEADER_SIZE 16

#define MIN_FAC 0.99999
#define MAX_FAC 1.00001

class SAL {
  private:
    char last_iq;
    double Lambda4sq;
    double ipol2(double x, double t, int pn);
    double ipol2d(double x, double t, int pn);
    double** fi[7];
    double** dfi[7];
    double** fi2[7];  //--- second derivatives of fi
    double** dfi2[7]; //--- second derivatives of dfi
    double* X;
    double* T;

    double BetheHeitler(double x, double QQ, double msq);
    double PSSuppr(double x, double QQ, double msq);
    int InLoc(double X[], int n, double x);
    void TriDiagSolve(double diag[], double above[], double below[], double b[], int n);
    double* InitSpline(double x[], double y[], int n, double yp0);
    double SplineValue(double x[], double y[], double M[], int npt, double xi);
    
    //==========================================
    void Aloc(double** a[7]) {
      int pn,iq;
      try {
        for(pn=0; pn<=Nf; pn++){
          a[pn] = new double*[nq];
          for(iq = 0; iq < nq; iq++) a[pn][iq] = new double[nx];
        }
      }
      catch(...) { cerr << "Could not allocate. Bye ..." <<endl; exit(1); }
    }

    //==========================================
    void Dealloc(double** a[7]) {
      int pn,iq;
      for(pn=0; pn<=Nf; pn++){
        for(iq = 0; iq < nq; iq++) delete[] a[pn][iq];
        delete[] a[pn];
      }
    }

    //==========================================
    double evSuppr(double QQ, int flav) {
      if(QQ <= m2[flav]) return 1;
      double tau = log(log(QQ/Lambda4sq)/
             log(m2[flav]/Lambda4sq));
      return tau > 1 ? 0 : 1 - tau;
    }

  public:
    unsigned char ver,nx,nq, Nf;
    unsigned char* Flav;
    double x_min, x_max;
    double QQ_min, QQ_max;
    double eq2[7];
    double m2[7];

    SAL(char* iname = NULL);
    
    ~SAL() {
      Dealloc(fi);
      Dealloc(dfi);
      Dealloc(fi2);
      Dealloc(dfi2);
    }
    
    
    void PDF(double x, double QQ, double f[]);
    double F2(double x, double QQ, int hflav=0);

};

class TxGrid {
  private:
    int Nx;
    double aNx, xlo, dx2;
  public:
    TxGrid(int nx) {
      Nx = nx;
      aNx = 22;
      xlo = 1e-5;
      double dx1 = exp(0.5*aNx);
      dx2 = (1-xlo)*(dx1+1)/(dx1-1);
    }
    double x(int ix) {
      return dx2/(1+exp(0.5*aNx*(1-(2.0*ix)/Nx))) + (1+xlo - dx2)/2;
    }
};

//-------  non-class functions  ---------------------
//---  defined in sal.cpp

void SALPDF(double x, double QQ, double f[]);
double SALF2(double x, double QQ);

#endif
