#define MAJOR_SAL_VERSION 1
#define MINOR_SAL_VERSION 3
#define SAL_FILE_VERSION 1

/*
  Format of the SAL.bin file:
  ---------------------------
  Header (char[HEADER_SIZE], cf. sal.h)
  Lambda4, mc, mb, mt, x_min, QQ_min, QQ_max (double)
  nx, nq   (char)
  x_values (nx*double)
  t_values (nq*double)
  active_flavors (nq*char)
  sfi[pn][iq][ix]  (double)
  sdfi[pn][iq][ix] (double)
  derivatives at x[0]:
  sfi1[pn][iq]  (double)
  sdfi1[pn][iq] (double)
*/

#include <cstring>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include "sal.h"

using namespace std;

enum {gluon, d_quark, u_quark, s_quark, c_quark, b_quark, t_quark};

//#define BIN_INPUT

#ifdef BIN_INPUT
  char* DEF_INPUT="SAL.bin";

  //===================================
  static int GetByte(FILE* f, unsigned char *v, int n=1) {
    return n == fread(v, 1, n, f);
  }

  //===================================
  static int GetDouble(FILE* f, double *v, int n=1) {
    return n == fread(v, sizeof(double), n, f);
  }
#else
  // TODO: make wrapper for the path to the SAL.dat file
  char* DEF_INPUT="/home/pmeinzinger/git/sherpa/PDF/SAL/SAL.dat";

  //===================================
  static int GetByte(FILE* f, unsigned char *v, int n=1) {
    unsigned short int u;
    while(n--) {
      if(1 != fscanf(f, "%hu", &u)) return 0;
      *v++ = u & 0xFF;
    }
    return 1;
  }

  //===================================
  static int GetDouble(FILE* f, double *v, int n=1) {
    while(n--) {
      if(1 != fscanf(f, "%lf", v++)) return 0;
    }
    return 1;
  }
#endif


//===================================
SAL::SAL(char* iname /* = NULL */) {
  unsigned char hdr[HEADER_SIZE];
  char Label[32];
  Nf = 6;
  FILE *infil;
  
  if(!iname) iname = DEF_INPUT;
  infil = fopen(iname, "rb");
  if(!infil) {
    cerr << "Cannot open '"<< iname <<"' for reading" <<endl;
    exit(1);
  }

#ifdef BIN_INPUT
  GetByte(infil, hdr, sizeof(hdr));
#else
  GetByte(infil, hdr, 2);
  fscanf(infil, "%s", Label);
#endif
  ver = hdr[0];
  double v[7];
  int i=0;
  GetDouble(infil, v, 7);
  Lambda4sq = pow(v[i++], 2);
  memset(m2, 0, sizeof(m2));
  m2[c_quark] = pow(v[i++], 2);
  m2[b_quark] = pow(v[i++], 2);
  m2[t_quark] = pow(v[i++], 2);
  x_min = v[i++] *MIN_FAC;
  x_max = 0.9999;
  QQ_min = v[i++] *MIN_FAC;
  QQ_max = v[i++] *MAX_FAC;

  GetByte(infil, &nx);
  GetByte(infil, &nq);
  X = new double[nx];
  T = new double[nq];
  GetDouble(infil, X, nx);
  GetDouble(infil, T, nq);
  Flav = new unsigned char[nq];
  GetByte(infil, Flav, nq);
  int iq, pn;
  Aloc(fi);
  for(iq=0; iq < nq; iq++) for(pn=0; pn <= Flav[iq]; pn++)
    GetDouble(infil, fi[pn][iq], nx);
  Aloc(dfi);
  for(iq=0; iq < nq; iq++) for(pn=0; pn <= Flav[iq]; pn++)
    GetDouble(infil, dfi[pn][iq], nx);
  double f1;

  for(pn=0; pn <= Nf; pn++) {
    eq2[pn] = ((pn & 1) ? 1 : 4)/9.0;
    fi2[pn] = new double*[nq];
    dfi2[pn] = new double*[nq];
  }
  eq2[0] = 0;

  for(iq=0; iq < nq; iq++) for(pn=0; pn <= Flav[iq]; pn++) {
    GetDouble(infil, &f1);
    fi2[pn][iq] = InitSpline(X, fi[pn][iq], nx, f1);
  }

  for(iq=0; iq < nq; iq++) for(pn=0; pn <= Flav[iq]; pn++) {
    GetDouble(infil, &f1);
    dfi2[pn][iq] = InitSpline(X, dfi[pn][iq], nx, f1);
  }

  fclose(infil);
}

//====================================================
void SAL::PDF(double x, double QQ, double f[]) {
  double fac0 = 1/(2*M_PI*x*x);
  double fac = fac0/((1-x)*(1-x));
  int pn;
  memset(f, 0, sizeof(double)*7);
  if(x >= x_max) return;
  if(x < x_min) {
    fprintf(stderr, "x = %g below the minimal allowed value %.4g.", x, x_min);
    exit(2);
  }
  if(QQ < QQ_min || QQ > QQ_max) {
    fprintf(stderr, "Q^2 = %g outside the allowed range [%.4g, %.4g]."
      , QQ, QQ_min, QQ_max);
    exit(2);
  }
  double t = log(QQ);
  last_iq = InLoc(T, nq, t);
  //cout << "last_iq = " << (int)last_iq << endl;
  int Aflav = Flav[last_iq];
  //cout << "Aflav " << Aflav << endl;
  f[0] = ipol2(x, t, 0) *fac0;
  for(pn=1; pn <= Aflav; pn++) {
  //cout << "pn " << pn << endl;
    f[pn] = ipol2(x, t, pn) *fac;
}
  for(pn=0; pn <= Aflav; pn++) if(f[pn] < 0 ) f[pn] = 0;
  for(pn=Aflav+1; pn <= Nf; pn++) f[pn] = 0.0;
}

//==========================================================
double SAL::BetheHeitler(double x, double QQ, double msq) {
  double moq2 =4*msq/QQ;
  if(x >= 1/(1+moq2)) return 0;
  double xb=1-x;
  double beta=sqrt(1-moq2*x/xb);
  return
    1.5/M_PI
    *(
      beta*(x*xb*(8-moq2) - 1)
      + log((1+beta)*(1+beta)/(moq2*x)*xb)
        *(x*x + xb*xb + moq2*x*(1-(3 + 0.5*moq2)*x))
    );
}

//==========================================================
double SAL::PSSuppr(double x, double QQ, double msq) {
//--- BH0
  double moq2 =4*msq/QQ;
  if(x >= 1/(1+moq2)) return 0;

  double xb=1-x;
  double beta=sqrt(1-moq2*x/xb);
  double BH0 = beta*(x*xb*8 - 1)
      + log((1+beta)*(1+beta)/(moq2*x)*xb)
        *(x*x + xb*xb);
//--- BHas
  double BHas = log(4*xb/x/moq2)* (x*x+xb*xb) +8*x*xb -1;

 return BH0/BHas;
}

//==========================================================
double SAL::F2(double x, double QQ, int hflav /* =0 */) {
/*
  hflav = 0 --- all
  hflav <= 3 --- all light: d+u+s
  hflav > 3 --- single heavy flavour
*/
  double t, val=0, f[7], df[7];
  double PSfac, Efac;
  int pn,pn0,pn1;

  if(x > 0.995) return 0;
  
  PDF(x,QQ,f);
  t = log(QQ);
  double fac0 = 1/(x*x);
  double fac = 1/(x*(1-x));

  int Aflav = Flav[last_iq];
  for(pn=gluon; pn <= Aflav; pn++) {
    df[pn] = ipol2d(x, t, pn) * (pn ? fac : fac0);
  }
  if(hflav < c_quark) {
    for(pn=d_quark; pn <= s_quark; pn++) val += 2*eq2[pn]*f[pn];
    val += 2*(df[1] + df[0]/3);
    if(hflav) return x*val; //--- light only
  }
  if(hflav == 0) {pn0 = c_quark; pn1 = Nf;} //--- all flavours
  else pn0 = pn1 = hflav; //--- one heavy flavour only
  for(pn=pn0; pn <= pn1; pn++) {
    Efac = evSuppr(QQ, pn);
    val += 2*Efac*eq2[pn]*eq2[pn]*BetheHeitler(x, QQ, m2[pn]);
    if(pn > Aflav) continue;
    PSfac = PSSuppr(x, QQ, m2[pn]);
    val += eq2[pn]*PSfac*( 2*(1-Efac)*(f[pn]+df[pn]) + df[0] );
  }
  return x*val;
}

//========================================================
double SAL::ipol2(double x, double t, int pn) {
  double g0;
  //cout << "fi2[pn][last_iq] " << fi2[pn][last_iq] << endl;
  //cout << "X " << X << endl;
  g0 = SplineValue(X, fi[pn][last_iq], fi2[pn][last_iq], nx, x);
  //cout << "g0 " << g0 << endl;
  return g0 + (SplineValue(X, fi[pn][last_iq+1], fi2[pn][last_iq+1], nx, x)-g0)
    /(T[last_iq+1]-T[last_iq])*(t-T[last_iq]);
}

//========================================================
double SAL::ipol2d(double x, double t, int pn) {
  double g0;
  g0 = SplineValue(X, dfi[pn][last_iq], dfi2[pn][last_iq], nx, x);
  return g0 + (SplineValue(X, dfi[pn][last_iq+1], dfi2[pn][last_iq+1], nx, x)-g0)
    /(T[last_iq+1]-T[last_iq])*(t-T[last_iq]);
}

//================================================
int SAL::InLoc(double X[], int n, double x) {
/* Assumed:
  n > 1, X[i] increasing, x inside [X[0], X[n-1]] or very close to the limits.
  returns ind:
     X[ind] < x <= X[ind+1] for x in (X[0], X[n-1]]
     ind = 0 for x <= X[0]
     ind =n-2 for x >= X[n-1]
*/
	int h,l,m;

  if(x <= X[0]) return 0;
  if(x >= X[n-1]) return n-2;
	l = -1;
	h = n;
	while (h-l > 1) {
		m = (h+l) >> 1;
		if (x > X[m]) l = m;
		else h = m;
	}
	return l;
}

//====================================================
void SAL::TriDiagSolve(double diag[], double above[], double below[], double b[], int n) {
  //--- solve Ax = b
  //--- on exit b is replaced by the solution x
  int j;
  double f;
  //--- forwards
  for(j=1; j < n; j++) {
    f = below[j]/diag[j-1]; //--- A_j,j-1/A_j-1,j-1
    diag[j] -= f*above[j];  //--- f*A_j-1,j
    b[j] -= f*b[j-1];
  }
  b[n-1] /= diag[n-1];
  //--- backwards
  for(j=n-2; j >= 0; j--) b[j] = (b[j] - above[j+1]*b[j+1])/diag[j];
}

//====================================================
double* SAL::InitSpline(double x[], double y[], int n, double yp0) {
  double *spp, *diag, *below, *above;
  try {
    spp = new double[n];
    diag = new double[n];
    below = new double[n];
    above = new double[n];
  }
  catch(...) {
    cerr << "InitSpline: Could not allocate spp etc. Bye ..." << endl;
    exit(1);
  }
  int N = n-1, j;
  double hj, hj1;
  for(j=1; j < N; j++) {
    hj = x[j] - x[j-1];
    hj1 = x[j+1] - x[j];
    below[j] = hj/6;
    diag[j] = (hj+hj1)/3;
    above[j+1] = hj1/6;
    spp[j] = (y[j+1]-y[j])/hj1 + (y[j-1]-y[j])/hj;
  }

  //--- use ds/dx at x[0]
  diag[0] = (x[1] - x[0])/3;
  above[1] = diag[0]/2;
  spp[0] = (y[1]-y[0])/(x[1] - x[0]) - yp0;

  #ifdef _USE_spN
  //--- use ds/dx at x[N]
  hj = x[N] - x[N-1];
  diag[N] = hj/3;
  below[N] = diag[N]/3;
  spp[N] = ypN - (y[N]-y[N-1])/hj;
  #else
  //--- use spp = 0 at x[N]
  diag[N] = 1;
  below[N] = spp[N] = 0;
  #endif

  TriDiagSolve(diag, above, below, spp, n);
  delete[] diag;
  delete[] below;
  delete[] above;
  return spp;
}

//====================================================
double SAL::SplineValue(double x[], double y[], double M[], int npt, double xi) {
  int n = InLoc(x,npt,xi)+1;
  double t,h;
  h = (x[n]-x[n-1]);
  t = (xi-x[n-1])/h;
  return
    y[n-1] + (y[n]-y[n-1])*t
    + h*h/6 * t *(((M[n]-M[n-1])*t + 3*M[n-1])*t -2*M[n-1]-M[n]);
}

//-------  non-class functions  ---------------------

static SAL* sal;

//====================================================
void SALPDF(double x, double QQ, double f[]) {
  if(!sal) sal = new SAL();
  sal->PDF(x, QQ, f);
}

//====================================================
double SALF2(double x, double QQ) {
  if(!sal) sal = new SAL();
  return sal->F2(x, QQ);
}

//====================================================
void ShowSALversion() {
  printf("SAL version %d.%02d\n", MAJOR_SAL_VERSION, MINOR_SAL_VERSION);
}

//-------  for FORTRAN calls  ---------------------

extern "C" {
  void salpdf_(double *x, double *QQ, double *f, int flen) {
    SALPDF(*x, *QQ, f);
  }
  double salf2_(double *x, double *QQ) {
    return SALF2(*x, *QQ);
  }
}
