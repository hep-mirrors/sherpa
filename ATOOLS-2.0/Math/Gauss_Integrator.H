#ifndef Gauss_Integrator_h
#define Gauss_Integrator_h

#include "Function_Base.H"    

namespace AMATOOLS {

  struct weightlist {
    int methode;
    int n;
    double * w;
    double * x;
    weightlist * next;
  };
  
  class Gauss_Integrator {
  protected:
    // weights or abscissas already calculated? available 
    static int n_gauleg,n_gaulag,n_gauherm,n_gaujac;  
    static weightlist* wlistroot;
    
    // status of one instance of the calculator
    int mode; int numberabsc;        // type of weights , number of abscissas
    Function_Base * func;
    weightlist* wlistact;
    void gauleg(double x[], double w[], int n);    // mode=1;
    void gaujac(double x[], double w[], int n, double alf = -0.5f, double bet= -0.5f);  // mode=5;
    double gammln(double xx);           // should be turned to static
  public:
    Gauss_Integrator(Function_Base *);
    double Legendre(double x1 , double x2,int n);
    double Jacobi(double x1, double x2, int n, double alf, double bet);
    double Chebyshev( double a, double b, double prec, int N_Max, int &I_Err );
    double Integrate(double x1, double x2, double prec, int mode=1, int nmax=65536 );
    //    double Integrate(double x1, double x2, int mode, int n);
  };
}

#endif

