#ifndef ATOOLS_Math_Gauss_Integrator_H
#define ATOOLS_Math_Gauss_Integrator_H

#include "ATOOLS/Math/Function_Base.H"    

namespace ATOOLS {

  struct Weight_Module {
    int      methode;
    int      n;
    double * w;
    double * x;
    Weight_Module * next;
  };
  
  class Gauss_Integrator {
  protected:
    // weights or abscissas already calculated ? available 
    static int s_ngauleg,s_ngaulag,s_ngauherm,s_ngaujac;  
    static Weight_Module * s_wlistroot;
    
    // status of one instance of the calculator
    int    m_numberabsc;        // number of abscissas
    Function_Base * m_func;
    Weight_Module * m_wlistact;

    void   GauLeg(double * x, double * w, int n);    // mode=1;
    void   GauJac(double * x, double * w, int n, double alf = -0.5, double bet= -0.5);  // mode=5;
  public:
    Gauss_Integrator(Function_Base *a =0);
    double Legendre(double x1 , double x2,int n);
    double Jacobi(double x1, double x2, int n, double alf, double bet);
    double Chebyshev( double a, double b, double prec, int n_max, int &i_err );
    double Integrate(double x1, double x2, double prec, int mode=1, int nmax=65536 );
  };


  /*!
    \file
    \brief contains the class Gauss_Integrator
  */

  /*!
    \class Gauss_Integrator
    \brief This class can be used to evaluate one dimensional integrals numerically.
    
    The Gauss_Integrator can calculate the integral of any given one dimensional
    analytical function (Function_Base) numerically.
    
    Several methods are avaliable:
     - Gauss - Legendre (default)
     - Gauss - Jacobi
     - Gauss - Chebyshev
    .
  */

  /*!
    \fn Gauss_Integrator::Gauss_Integrator(Function_Base *)
    \brief Constructor taking the function to be integrated.
  */

  /*!
    \fn  double Gauss_Integrator::Legendre(double x1 , double x2,int n)
    \brief Gauss - Legendre integration

    Calculates the integral 
    \f[
       \int_{x_1}^{x_2} dx f(x) = \sum^n w_i f(\xi_i)
    \f]
    The weights \f$w_i\f$ and abscissas \f$\xi_i\f$ are such, that in case 
    the function \f$f(x)\f$ is a polynomial of order k < n 
    the result is exact!
  */

  /*!
    \fn double Gauss_Integrator::Jacobi(double x1, double x2, int n, double alf, double bet)
    \brief Gauss-Jacobi integration
    
    This integrator is optimized to integrate function of the Form
    \f[
      f(x) = (1-x)^\alpha (1+x)^\beta  P(x)
    \f]
    where \f$P(x)\f$ is a polynominal.
  */


  /*!
    \fn double  Gauss_Integrator::Chebyshev( double a, double b, double prec, int N_Max, int &I_Err )
    \brief Gauss-Chebyshev integration

    This integrator is optimized to integrate function of the Form
    \f[
      f(x) = (1-x^2)^{-1/2} P(x)
    \f]
    where \f$P(x)\f$ is a polynominal.
  */

  /*!
    \fn   double  Gauss_Integrator::Integrate(double x1, double x2, double prec, int mode=1, int nmax=65536 )
    \brief Integrates a function between x1 and x2 nummerically up to precission prec.

    Calculates the integral 
    \f[
       \int_{x_1}^{x_2} dx f(x)
    \f]
    utilizing either a Legendre, a Chebyshev or a Jacobi method.
  */


}

#endif

