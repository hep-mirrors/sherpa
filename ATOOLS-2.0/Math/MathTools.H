/*  Declarations for diskrete functions  */
#ifndef mathtools1_h
#define mathtools1_h

#ifdef __GNUC__
// GNU C++ Compiler       
#include <cmath>        

/*
  if __GNUC__ == 3 && __GNUC_MINOR__ == 0.
  if defined __GNUC__ && defined __cplusplus && __GNUC_MINOR__ >= 8
  if !defined __GNUC__ || __GNUC__ < 2 || __GNUC_MINOR__ < 7 
  #define GCC_VERSION (__GNUC__ * 10000 \
  + __GNUC_MINOR__ * 100 \
  + __GNUC_PATCHLEVEL__)
  ...
  Test for GCC > 3.2.0 
  #if GCC_VERSION > 30200 
*/

#endif


#if defined(__sgi) && !defined(__GNUC__)
// SGI IRIX C++ Compiler, complex but not double methods need "std::", e.g. Abs() exp()
#include <math.h> 
#endif


#include "MyComplex.H"

namespace ATOOLS {

  template <class Type>
  inline Type Min(Type a, Type b) { return a<b ? a : b; }
  template <class Type>
  inline Type Max(Type a, Type b) { return a>b ? a : b; }

  /*
  inline int Min(int a, int b) { return a<b ? a : b; }
  inline int Max(int a, int b) { return a>b ? a : b; }
  inline int Min(long int a, long int b) { return a<b ? a : b; }
  inline int Max(long int a, long int b) { return a>b ? a : b; }

  inline double Min(double a, double b) { return a<b ? a : b; }
  inline double Max(double a, double b) { return a>b ? a : b; }
  inline long double Min(long double a, long double b) { return a<b ? a : b; }
  inline long double Max(long double a, long double b) { return a>b ? a : b; }
  */

  inline int         Sign(const int& a) {return (a<0) ? -1 : 1;}
  inline int         iabs(const int& a) {return a>0 ? a : -a;} 
  inline double      dabs(const double& a) {return a>0 ? a : -a;} 
  inline double      sqr(double x) {return x*x;} 
  inline double      Accu() {return 1.e-12;};
  inline int IsZero(const double a) {
    return (dabs(a)<Accu()) ? 1 : 0;
  }
  inline int IsZero(const Complex& a) {
    return (std::abs(a)<Accu()) ? 1 : 0;
  }
  inline int IsEqual(const double a,const double b) {
    if (a==0. && b==0.) return 1;
    return (dabs(a-b)/(dabs(a)+dabs(b))<Accu()) ? 1 : 0;
  }
  inline int IsEqual(const Complex& a,const Complex& b) {
    if (a==Complex(0.,0.) && b==Complex(0.,0.)) return 1;
    return (std::abs(a-b)/(std::abs(a)+std::abs(b))<Accu()) ? 1 : 0;
  } 
  inline Complex csqrt(const double d)
  {
    if (d<0) return Complex(0.,sqrt(-d));
    return sqrt(d);
  }
  inline Complex csqr(Complex x) {return x*x;} 

#define GAMMA_E 0.5772156649015328606

  double Gammln(double xx);


  /*!
    \file
    \brief contains a collection of simple mathematical functions
  */

  /*!
    \fn inline Type Min(Type a, Type b)
    \brief returns the minimum of two numbers  
  */

  /*!
    \fn inline Type Max(Type a, Type b) 
    \brief returns the maximum of two numbers
  */

  /*! 
    \fn  inline int         Sign(const int& a) {return (a<0) ? -1 : 1;}
    \brief returns the sign of the argument
  */

  /*! 
    \fn  inline int         iabs(const int& a) {return a>0 ? a : -a;} 
    \brief returns the absolute value of the argument
  */

  /*! 
    \fn  inline double      dabs(const double& a) {return a>0 ? a : -a;} 
    \brief returns the absolute value of the argument
  */

  /*! 
    \fn  inline double      sqr(double x) {return x*x;} 
    \brief returns the argument squared 
  */

  /*! 
    \fn  inline double      Accu() {return 1.e-12;};
    \brief returns a (platform dependent) precission, default is \f$1^{-12}\f$
  */

  /*! 
    \fn  inline int IsZero(const double a) 
    \brief returns \em true if argument is smaller than Accu()
  */

  /*! 
    \fn  inline int IsZero(const Complex& a)
    \brief  returns \em true if argument is smaller than Accu()
  */

  /*! 
    \fn  inline int IsEqual(const double a,const double b)
    \brief  returns \em true if arguments are equal (compared to Accu())
  */

  /*! 
    \fn  inline int IsEqual(const Complex& a,const Complex& b)
    \brief  returns \em true if arguments are equal (compared to Accu())
  */

  /*! 
    \fn  inline Complex csqrt(const double d)
    \brief returns the complex root of a (possibly negative) float or double variable
  */

  /*! 
    \fn  inline Complex csqr(Complex x)
    \brief returns the argument squared
  */

  /*!
    \fn    double Gammln(double xx)
    \brief calculates the logarithm of the Gammafunction
  */

}

#endif
