#ifndef ATOOLS_Math_MathTools_H
#define ATOOLS_Math_MathTools_H
/*  Declarations for discrete functions  */

#ifdef __GNUC__
// GNU C++ Compiler       
#include <cmath>        
#include <cstdlib>

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
#include <iostream>
#include <math.h> 
#endif


#include "ATOOLS/Math/MyComplex.H"

namespace ATOOLS {

  template <class Type> const Type &Min(const Type &a,const Type &b) 
  { return a<b?a:b; }
  template <class Type> const Type &Max(const Type &a,const Type &b) 
  { return a>b?a:b; }

  template <class Type> Type &Min(Type &a,Type &b) 
  { return a<b?a:b; }
  template <class Type> Type &Max(Type &a,Type &b) 
  { return a>b?a:b; }

  inline double Accu() {return 1.e-12;}

  inline int    Sign(const int& a)    { return a<0?-1:1;       }
  inline double Sign(const double& a) { return a<0.0?-1.0:1.0; }

  inline double Theta(const double &a) { return a<0.0?0.0:1.0; }

  inline int    iabs(const int& a)    { return a<0?-a:a;   } 
  inline double dabs(const double& a) { return a<0.0?-a:a; } 

  template <typename Scalar>
  inline Scalar sqr(const Scalar &x)   { return x*x; } 
  template <typename Scalar> inline std::complex<Scalar> 
  csqr(const std::complex<Scalar> &x) { return x*x; } 


  inline int IsZero(const double &a,const double &crit) 
  { return dabs(a)<crit?1:0; }
  inline int IsZero(const Complex &a,const double &crit) 
  { return std::abs(a)<crit?1:0; }

  inline int IsEqual(const double &a,const double &b) 
  {
    if (a==0. && b==0.) return 1;
    return (dabs(a-b)/(dabs(a)+dabs(b))<Accu()) ? 1 : 0;
  }
  inline int IsEqual(const double &a,const double &b,const double &crit) 
  {
    if (a==0. && b==0.) return 1;
    return (dabs(a-b)/(dabs(a)+dabs(b))<crit) ? 1 : 0;
  }
  inline int IsEqual(const Complex &a,const Complex &b) 
  {
    if (a==Complex(0.,0.) && b==Complex(0.,0.)) return 1;
    return (std::abs(a-b)/(std::abs(a)+std::abs(b))<Accu()) ? 1 : 0;
  } 
  inline Complex csqrt(const double &d)
  {
    if (d<0) return Complex(0.,sqrt(-d));
    return sqrt(d);
  }

#define GAMMA_E 0.5772156649015328606

  double Gammln(double xx);

  double ReIncompleteGamma0(double x,double prec=1.e-6);

  double  DiLog(double x);
  Complex DiLog(const Complex& x);

  int Factorial(const int n);

  double ExpIntegral(int n, double x);

  template<typename Scalar> inline bool IsNan(const Scalar& x);
  template<typename Scalar> inline bool IsBad(const Scalar& x);
  template<typename Scalar> inline bool IsZero(const Scalar& x);
  template<typename Scalar> inline Scalar Abs(const Scalar& x);
  template<typename Scalar> inline Scalar Abs(const std::complex<Scalar>& x);

  template<> inline bool IsNan<double>(const double& x) { 
    return std::isnan(x)||std::isnan(-x);
  }
  template<> inline bool IsBad<double>(const double& x) { 
    return IsNan(x)||std::isinf(x)||std::isinf(-x);
  }
  template<> inline bool IsZero<double>(const double& x) {
    return dabs(x)<Accu()?1:0;
  }
  template<> inline double Abs<double>(const double& x) {
    return x>0.0?x:-x;
  }

  template<> inline bool IsNan<long double>(const long double& x) { 
    return std::isnan(x)||std::isnan(-x);
  }
  template<> inline bool IsBad<long double>(const long double& x) { 
    return IsNan(x)||std::isinf(x)||std::isinf(-x);
  }
  template<> inline bool IsZero<long double>(const long double& x) {
    return dabs(x)<Accu()?1:0;
  }
  template<> inline long double Abs<long double>(const long double& x) {
    return x>0.0?x:-x;
  }


  template<> inline bool IsNan<Complex>(const Complex& x) {
    return (std::isnan(real(x)) || std::isnan(imag(x)) ||
	    std::isnan(-real(x)) || std::isnan(-imag(x))); 
  }
  template<> inline bool IsBad<Complex>(const Complex& x) { 
    return IsNan(x)||std::isinf(real(x))||std::isinf(imag(x))
      ||std::isinf(-real(x))||std::isinf(-imag(x));
  }
  template<> inline bool IsZero<Complex>(const Complex& x) {
    return std::abs(x)<Accu()?1:0;
  }
  template<> inline double Abs<double>(const Complex& x) {
    return std::abs(x);
  }

  template<> inline bool IsNan<std::complex<long double> >
  (const std::complex<long double>& x) {
    return (std::isnan(real(x)) || std::isnan(imag(x)) ||
	    std::isnan(-real(x)) || std::isnan(-imag(x))); 
  }
  template<> inline bool IsBad<std::complex<long double> >
  (const std::complex<long double>& x) { 
    return IsNan(x)||std::isinf(real(x))||std::isinf(imag(x))
      ||std::isinf(-real(x))||std::isinf(-imag(x));
  }
  template<> inline bool IsZero<std::complex<long double> >
  (const std::complex<long double>& x) {
    return std::abs(x)<Accu()?1:0;
  }
  template<> inline long double Abs<long double>
  (const std::complex<long double>& x) {
    return std::abs(x);
  }



  template<class T1, class T2>
  struct promote_trait {
  };

#define DECLARE_PROMOTE(A,B,C)           \
  template<> struct promote_trait<A,B> { \
    typedef C T_promote;                 \
  }

  DECLARE_PROMOTE(double,Complex,Complex);
  DECLARE_PROMOTE(Complex,double,Complex);
  DECLARE_PROMOTE(int,double,double);
  DECLARE_PROMOTE(double,int,double);
  DECLARE_PROMOTE(int,Complex,Complex);
  DECLARE_PROMOTE(Complex,int,Complex);
  
  DECLARE_PROMOTE(double,double,double);
  DECLARE_PROMOTE(Complex,Complex,Complex);

  DECLARE_PROMOTE(long double,std::complex<long double>,
		  std::complex<long double>);
  DECLARE_PROMOTE(std::complex<long double>,long double,
		  std::complex<long double>);
  DECLARE_PROMOTE(int,long double,long double);
  DECLARE_PROMOTE(long double,int,long double);
  DECLARE_PROMOTE(int,std::complex<long double>,std::complex<long double>);
  DECLARE_PROMOTE(std::complex<long double>,int,std::complex<long double>);
  
  DECLARE_PROMOTE(long double,long double,long double);
  DECLARE_PROMOTE(std::complex<long double>,std::complex<long double>,
		  std::complex<long double>);

#define PROMOTE(TYPE1,TYPE2) typename promote_trait<TYPE1,TYPE2>::T_promote

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

  /*!
    \fn    double ReIncomplietGamma0(double xx)
    \brief calculates the real part of the incomplete Gammafunction.
  */

  /*!
    \fn    double DiLog(double x)
    \brief calculates the real part of Li_2(x).
  */

}

#endif
