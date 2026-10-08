#ifndef ATOOLS_Math_MyComplex_H
#define ATOOLS_Math_MyComplex_H
//Complex header
#include <complex>
typedef std::complex<double> Complex;
inline Complex complexconjugate(const Complex& c) {return Complex(c.real(), -c.imag());}
inline double ToDouble(const Complex& c) {return c.real();}

#endif
