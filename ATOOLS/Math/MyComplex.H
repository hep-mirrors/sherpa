#ifndef ATOOLS_Math_MyComplex_H
#define ATOOLS_Math_MyComplex_H
//Complex header
#include <complex>
typedef std::complex<double> Complex;
inline std::ostream& operator<<(std::ostream& ss, const Complex& c) {
    double a(c.real());
    double b(c.imag());
    if (a == 0.) ss << "i*" << b;
    else if (b == 0.) ss << a;
    else ss << a << " + i*" << b;
    return ss;
}
inline Complex complexconjugate(const Complex& c) {return Complex(c.real(), -c.imag());}
inline double ToDouble(const Complex& c) {return c.real();}

#endif
