#ifndef ATOOLS_Math_MyComplex_H
#define ATOOLS_Math_MyComplex_H
//Complex header
#include <complex>
typedef std::complex<double> Complex;
inline Complex complexconjugate(const Complex& c) {return Complex(c.real(), -c.imag());}
inline double ToDouble(const Complex& c) {return c.real();}

namespace ATOOLS
{
    //proper output redefined
    inline std::ostream& operator<<(std::ostream& ss, const Complex& c) {
        double a(c.real());
        double b(c.imag());
        if (a != 0.) {
            if (b != 0) ss << a << "+i*" << b;
            else ss << a;
        }
        else {
            if (b != 0) ss << "i*" << b;
            else ss << "0";
        }
        return ss;    
    }
} // namespace ATOOLS

#endif
