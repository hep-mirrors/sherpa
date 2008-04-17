//Complex header

#ifndef mycomplex_h
#define mycomplex_h

#include <complex>

#ifdef __GNUC__
typedef std::complex<double> Complex;
#else
typedef std::complex<double> Complex;
#endif

#endif
