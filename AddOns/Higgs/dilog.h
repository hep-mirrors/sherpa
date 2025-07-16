#ifndef DILOG_H
#define DILOG_H

#include "ATOOLS/Math/MyComplex.H"

const double PI = M_PI;
const double E = M_E;

const double SMALL = 1.0e-15;
const double TINY =  1.0e-30;

const double PI_DEF = 3.1415926535897932385;
const double PISQ   =  9.86960440108935861883;
const double ZETA2  =  1.64493406684822643647;
const double PISQ6  =  1.64493406684822643647;
const double PI4    = 97.409091034002437236440;
const double ZETA3  =  1.2020569031595942855;
const double ZETA4  =  1.082323233711138191516;
const double ZETA5  =  1.036927755143369926331;
const double SQRT2  =  1.4142135623730950488;
const double SQRTPI =  1.7724538509055160273;

double li2(double x);
double dilog(double x);
double ReLi2(Complex x);
double ImLi2(Complex x);
Complex CLi2(Complex x);
double li3(double x);
double S12(double x);  // only gets real part correct for x > 1.
double li4(double x);
double myli2(double x);
double mydilog(double x);
Complex i3_3m(double x, double y, double z);
double fastCl(double x);

Complex Clog(double s1, double s2); // Complex log of a ratio
Complex Clog1(double s);   // Complex log of a single s
Complex L0(double s1, double s2);
Complex L0(double s1, double s2);
Complex L1(double s1, double s2);
Complex L2(double s1, double s2);
// Li_2(1-r), for r = s1/s2:
Complex CLi2r(double s1, double s2);
Complex Lsm1(double s1, double s2, double s3, double s4);
Complex Ls0(double s1, double s2, double s3, double s4);
Complex Ls1(double s1, double s2, double s3, double s4);
Complex Ls2(double s1, double s2, double s3, double s4);
Complex Ls3(double s1, double s2, double s3, double s4);

# endif    /*  DILOG_H   */
