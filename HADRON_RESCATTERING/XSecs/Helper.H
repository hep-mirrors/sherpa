#ifndef HADRON_RESCATTERING_HELPER_H
#define HADRON_RESCATTERING_HELPER_H

#include <math.h>


// Breit-Wigner mass distribution A(m), not used for narrow dist.
double A(double m, double m0, double Gamma) ;
//center of mass p.
double p_CM(double m, double mA, double mB) ;
//mass dependent partial width, equation 24
double partialWidth_massDep(double m, double m0,double Gamma0,int l,double mA, double mB);
//total width
double totalWidth(double m,double RES_M0, double GAMMA0, double L, double m1, double m2);
// Breit-Wigner mass distribution A(m) per Eq. (27) (kept for completeness)
double A_of_m(double m, double ResonanceMass);






#endif