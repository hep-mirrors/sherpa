#ifndef HADRON_RESCATTERING_HELPER_H
#define HADRON_RESCATTERING_HELPER_H


#include "ATOOLS/Math/Vec4.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Vector.H"
#include <iostream>
#include <math.h>
#include <vector>


    ////////https://arxiv.org/pdf/2005.05658
    ////////Equation 24,25,26,27 are here. They are used to calculate equation 23.

// Breit-Wigner mass distribution A(m), not used for narrow dist.
double A(double m, double m0, double Gamma) ;
//center of mass p.
double p_CM(double m, double mA, double mB) ;
//mass dependent partial width, equation 24
double partialWidth_massDep(double m, double m0,double Gamma0,double l,double mA, double mB);
//total width


double totalWidth(double m,double RES_M0, double GAMMA0, std::vector<double> br, double L, double m1, double m2);

double totalWidth(double m,double RES_M0, double GAMMA0, double L, double m1, double m2);

double totalWidth(double m, double RES_M0, double GAMMA0, 
                  const std::vector<double> br, 
                  double L, 
                  const std::vector<double> m1, 
                  const std::vector<double> m2, 
                  const std::vector<int> LVector);
// Breit-Wigner mass distribution A(m) per Eq. (27) (kept for completeness)
double A_of_m(double m, double ResonanceMass);

#endif
