#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include <iostream>


double A(double m, double m0, double Gamma) {
    double denom = std::pow(m*m - m0*m0, 2) + 0.25 * Gamma * Gamma;
    if (denom <= 0.0) return 0.0;
    return (1.0 / (2.0 * M_PI)) * (Gamma / denom);
}

double p_CM(double m, double mA, double mB) {
    if (m <= mA + mB) return 0.0;
    double term1 = m*m - (mA + mB)*(mA + mB);
    double term2 = m*m - (mA - mB)*(mA - mB);
    double val = std::sqrt(term1 * term2) / (2.0 * m);
    return val;
}

double partialWidth_massDep(double m, double m0,
                            double Gamma0,
                            int l,
                            double mA, double mB)
{
    double p_m   = p_CM(m,  mA, mB);
    double p_m0  = p_CM(m0, mA, mB);
    if (p_m0 <= 0.0) return 0.0; // should not happen if m0 is above threshold
    if (p_m <= 0.0) return 0.0;  // below threshold -> no decay

    double pow_ratio = std::pow(p_m, 2*l + 1) / std::pow(p_m0, 2*l + 1);

    double damping = 1.2 / (1.0 + 0.2 * (p_m / p_m0));

    double result = Gamma0 * (m0 / m) * pow_ratio * damping;
    return result;
}

double totalWidth(double m,double RES_M0, double GAMMA0, double L, double m1, double m2)
{
    //m1 proton
    //m2 pion
    double g_Npi = partialWidth_massDep(m, RES_M0, GAMMA0, L, m1, m2);
    double Gamma_tot = g_Npi; // add other channels when necessary
    return Gamma_tot;
}

// double A_of_m(double m, double ResonanceMass)
// {
//     double Gamma = totalWidth(m);
//     double denom = (m*m - ResonanceMass*ResonanceMass);
//     denom = denom*denom + 0.25 * Gamma * Gamma;
//     if (denom <= 0.0) return 0.0;
//     return (1.0 / (2.0 * M_PI)) * (Gamma / denom);
// }