#include "HADRON_RESCATTERING/XSecs/NPion.H"
#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include <iostream>


NPion::NPion()
{
//default ctor
}

NPion::~NPion()
{
//default dctor. 
}

double NPion::pPiPlus(double const &s)
{

    const double RES_M0 = 1.23055;     // Nominal resonance mass (GeV)
    const double GAMMA0 = 0.1122;      // Nominal total width at m0 (GeV)
    const int    L = 0;                // angular momentum (s-wave)
    const double PROTON_MASS = 0.938272;
    const double PION_MASS = 0.13957;
    const double GevToMB = 0.389379;   // GeV^-2 -> mb

    const double fullDecayWidth = GAMMA0;     // GeV (nominal)
    const double ResonanceMass  = RES_M0;     // GeV

    if (s <= 0.0) return 0.0;

    double sqrt_s = std::sqrt(s);

    // HPR1R2 parametrization
    if ( sqrt_s > 5.0 ) {
        return m_hpr1r2.xs_tot(HADRON_RESCATTERING::hpr1r2::pPiPlus, s);
    }

    // Compute CM momentum for incoming p + pi at invariant mass m = sqrt(s)
    double m = sqrt_s;
    double pcm = p_CM(m, PROTON_MASS, PION_MASS);

    // If below threshold no resonance
    if (pcm <= 0.0) {
        return 0.0;
    }

    // partial and total widths at the current mass m
    double Gamma_partial = partialWidth_massDep(m, RES_M0, GAMMA0, L, PROTON_MASS, PION_MASS);
    double Gamma_tot     = totalWidth(m,RES_M0, GAMMA0, L, PROTON_MASS,PION_MASS);

    // If widths are zero ret 0
    if (Gamma_partial <= 0.0 || Gamma_tot <= 0.0) {
        return 0.0;
    }

    // spin factor: (2S_R+1) / ((2S_A+1)(2S_B+1)) from the paper
    // For Delta(1232) S_R = 3/2 -> 2S_R+1 = 4; proton S_A=1/2 -> 2S_A+1=2; pion S_B=0 -> 1
    const double spinFactor = 4.0 / (2.0 * 1.0);

    // sigma = pi / p_cm^2 * spinFactor * (Gamma_partial * Gamma_tot) / ((mR - m)^2 + 1/4 Gamma_tot^2)
    double denom = (ResonanceMass - m)*(ResonanceMass - m) + 0.25 * Gamma_tot * Gamma_tot;
    if (denom <= 0.0) return 0.0;

    double sigmaGeV2 = M_PI / (pcm * pcm) * spinFactor * (Gamma_partial * Gamma_tot) / denom;

    // convert to mb
    double sigma_mb = sigmaGeV2 * GevToMB;

    // include non-resonant background contributions here.
    return sigma_mb;

}