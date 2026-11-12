#ifndef HADRON_RESCATTERING_NKaon_H
#define HADRON_RESCATTERING_NKaon_H

#include <math.h>
#include <vector>
#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "HADRON_RESCATTERING/XSecs/HeraParam.H"


//resonance implementation follows from 9803035

using namespace HADRON_RESCATTERING;

class NKaon 
{
    private:

    const double PROTON_MASS     = 0.938272;
    const double NEUTRON_MASS    = 0.939565;
    const double XI_MASS = 1.317 ;
    const double KAON_MASS       = 0.493677; // Charged Kaon (K-)
    const double PION_MASS       = 0.139570; // Charged Pion (pi+/-)
    const double LAMBDA_MASS     = 1.11568;  // Ground state Lambda
    const double SIGMA_MASS      = 1.18937;  // Ground state Sigma
    const double SIGMA_STAR_MASS = 1.38500;  // Sigma(1385) resonance
    const double ETA_MASS        = 0.54786;  // Eta meson
    const double KSTAR_MASS      = 0.89166;  // K*(892) meson
    const double DELTA1232_MASS = 1.232;
    const double GevToMB     = 0.389379338; // GeV^-2 -> mb

    struct ResonanceParams
    {
        double mass;           // Nominal resonance mass (GeV)
        double width0;         // Nominal total width at m0 (GeV)
        std::vector<double> branchingRatio; // The factor multiplied by the Breit-Wigner formula
        double twoSRplusOne;   // 2*S_R + 1
        double L;
        double probFactor;
        std::vector<double> m1;
        std::vector<double> m2;
        std::vector<int> LFactors; 
    };

    std::vector<ResonanceParams> m_Resonances;
    std::vector<ResonanceParams> m_Resonances_n_KaonMinus;

    double calculateResonanceScatteringXSec(double sqrt_s, const ResonanceParams& res, double massProton, double massPion) const;

    public:
    NKaon();
    HPR1R2  m_hpr1r2;

    double pKaonMinus(const double &s);
    double nKaonMinus(const double &s);


};








#endif