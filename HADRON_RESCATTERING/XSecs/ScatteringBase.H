#ifndef HADRON_RESCATTERING_XSecs_ScatteringBase_H
#define HADRON_RESCATTERING_XSecs_ScatteringBase_H

#include "ATOOLS/Phys/Flavour.H"
#include <vector>

class ScatteringBase 
{
    public:
    ScatteringBase();
    virtual ~ScatteringBase() = default;
    virtual double XStot(const ATOOLS::Flavour & A, const ATOOLS::Flavour & B, const double & s) = 0;
    static bool initialized;  
    void Initialize();

    protected:
    const double PROTON_MASS     = 0.938272;
    const double NEUTRON_MASS    = 0.939565;
    const double XI_MASS         = 1.317;
    const double KAON_MASS       = 0.493677; // Charged Kaon (K-)
    const double PION_MASS       = 0.139570; // Charged Pion (pi+/-)
    const double LAMBDA_MASS     = 1.11568;  // Ground state Lambda
    const double SIGMA_MASS      = 1.18937;  // Ground state Sigma
    const double SIGMA_STAR_MASS = 1.38500;  // Sigma(1385) resonance
    const double ETA_MASS        = 0.54786;  // Eta meson
    const double KSTAR_MASS      = 0.89166;  // K*(892) meson
    const double DELTA1232_MASS  = 1.232;
    const double CHARGEDKAON_MASS   = 0.49368;
    const double KAON0_MASS         = 0.47761;
    const double OMEGA_MASS         = 0.782;
    const double RHO_MASS           = 0.769;
    const double KAONSTAR_MASS      = 0.893;
    const double N1440_MASS         = 1.440;
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

    static std::vector<ResonanceParams> m_Resonances;
    static std::vector<ResonanceParams> m_Resonances_n_KaonMinus;

    static std::vector<ResonanceParams> m_PiPlusPiMinusRes;
    static std::vector<ResonanceParams> m_PiPlusPiZeroRes;
    static std::vector<ResonanceParams> m_PiZeroPiZeroRes;
    static std::vector<ResonanceParams> m_KPlusPiMinusRes;

    static std::vector<ResonanceParams> pPiPlus_nPiMinus;
    static std::vector<ResonanceParams> pPiMZ_nPiPZ;


};





#endif