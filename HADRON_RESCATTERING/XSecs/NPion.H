#ifndef HADRON_RESCATTERING_NPION_H
#define HADRON_RESCATTERING_NPION_H

#include "HADRON_RESCATTERING/XSecs/Helper.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "HADRON_RESCATTERING/XSecs/HeraParam.H"
#include <cmath>
#include <iostream>
#include <math.h>
#include <complex>
#include <vector>

using namespace HADRON_RESCATTERING;

class NPion 
{
    private:
    const double PROTON_MASS = 0.938272;
    const double PION_MASS   = 0.13957;  
    const double DELTA1232_MASS = 1.234;
    const double N1440_MASS =1.440;
    const double RHO_MASS = 0.77;
    const double ETA_MASS = 0.547;
    const double KAON_MASS = 0.493;
    const double NEUTRON_MASS = 0.9396;
    const double GevToMB     = 0.389379338; // GeV^-2 -> mb
    HADRON_RESCATTERING::HPR1R2  m_hpr1r2;
    HeraParam m_hera;

    struct ResonanceParams {
        double mass;           // Nominal resonance mass (GeV)
        double width0;         // Nominal total width at m0 (GeV)
        std::vector<double> couplingFactor; // The factor multiplied by the Breit-Wigner formula
        double twoSRplusOne;   // 2*S_R + 1
        double L;
        double cgCoeff;
        std::vector<double> m1;
        std::vector<double> m2;
        std::vector<int> LFactors; 
        //for L factors, P_R = P_1 P_2 (-1)^L . 
    };
    //m_pPiPlusResonances

    const std::vector<ResonanceParams> pPiPlus_nPiMinus;
    std::vector<ResonanceParams> pPiMZ_nPiPZ;

    void initForProtonPion(); 
    double calculateResonanceSigma(double sqrt_s, const ResonanceParams& res, double massNeutron, double massPion) const;
    
    double calculateResonanceSigmaForProtonPion(double sqrt_s, const ResonanceParams& res, double massProton, double massPion, int channelIndex) const;

    double calculateResonanceSigmaForPionNought(double sqrt_s, const ResonanceParams& res, double m1, double m2) const;
    
    double calculateResonanceSigmaForNeutronPion(double sqrt_s, const ResonanceParams& res, double m1, double m2,int channelIndex) const;

    public:

    NPion();
    ~NPion();

    //total scattering cross sections
    //pPi
    double pPiMinus(const double & s);
    double pPiPlus(const double & s);
    double pPiZero(const double &s);
    double pPiZeroWignerEckart(const double &s);
    
    //npi
    double nPiMinus(const double & s);
    double nPiPlus(const double & s);
    double nPiZero(const double & s);


};

#endif