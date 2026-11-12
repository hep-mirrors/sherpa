#include "HADRON_RESCATTERING/XSecs/NPion.H"

using namespace ATOOLS;

NPion::NPion() : 
    pPiPlus_nPiMinus({
            {1.232, 0.115, {1.0}, 4.0, 1, (2./3.), 
            {PROTON_MASS}, 
            {PION_MASS},
            {1}}, // RES_1232 (Delta(1232) P33)
            {1.700, 0.200, {0.15, 0.55, 0.30}, 4.0, 1, (2./3.), 
            {PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,PION_MASS,PION_MASS},
            {1, 1, 1}}, // RES_1600 (Delta(1600) P33)
            {1.675, 0.180, {0.25, 0.60, 0.15}, 2.0, 0, (2./3.), 
            {PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,PION_MASS,PION_MASS},
            {0, 2, 0}}, // RES_1620 (Delta(1620) S11)
            {1.750, 0.300, {0.20, 0.10, 0.55, 0.15}, 4.0, 2, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {2, 0, 0, 2}}, // RES_1700 (Delta(1700) D33)
            {1.850, 0.240, {0.30,0.15, 0.30, 0.25}, 2.0, 1, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {0, 0, 2, 0}}, // RES_1900 (Delta(1900) S31)
            {1.880, 0.280, {0.20, 0.60, 0.10, 0.10}, 6.0, 3, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},   
            {3, 1, 1, 3}}, // RES_1905 (Delta(1905) F15)
            {1.900, 0.250, {0.35, 0.40, 0.15, 0.10}, 2.0, 1, (2./3.),
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {1,1,1,1}}, // RES_1910 (Delta(1910) P11)
            {1.920, 0.150, {0.15, 0.30, 0.30, 0.25}, 4.0, 1, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {1, 1, 1, 1}}, // RES_1920 (Delta(1920) P33)
            {1.930, 0.250, {0.20,0.25, 0.25, 0.30}, 6.0, 2, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {2, 2, 2, 2}}, // RES_1930 (Delta(1930) D15)
            {1.950, 0.250, {0.45, 0.15, 0.20,0.20}, 8.0, 3, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {3, 3, 3, 3}} // RES_1950 (Delta(1950) F37)
            //for L calculations: 
            /*
            parity : P = Pn P_pi(-1)^L = -(-1)^L --> if P = + then L is odd, if P = - then L is even.
            1232 --> J^P = 3/2^+ --> P = + --> smallest odd = 1
            1620 --> J^P = 1/2^- --> P = - --> smallest even = 0
            1720 --> J^P = 3/2^- --> P = - --> L even --> smallest even gives J = 1/2 --> then L = 2. 
            1900 --> j^p = 1/2^-
            1905 --> J^P = 5/2^+ --> L =3 
            1910 --> J^P = 1/2^+
            1920 --> J^P = 3/2^+
            1930 --> J^P = 5/2^- --> L=2 contains 5/2
            1950 --> J^p = 7/2^+ --> L=3
            */
    }), pPiMZ_nPiPZ 
    ({
        // for pion^0, I_m = 1, m_m = 0 --> +1 for pion+, 0 for pion^0, -1 for pion^-
        // for proton I_b = 1/2, m_b = 1/2 --> -1/2 is for neutron
        // for delta family we have I_delta = 3/2, M_delta = 1/2 // for N family we have I_R = 1/2, M_R = 1/2 
        // { Mass(GeV),Width0(GeV), BR, 2S_R+1, L , CGCoeff, probFactor, exitChannelParticleMasses,Exit L's }
            {1.232, 0.131, {1.0}, 4.0, 1, (2./3.), 
            {PROTON_MASS}, 
            {PION_MASS},
            {1}}, // RES_1232 (Delta(1232) P33)
            {1.570, 0.250, {0.15, 0.55, 0.30}, 4.0, 1, (2./3.), 
            {PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,PION_MASS,PION_MASS},
            {1, 1, 1}}, // RES_1600 (Delta(1600) P33)
            {1.610, 0.130, {0.25, 0.60, 0.15}, 2.0, 0, (2./3.), 
            {PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,PION_MASS,PION_MASS},
            {0, 2, 0}}, // RES_1620 (N(1620) S11)
            {1.710, 0.300, {0.20, 0.10, 0.55, 0.15}, 4.0, 2, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {2, 0, 0, 2}}, // RES_1700 (Delta(1700) D33)
            {1.860, 0.250, {0.30,0.15, 0.30, 0.25}, 2.0, 1, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {0, 0, 2, 0}}, // RES_1900 (Delta(1900) S31)
            {1.880, 0.330, {0.20, 0.60, 0.10, 0.10}, 6.0, 3, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {3, 1, 1, 3}}, // RES_1905 (N(1905) F15)
            {1.900, 0.300, {0.35, 0.40, 0.15, 0.10}, 2.0, 1, (2./3.),
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {1,1,1,1}}, // RES_1910 (N(1910) P11)
            {1.920, 0.300, {0.15, 0.30, 0.30, 0.25}, 4.0, 1, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {1, 1, 1, 1}}, // RES_1920 (Delta(1920) P33)
            {1.930, 0.300, {0.20,0.25, 0.25, 0.30}, 6.0, 2, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {2, 2, 2, 2}}, // RES_1930 (N(1930) D15)
            {1.930, 0.285, {0.45, 0.15, 0.20,0.20}, 8.0, 3, (2./3.), 
            {PROTON_MASS,PROTON_MASS,DELTA1232_MASS,N1440_MASS}, 
            {PION_MASS,RHO_MASS,PION_MASS,PION_MASS},
            {3, 3, 3, 3}}, // RES_1950 (Delta(1950) F37)
            {1.440, 0.350, {0.70, 0.05, 0.25}, 2, 1, (1./3.), 
            {PROTON_MASS, DELTA1232_MASS, PROTON_MASS},
            {PION_MASS, PION_MASS, RHO_MASS},
            {1, 1, 1}}, // N(1440) 
            {1.515, 0.110, {0.6,0.15, 0.25}, 4, 2, (1./3.), 
            {PROTON_MASS, PROTON_MASS, DELTA1232_MASS}, // mA
            {PION_MASS, RHO_MASS, PION_MASS},
            {2, 0, 0}}, // N(1520) 
            {1.530, 0.150, {0.55,0.35, 0.05,0.05}, 2, 0, (1./3.), 
            {PROTON_MASS, PROTON_MASS, RHO_MASS, N1440_MASS}, // mA
            {PION_MASS, ETA_MASS, PION_MASS, PION_MASS},
            {0, 0, 0, 0}}, // N(1535) 
            {1.650, 0.125, {0.65, 0.05, 0.05, 0.10, 0.05, 0.10, }, 2, 0, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS, 1.116}, // mA (1.116 for Lambda)
            {PION_MASS, ETA_MASS, RHO_MASS, PION_MASS, PION_MASS, KAON_MASS},
            {0, 0, 0, 0, 0, 0}}, // N(1650) 
            {1.675, 0.140, {0.45, 0.55}, 6, 2, (1./3.), 
            {PROTON_MASS, DELTA1232_MASS}, // mA
            {PION_MASS, PION_MASS},
            {2, 2}}, // N(1675) 
            {1.675, 0.145, {0.65, 0.20,0.15}, 6, 3, (1./3.), 
            {PROTON_MASS, PROTON_MASS, DELTA1232_MASS}, // mA
            {PION_MASS, RHO_MASS, PION_MASS},
            {3, 1, 1}}, // N(1680) 
            {1.720, 0.200, {0.10,0.05, 0.05, 0.45, 0.35}, 4, 2, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS}, // mA
            {PION_MASS, ETA_MASS, RHO_MASS, PION_MASS, PION_MASS},
            {2, 2, 0, 0, 2}}, // N(1700) 
            {1.710, 0.140, {0.15, 0.20, 0.05, 0.20, 0.20}, 2, 1, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, 1.116}, // mA
            {PION_MASS, ETA_MASS, RHO_MASS, PION_MASS, KAON_MASS},
            {1, 1, 1, 1, 1}}, // N(1710) 
            {1.720, 0.250, {0.15, 0.25, 0.45, 0.10, 0.05}, 4, 1, (1./3.), 
            {PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS, 1.116}, // mA
            {PION_MASS, RHO_MASS, PION_MASS, PION_MASS, KAON_MASS},
            {1, 1, 1, 1, 1}}, // N(1720) 
            {1.920, 0.200, {0.35, 0.55, 0.05, 0.05}, 4, 1, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS}, // mA
            {PION_MASS, 0.783, RHO_MASS, PION_MASS}, 
            {1, 1, 1, 1}}, // N(1900) 
            {2.020, 0.300, {0.05, 0.15, 0.25, 0.30, 0.15, 0.10}, 8, 3, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS, 1.116}, // mA
            {PION_MASS, ETA_MASS, RHO_MASS, PION_MASS, PION_MASS, KAON_MASS},
            {3, 3, 3, 3, 3, 3}}, // N(1990) 
            // {1.875, 0.200, {0.60, 0.05, 0.25, 0.05, 0.05}, 4, 2, (1./3.), 
            // {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS}, // mA
            // {PION_MASS, ETA_MASS, RHO_MASS, PION_MASS, PION_MASS},
            // {2, 2, 1, 1, 1}}, // N(2080) 
            {2.180, 0.400, {0.35, 0.30, 0.15, 0.15, 0.05}, 8, 4, (1./3.), 
            {PROTON_MASS, PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS}, // mA
            {PION_MASS, RHO_MASS, PION_MASS, PION_MASS, KAON_MASS},
            {4, 2, 2, 2, 2}}, // N(2190) 
            {2.250, 0.400, {0.35,0.25, 0.20, 0.20}, 10, 5, (1./3.), 
            {PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS}, // mA
            {PION_MASS, RHO_MASS, PION_MASS, PION_MASS},
            {5, 4, 4, 4}}, // N(2220) 
            {2.280, 0.500, {0.30, 0.25, 0.20, 0.20, 0.05}, 10, 4, (1./3.), 
            {PROTON_MASS, PROTON_MASS, DELTA1232_MASS, N1440_MASS, 1.116}, // mA
            {PION_MASS, RHO_MASS, PION_MASS, PION_MASS, KAON_MASS},
            {4, 3, 3, 3, 3}} // N(2250) // throw it out.
    })
    {
        //default ctor, alter if necessary. 
    }

NPion::~NPion()
{
    //default dtor
}

double NPion::calculateResonanceSigma(double sqrt_s, const ResonanceParams& res, double m1, double m2) const {
    double m = sqrt_s;
    double pcm = p_CM(m, m1, m2);

    if (pcm <= 0.0) 
    {
        return 0.0;
    }

    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.couplingFactor, res.L, res.m1, res.m2, res.LFactors);

    double spinFactor = res.twoSRplusOne / (2.0 * 1.0);

    double massDiff = res.mass - m;
    double denom = massDiff * massDiff + 0.25 * Gamma_tot * Gamma_tot;
    if (denom <= 0.0) return 0.0;

    // Breit-Wigner Cross-section formula (in GeV^-2)
    double sigmaGeV = res.couplingFactor[0] * M_PI / (pcm * pcm) * spinFactor * (Gamma_tot * Gamma_tot) / denom;

    return sigmaGeV;
}

double NPion::calculateResonanceSigmaForProtonPion(double sqrt_s, const ResonanceParams& res, double m1, double m2,int channelIndex) const
{
    double pcm = p_CM(sqrt_s, m1, m2);

    if (pcm <= 1e-3 || sqrt_s <= m1 + m2 + 1e-3) 
    {
        return 0.0;
    }

    // Use the full width and apply isospin coupling factor via cgCoeff
    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.couplingFactor, res.L, res.m1, res.m2, res.LFactors);

    double spinFactor = res.twoSRplusOne / (2.0);

    double massDiff = res.mass - sqrt_s;
    double denom = ((massDiff * massDiff) + (Gamma_tot * Gamma_tot) *0.25);
    double sigmaGeV = 0.0; 
    if (channelIndex == 0) // proton-pion^- channel
    {
        sigmaGeV = (1.-res.cgCoeff) // prob is 1/3 for proton-pion^-. init function can also be used for proton-pion^0
                * spinFactor 
                * M_PI / (pcm * pcm) 
                * (res.couplingFactor[0]*Gamma_tot * Gamma_tot) / denom; 
    }

    if (channelIndex == 1) // proton-pion^0 channel
    {
        sigmaGeV = (res.cgCoeff) // prob is 2/3 for proton-pion^0. init function can also be used for proton-pion^0
                * spinFactor 
                * (M_PI / (pcm * pcm)) 
                * (res.couplingFactor[0]*Gamma_tot * Gamma_tot) / denom; // this is checked with Wigner-Eckart Theorem
    }
    return sigmaGeV;
}

double NPion::nPiMinus(double const &s)
{
    if (s <= 0.0) return 0.0;

    double sqrt_s = std::sqrt(s);
    if ( sqrt_s >= 2.) {
        return m_hpr1r2.xs_tot(HADRON_RESCATTERING::hpr1r2::nPiMinus, s);
    }

    else if (sqrt_s < 2)
    {
        if (sqrt_s <= (NEUTRON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
    
        for (const auto& res : pPiPlus_nPiMinus) 
        {
            //not sure about this, but assuming the same resonances as pPi+ channel.
        totalSigmaGeV2 += calculateResonanceSigma(sqrt_s, res, NEUTRON_MASS, PION_MASS);
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB;
        return sigma_mb;
    }
    return 0;
}

double NPion::pPiPlus(double const &s)
{
    if (s <= 0.0) return 0.0;
    msg_Out()<<"PpiPlus Called"<<std::endl;
    double sqrt_s = std::sqrt(s);
    if ( sqrt_s >= 2.) {
        return m_hpr1r2.xs_tot(HADRON_RESCATTERING::hpr1r2::pPiPlus, s);
    }
 
    else if (sqrt_s < 2)
    {
        if (sqrt_s <= (PROTON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
    
        for (const auto& res : pPiPlus_nPiMinus) 
        {
        totalSigmaGeV2 += calculateResonanceSigma(sqrt_s, res,PROTON_MASS,PION_MASS);
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB;
        return sigma_mb;
    }
    return 0;
}

double NPion::pPiMinus(const double &s)
{

    double sqrt_s = std::sqrt(s);

    if (s <= 0.0) return 0.0;

    if (sqrt_s < 2)
    {
        if (sqrt_s <= (PROTON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
        msg_Out()<<"PpiMinus Called"<<std::endl;
        for (const auto& res : pPiMZ_nPiPZ) 
        {
            totalSigmaGeV2 +=  calculateResonanceSigmaForProtonPion(sqrt_s, res,PROTON_MASS,PION_MASS,0); //0 for proton-pion^-
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB ;
        return sigma_mb;
    }
    return 0;
}

double NPion::pPiZero(const double &s)
{

    double sqrt_s = std::sqrt(s);

    if (s <= 0.0) return 0.0;

    if (sqrt_s < 2)
    {
        if (sqrt_s <= (PROTON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
    
        for (const auto& res : pPiMZ_nPiPZ) 
        {
            totalSigmaGeV2 +=  calculateResonanceSigmaForProtonPion(sqrt_s, res,PROTON_MASS,PION_MASS,1);//1 for proton-pion^0
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB ;
        return sigma_mb;
    }
    return 0;
}

double NPion::pPiZeroWignerEckart(const double &s) 
{
    double sigMinus = pPiMinus(s);
    double sigPlus = pPiPlus(s);

    return (1./2.)*sigPlus + (1./2.)*sigMinus; 
}


double NPion::calculateResonanceSigmaForNeutronPion(double sqrt_s, const ResonanceParams& res, double m1, double m2,int channelIndex) const
{
    double pcm = p_CM(sqrt_s, m1, m2);

    if (pcm <= 1e-3 || sqrt_s <= m1 + m2 + 1e-3) 
    {
        return 0.0;
    }

    // Use the full width and apply isospin coupling factor via cgCoeff
    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.couplingFactor, res.L, res.m1, res.m2, res.LFactors);

    double spinFactor = res.twoSRplusOne / (2.0);

    double massDiff = res.mass - sqrt_s;
    double denom = ((massDiff * massDiff) + (Gamma_tot * Gamma_tot) *0.25);
    double sigmaGeV = 0.0; 
    if (channelIndex == 0) // neutron-pion^+ channel
    {
        sigmaGeV = (1.-res.cgCoeff) // prob is 1/3 for delta, 2/3 or N: neutron-pion^+
                * spinFactor 
                * M_PI / (pcm * pcm) 
                * (res.couplingFactor[0]*Gamma_tot * Gamma_tot) / denom; 
    }

    if (channelIndex == 1) // neutron-pion^0 channel
    {
        sigmaGeV = (res.cgCoeff) // prob is 2/3 for delta, 1/3 or N: neutron-pion^0
                * spinFactor 
                * (M_PI / (pcm * pcm)) 
                * (res.couplingFactor[0]*Gamma_tot * Gamma_tot) / denom; // this is checked with Wigner-Eckart Theorem
    }
    return sigmaGeV;
}


double NPion::nPiPlus(const double & s)
{
    double sqrt_s = std::sqrt(s);

    if (s <= 0.0) return 0.0;
    msg_Out()<<"NPiPlus Called"<<std::endl;
    if (sqrt_s < 2)
    {
        if (sqrt_s <= (NEUTRON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
    
        for (const auto& res : pPiMZ_nPiPZ) 
        {
            totalSigmaGeV2 +=  calculateResonanceSigmaForProtonPion(sqrt_s, res,PROTON_MASS,PION_MASS,0); //0 for neutron-pion^+ channel
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB ;
        return sigma_mb;
    }
    return 0;
}

double NPion::nPiZero(const double & s)
{
    double sqrt_s = std::sqrt(s);

    if (s <= 0.0) return 0.0;

    if (sqrt_s < 2)
    {
        if (sqrt_s <= (NEUTRON_MASS + PION_MASS)) 
        {
            return 0.0;
        }

    // 2. Loop through all resonances and sum their cross-sections
        double totalSigmaGeV2 = 0.0;
    
        for (const auto& res : pPiMZ_nPiPZ) 
        {
            totalSigmaGeV2 +=  calculateResonanceSigmaForProtonPion(sqrt_s, res,PROTON_MASS,PION_MASS,1); //1 for neutron-pion^0 channel
        }
    
        double sigma_mb = totalSigmaGeV2 * GevToMB ;
        return sigma_mb;
    }
    return 0;
}

