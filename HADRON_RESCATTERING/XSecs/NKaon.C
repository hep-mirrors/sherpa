#include "NKaon.H"

NKaon::NKaon() : m_Resonances(
{
    // --- Lambda Resonances (I=0) ---
    {1.405, 0.050, {1.00}, 2.0, 0, 0.5, 
    {SIGMA_MASS}, {PION_MASS}, {0}}, //lambda 1405
    
    {1.520, 0.016, {0.45, 0.43, 0.11}, 4.0, 2, 0.5, 
    {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 0}}, //lambda 1520

    {1.600, 0.150, {0.35, 0.65}, 2.0, 1, 0.5, 
    {PROTON_MASS, SIGMA_MASS}, 
    {KAON_MASS, PION_MASS}, {1, 1}}, //lambda 1600

    {1.670, 0.035, {0.20, 0.50, 0.30}, 2.0, 0, 0.5, 
    {PROTON_MASS, SIGMA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, ETA_MASS}, {0, 0, 0}}, //lambda 1670

    {1.690, 0.060, {0.25, 0.45, 0.30}, 4.0, 2, 0.5, 
    {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 0}}, //lambda 1690

    {1.800, 0.300, {0.40, 0.20, 0.20, 0.20}, 2.0, 0, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS}, {0, 1, 0, 1}}, //lambda 1800

    {1.810, 0.150, {0.35, 0.45, 0.15, 0.05}, 2.0, 1, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS}, {1, 1, 1, 1}}, //lambda 1810

    {1.820, 0.080, {0.73, 0.16, 0.11}, 6.0, 2, 0.5, 
    {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 2}}, //lambda 1820

    {1.830, 0.095, {0.10, 0.70, 0.20}, 6.0, 2, 0.5, 
    {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 2}}, //lambda 1830

    {1.890, 0.100, {0.37, 0.21, 0.11, 0.31}, 4.0, 1, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS}, {1, 1, 1, 1}}, //lambda 1890

    {2.100, 0.200, {0.35, 0.20, 0.05, 0.30, 0.02, 0.08}, 8.0, 3, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS, ETA_MASS, 0.782}, {3,2,3,2,3,3}}, //lambda 2100

    {2.110, 0.200, {0.25, 0.45, 0.30}, 6.0, 3, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS}, {3, 2, 3}}, //lambda 2110

    // --- Sigma Resonances (I=1) ---
    {1.384, 0.036, {0.12, 0.88}, 4.0, 1, 0.5, {SIGMA_MASS, LAMBDA_MASS}, {PION_MASS, PION_MASS}, {1, 1}}, //sigma 1385

    {1.660, 0.100, {0.30, 0.35, 0.35}, 2.0, 1, 0.5, {PROTON_MASS, SIGMA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {1, 1, 1}}, //sigma 1660

    {1.670, 0.060, {0.15, 0.70, 0.15}, 4.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 2}}, //sigma 1670

    {1.750, 0.090, {0.40, 0.05, 0.55}, 2.0, 0, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_MASS}, 
    {KAON_MASS, PION_MASS, ETA_MASS}, {0, 0, 1}}, //sigma 1750

    {1.775, 0.120, {0.40, 0.04, 0.10, 0.23, 0.23}, 6.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405}, 
    {KAON_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS}, {2, 2, 0, 2, 2}}, //sigma 1775

    {1.915, 0.120, {0.15, 0.40, 0.05, 0.40}, 6.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS, PION_MASS}, {2,2,2,2}}, //sigma 1915

    {1.940, 0.220, {0.10, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15}, 4.0, 2, 0.5, {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405, DELTA1232_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS, KAON_MASS}, {2, 1, 2, 1, 2, 2, 2}},//sigma 1940
    
    {2.030, 0.180, {0.20, 0.04, 0.10, 0.10, 0.20, 0.18, 0.18}, 8.0, 3, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405, 1.232}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS, KAON_MASS}, {3, 2, 3, 2, 3,3,3}},//sigma 2030
}) 
{
    // Def. Constructor 
    m_Resonances_n_KaonMinus = { 
        {1.384, 0.036, {0.12, 0.88}, 4.0, 1, 0.5, {SIGMA_MASS, LAMBDA_MASS}, {PION_MASS, PION_MASS}, {1, 1}}, //sigma 1385

    {1.660, 0.100, {0.30, 0.35, 0.35}, 2.0, 1, 0.5, {PROTON_MASS, SIGMA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {1, 1, 1}}, //sigma 1660

    {1.670, 0.060, {0.15, 0.70, 0.15}, 4.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS}, {2, 2, 2}}, //sigma 1670

    {1.750, 0.090, {0.40, 0.05, 0.55}, 2.0, 0, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_MASS}, 
    {KAON_MASS, PION_MASS, ETA_MASS}, {0, 0, 1}}, //sigma 1750

    {1.775, 0.120, {0.40, 0.04, 0.10, 0.23, 0.23}, 6.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405}, 
    {KAON_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS}, {2, 2, 0, 2, 2}}, //sigma 1775

    {1.915, 0.120, {0.15, 0.40, 0.05, 0.40}, 6.0, 2, 0.5, {PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS}, 
    {KAON_MASS, PION_MASS, PION_MASS, PION_MASS}, {2,2,2,2}}, //sigma 1915

    {1.940, 0.220, {0.10, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15}, 4.0, 2, 0.5, {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405, DELTA1232_MASS}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS, KAON_MASS}, {2, 1, 2, 1, 2, 2, 2}},//sigma 1940
    
    {2.030, 0.180, {0.20, 0.04, 0.10, 0.10, 0.20, 0.18, 0.18}, 8.0, 3, 0.5, 
    {PROTON_MASS, PROTON_MASS, SIGMA_MASS, SIGMA_STAR_MASS, LAMBDA_MASS, 1.405, 1.232}, 
    {KAON_MASS, KSTAR_MASS, PION_MASS, PION_MASS, PION_MASS, PION_MASS, KAON_MASS}, {3, 2, 3, 2, 3,3,3}},//sigma 2030
    };
}


double NKaon::pKaonMinus(const double &s)
{
    double sqrt_s = std::sqrt(s);

    if (sqrt_s <= (PROTON_MASS + KAON_MASS)) 
    {
        return 0.0;
    }

    if (sqrt_s < 2.16)
    {

        double totalSigmaGeV = 0.0;
    
        for (const auto& res : m_Resonances) 
        {
            totalSigmaGeV +=  calculateResonanceScatteringXSec(sqrt_s, res, PROTON_MASS, KAON_MASS); 
        }
        
        double e0 = 1.433;

        if (sqrt_s < 1.4738188)
        totalSigmaGeV += 5.93763355 / pow((sqrt_s - 1.251377),2);
        else if (sqrt_s < 1.485215)
        totalSigmaGeV += -1.296457765e7 * pow((sqrt_s - e0),4)
                + 2.160975431e4 * pow((sqrt_s - e0),2) + 120.;
        else if (sqrt_s < 1.977)
        totalSigmaGeV += 3. + 1.0777e6 * exp(-6.4463 * sqrt_s)
                - 10. * exp(-pow((sqrt_s - 1.644),2) / 0.004)
                + 10. * exp(-pow((sqrt_s - 1.977),2) / 0.004);
        else
        totalSigmaGeV += 1.0777e6 * exp(-6.44463 * sqrt_s) + 12.5;

    
        double sigma_mb = totalSigmaGeV  ;
        msg_Out()<<sigma_mb<<std::endl;
        return sigma_mb;
    }
    
    else
    {
        return m_hpr1r2.xs_tot(hpr1r2::pKMinus,s);
    }
}

double NKaon::calculateResonanceScatteringXSec(double sqrt_s, const ResonanceParams& res, double m1, double m2) const
{
    double pcm = p_CM(sqrt_s, m1, m2);

    if (pcm <= 1e-3 || sqrt_s <= m1 + m2 + 1e-3) 
    {
        return 0.0;
    }
    double nominalPartialWidth = res.width0;
    double Gamma_tot = totalWidth(sqrt_s, res.mass, res.width0, res.branchingRatio, res.L, res.m1, res.m2, res.LFactors);

    if (Gamma_tot < 1e-9) 
    {
        return 0.0;
    }

    double spinFactor = res.twoSRplusOne / (2.0);

    double massDiff = res.mass - sqrt_s;
    double denom = ((massDiff * massDiff) + (Gamma_tot * Gamma_tot) *0.25);
    double sigma = 0.0; 

    sigma = GevToMB* ((res.probFactor)) 
                * spinFactor 
                * M_PI / (pcm * pcm) 
                * (res.branchingRatio[0]*Gamma_tot * Gamma_tot) / denom; 
    return sigma;/*  */
}


// Channel,Isospin Composition,Resonances to Include
// pK−,$\frac{1}{2},I=0\rangle + \frac{1}{2}
// nK−,$1,I=1\rangle$
// pK+,$1,"I=1, S=+1\rangle$"
// nK+,$\frac{1}{2},I=0\rangle + \frac{1}{2}


double NKaon::nKaonMinus(const double &s)
{
    double sqrt_s = std::sqrt(s);

    if (sqrt_s <= (PROTON_MASS + KAON_MASS)) 
    {
        return 0.0;
    }

    if (sqrt_s < 2.16)
    {

        double totalSigmaGeV = 0.0;
    
        for (const auto& res : m_Resonances_n_KaonMinus) 
        {
            totalSigmaGeV +=  calculateResonanceScatteringXSec(sqrt_s, res, NEUTRON_MASS, KAON_MASS); 
        }
        
        double e0 = 1.433;

        if (sqrt_s < 1.4738188)
        totalSigmaGeV += 5.93763355 / pow((sqrt_s - 1.251377),2);
        else if (sqrt_s < 1.485215)
        totalSigmaGeV += -1.296457765e7 * pow((sqrt_s - e0),4)
                + 2.160975431e4 * pow((sqrt_s - e0),2) + 120.;
        else if (sqrt_s < 1.977)
        totalSigmaGeV += 3. + 1.0777e6 * exp(-6.4463 * sqrt_s)
                - 10. * exp(-pow((sqrt_s - 1.644),2) / 0.004)
                + 10. * exp(-pow((sqrt_s - 1.977),2) / 0.004);
        else
        totalSigmaGeV += 1.0777e6 * exp(-6.44463 * sqrt_s) + 12.5;

    
        double sigma_mb = totalSigmaGeV  ;
        msg_Out()<<sigma_mb<<std::endl;
        return sigma_mb;
    }
    
    else
    {
        msg_Out()<<m_hpr1r2.xs_tot(hpr1r2::nKMinus,s)<<std::endl;

        return m_hpr1r2.xs_tot(hpr1r2::nKMinus,s);
    }
}