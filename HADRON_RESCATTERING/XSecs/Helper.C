#include "HADRON_RESCATTERING/XSecs/Helper.H"


double p_CM(double m, double mA, double mB) {
    // Threshold check: if energy is below sum of masses, momentum is 0
    if (m <= mA + mB) return 0.0;

    double s = m * m;
    double term1 = s - (mA + mB) * (mA + mB);
    double term2 = s - (mA - mB) * (mA - mB);

    // Safety check to avoid sqrt of negative number due to precision
    if (term1 < 0.0 || term2 < 0.0) return 0.0;

    return std::sqrt(term1 * term2) / (2.0 * m);
}

// ============================================================================
// 2. Mass Dependent Partial Width
// Formula: Gamma(m) = Gamma0 * (m0/m) * (p/p0)^(2L+1) * 1.2 / (1 + 0.2*(p/p0)^2L)
// ============================================================================
double partialWidth_massDep(double m, double m0, double Gamma0, double l, double mA, double mB) {
    // 1. Threshold checks
    if (m <= mA + mB) return 0.0;
    
    // 2. Calculate momenta at current mass (p) and nominal mass (p0)
    double p_m   = p_CM(m,  mA, mB);
    double p_m0  = p_CM(m0, mA, mB);

    // If resonance is below its own decay threshold (rare/impossible for valid resonances), return 0
    if (p_m0 <= 0.0) return 0.0; // should not happen if m0 is above threshold
    if (p_m <= 0.0) return 0.0;  // below threshold -> no decay

    double pow_ratio = std::pow(p_m, 2*l + 1) / std::pow(p_m0, 2*l + 1);

    double damping = 1.2 / ( 1.0 + 0.2 * ( std::pow(p_m,2*l) / std::pow(p_m0,2*l)));

    double result = Gamma0 * (m0 / m) * pow_ratio * damping;
    return result;
}

// ============================================================================
// 3. Total Mass Dependent Width
// ============================================================================
double totalWidth(double m, double RES_M0, double GAMMA0, 
                  const std::vector<double> br, 
                  double L_entrance, 
                  const std::vector<double> m1, 
                  const std::vector<double> m2, 
                  const std::vector<int> LVector) 
{
    double total_gamma = 0.0;

    // Loop over all defined decay channels
    for (size_t i = 0; i < br.size(); ++i) 
    {
        // 1. Get Nominal Partial Width for this channel
        // Gamma0_partial = Total_Gamma0 * Branching_Ratio
        double nominal_partial = GAMMA0 * br[i];

        double current_L = LVector[i];
        double current_mA = m1[i];
        double current_mB = m2[i];
        // 2. Calculate Mass-Dependent Partial Width
        // We use the specific masses (m1[i], m2[i]) and Angular Momentum (LVector[i]) 
        // for this specific channel.
        double partial = partialWidth_massDep(m, 
                                                 RES_M0, 
                                                 nominal_partial, // Use Nominal PARTIAL Width
                                                 current_L, 
                                                 current_mA, 
                                                 current_mB);
        // double partial = partialWidth_massDep(m, RES_M0, nominal_partial, 
        //                                       (LVector[i]), 
        //                                       m1[i], m2[i]);
        
        // 3. Add to total
        total_gamma += partial;
    }

    return total_gamma;
}







