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


RotationAndBoost::RotationAndBoost() { reset(); }

void RotationAndBoost::reset() {
    for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
        M[i][j] = (i==j) ? 1. : 0.;
}

ATOOLS::Vec4D RotationAndBoost::Transform(const ATOOLS::Vec4D& p) const {
    double t=p[0], x=p[1], y=p[2], z=p[3];
    ATOOLS::Vec4D out;
    out[0] = M[0][0]*t + M[0][1]*x + M[0][2]*y + M[0][3]*z;
    out[1] = M[1][0]*t + M[1][1]*x + M[1][2]*y + M[1][3]*z;
    out[2] = M[2][0]*t + M[2][1]*x + M[2][2]*y + M[2][3]*z;
    out[3] = M[3][0]*t + M[3][1]*x + M[3][2]*y + M[3][3]*z;
    return out;
}

void RotationAndBoost::Boost(double betaX, double betaY, double betaZ) {
    double b2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
    double gm = 1.0 / sqrt(1.0 - b2);
    double gf = gm*gm / (1.0 + gm);

    double Mb[4][4] = {
        {gm,        gm*betaX,           gm*betaY,           gm*betaZ        },
        {gm*betaX,  1.+gf*betaX*betaX,  gf*betaX*betaY,     gf*betaX*betaZ  },
        {gm*betaY,  gf*betaY*betaX,     1.+gf*betaY*betaY,  gf*betaY*betaZ  },
        {gm*betaZ,  gf*betaZ*betaX,     gf*betaZ*betaY,     1.+gf*betaZ*betaZ}
    };

    double Mtmp[4][4];
    for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
        Mtmp[i][j] = M[i][j];

    for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
        M[i][j] = Mb[i][0]*Mtmp[0][j] + Mb[i][1]*Mtmp[1][j]
                + Mb[i][2]*Mtmp[2][j] + Mb[i][3]*Mtmp[3][j];
}

void RotationAndBoost::Rot(double theta, double phi) {
    double cthe=cos(theta), sthe=sin(theta);
    double cphi=cos(phi),   sphi=sin(phi);

    double Mr[4][4] = {
        {1., 0.,          0.,    0.         },
        {0., cthe*cphi,  -sphi,  sthe*cphi  },
        {0., cthe*sphi,   cphi,  sthe*sphi  },
        {0., -sthe,       0.,    cthe       }
    };

    double Mtmp[4][4];
    for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
        Mtmp[i][j] = M[i][j];

    for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
        M[i][j] = Mr[i][0]*Mtmp[0][j] + Mr[i][1]*Mtmp[1][j]
                + Mr[i][2]*Mtmp[2][j] + Mr[i][3]*Mtmp[3][j];
}

void RotationAndBoost::toCenterOfMomentumFrame(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2) {
    reset();

    ATOOLS::Vec4D pSum = p1 + p2;

    // boost to CM frame - boost dir and inv to get theta/phi
    double betaX = -pSum[1]/pSum[0];
    double betaY = -pSum[2]/pSum[0];
    double betaZ = -pSum[3]/pSum[0];

    // Boost dir and inv manually to get theta/phi with a fancy lambda function
    //capture clause--> capture everything by reference
    auto boostBack = [&](ATOOLS::Vec4D v) -> ATOOLS::Vec4D {
        double b2 = betaX*betaX + betaY*betaY + betaZ*betaZ;
        double gm = 1./sqrt(1.-b2);
        double p1 = betaX*v[1] + betaY*v[2] + betaZ*v[3];
        double p2 = gm*(gm*p1/(1.+gm) + v[0]);
        v[1] += p2*betaX; v[2] += p2*betaY; v[3] += p2*betaZ;
        v[0]  = gm*(v[0] + p1);
        return v;
    };

    ATOOLS::Vec4D dir = boostBack(p1);
    ATOOLS::Vec4D inv = boostBack(p2);

    double theta = atan2(sqrt(dir[1]*dir[1] + dir[2]*dir[2]), dir[3]);
    double phi   = atan2(dir[2], dir[1]);

    Boost(betaX, betaY, betaZ);
    Rot(0., -phi);
    Rot(-theta, phi);

    //  final equal-velocity boost if masses differ
    double sDir = p1[0]*p1[0] - p1[1]*p1[1] - p1[2]*p1[2] - p1[3]*p1[3];
    double sInv = p2[0]*p2[0] - p2[1]*p2[1] - p2[2]*p2[2] - p2[3]*p2[3];

    if (std::abs(sDir - sInv) > 1e-6 * (sDir + sInv)) {
        // Apply the two rotations to dir so it lies along z
        dir = Transform(p1);  // reuse the full matrix built so far

        double pAbs = sqrt(dir[1]*dir[1] + dir[2]*dir[2] + dir[3]*dir[3]);
        ATOOLS::Vec4D invT  = Transform(p2);

        double beta = (dir[0]*invT[0] - pAbs*pAbs - sqrt(sDir*sInv))
                    * (dir[0] + invT[0]) / (pAbs * (sDir - sInv));

        Boost(0., 0., beta);  // pure z-boost
    }
}






