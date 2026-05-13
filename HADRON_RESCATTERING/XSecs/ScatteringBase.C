#include "HADRON_RESCATTERING/XSecs/ScatteringBase.H"


 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_Resonances;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_Resonances_n_KaonMinus;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_PiPlusPiMinusRes;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_PiPlusPiZeroRes;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_PiZeroPiZeroRes;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::m_KPlusPiMinusRes;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::pPiPlus_nPiMinus;
 std::vector<ScatteringBase::ResonanceParams> ScatteringBase::pPiMZ_nPiPZ;
bool ScatteringBase::initialized = false; //global static data to ensure init runs only once during runtime. 

ScatteringBase::ScatteringBase()
{
    if ( !initialized ) 
    {
        Initialize();
        initialized = true;
    }
}

void ScatteringBase::Initialize()
{
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

    m_Resonances = 
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
    };


    m_PiPlusPiMinusRes = {
    {0.782, 0.008, {0.02, 0.09, 0.89}, 3.0, 1.0, 1.0/3.0, {PION_MASS, 0, PION_MASS}, {PION_MASS, PION_MASS, 2* PION_MASS}, {1, 1, 1}}, //w(770 ish)   
    {0.769, 0.151, {1.0}, 3.0, 1.0, 1.0, {PION_MASS}, {PION_MASS}, {1}}, //rho(770)
    {0.980, 0.1, {0.7, 0.3}, 1.0, 0., (2.0/3.0), {PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, CHARGEDKAON_MASS}, {0, 0}}, //f0(980)
    {1.275, 0.185, {0.84, 0.11, 0.05}, 5.0, 2.0, (2.0/3.0), {PION_MASS, 2 * PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, 2 * PION_MASS, CHARGEDKAON_MASS}, {2, 2, 2}}, //f2(1270)
    {1.370, 0.2, {0.1, 0.7, 0.2}, 1.0, 0., 2.0/3.0, {PION_MASS, 2*PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, 2*PION_MASS, CHARGEDKAON_MASS}, {0, 0, 0}}, // f0(1370)
    //{0.475, 0.5, {1.0}, 1.0, 0., (2.0/3.0), {PION_MASS}, {PION_MASS}, {0}}, //sigma or f0(500)
    {1.525, 0.076, {0.01, 0.1, 0.89}, 5.0, 2.0, 1.0/3.0, {PION_MASS, ETA_MASS, CHARGEDKAON_MASS}, {PION_MASS, ETA_MASS, CHARGEDKAON_MASS}, {2, 2, 2}}, // f'(1525)
    {1.465, 0.31, {0.5, 0.5}, 3.0, 1.0, 1.0, {PION_MASS, 2 * PION_MASS}, {PION_MASS, 2 * PION_MASS}, {1, 1}}, // rho(1465)
    {1.7, 0.235, {0.1, 0.2}, 3.0, 1.0, 1.0, {PION_MASS, 2 * PION_MASS, RHO_MASS}, {PION_MASS, 2 * PION_MASS, 2 * PION_MASS}, {1, 1, 1}}}; // rho (1700)

    m_PiPlusPiZeroRes = {
    {0.769, 0.151, {1.0}, 3.0, 1.0, 1.0, {PION_MASS}, {PION_MASS}, {1}}, //rho(770)
    {1.465, 0.31, {0.5, 0.5}, 3.0, 1.0, 1.0, {PION_MASS, 2 * PION_MASS}, {PION_MASS, 2 * PION_MASS}, {1, 1}}, // rho(1465)
    {1.7, 0.235, {0.1, 0.2, 0.7}, 3.0, 1.0, 1.0, {PION_MASS, 2 * PION_MASS, RHO_MASS}, {PION_MASS, 2 * PION_MASS, 2 * PION_MASS}, {1, 1, 1}}}; // rho (1700)

    m_PiZeroPiZeroRes = {
    {0.782, 0.008, {0.02, 0.09, 0.89}, 3.0, 1.0, 1.0/3.0, {PION_MASS, 0, PION_MASS}, {PION_MASS, PION_MASS, 2 * PION_MASS}, {1, 1, 1}}, //w(770 ish)
    {0.980, 0.1, {0.7, 0.3}, 1.0, 0., (1.0/3.0), {PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, CHARGEDKAON_MASS}, {0, 0}}, //f0(980)
    {1.370, 0.2, {0.1, 0.7, 0.2}, 1.0, 0., 1.0, {PION_MASS, 2*PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, 2*PION_MASS, CHARGEDKAON_MASS}, {0, 0, 0}}, // f0(1370) 
    {1.275, 0.185, {0.84, 0.11, 0.05}, 5.0, 2.0, (1.0), {PION_MASS, 2 * PION_MASS, CHARGEDKAON_MASS}, {PION_MASS, 2 * PION_MASS, CHARGEDKAON_MASS}, {2, 2, 2}}, //f2(1270) 
    {1.525, 0.076, {0.01, 0.1, 0.89}, 5.0, 2.0, 1.0, {PION_MASS, ETA_MASS, CHARGEDKAON_MASS}, {PION_MASS, ETA_MASS, CHARGEDKAON_MASS}, {2, 2, 2}}}; // f'(1525)  

    m_KPlusPiMinusRes = {
    {0.893, 0.050, {1.0}, 3.0, 1.0, 2.0/3.0, {CHARGEDKAON_MASS}, {PION_MASS}, {1}}, //K*(890)
    {1.429, 0.287, {1.0}, 1.0, 0., 2.0/3.0, {CHARGEDKAON_MASS}, {PION_MASS}, {0}}, //K*0
    {1.430, 0.10, {0.5, 0.25, 0.09, 0.03, 0.13}, 5.0, 2.0, 2.0/3.0, {CHARGEDKAON_MASS, KAONSTAR_MASS, CHARGEDKAON_MASS, CHARGEDKAON_MASS, CHARGEDKAON_MASS}, {PION_MASS, PION_MASS, RHO_MASS, OMEGA_MASS, 2 * PION_MASS}, {2, 2, 2, 2, 1}}, //K*2(1430)
    {1.41, 0.227, {0.3, 0.65, 0.05}, 3.0, 1.0, 2.0/3.0, {CHARGEDKAON_MASS, KAONSTAR_MASS, CHARGEDKAON_MASS}, {PION_MASS, PION_MASS, RHO_MASS}, {1, 1, 1}}, //K*(1410)
    {1.68, 0.323, {0.4, 0.3, 0.3}, 3.0, 1.0, 2.0/3.0, {CHARGEDKAON_MASS, KAONSTAR_MASS, CHARGEDKAON_MASS}, {PION_MASS, PION_MASS, RHO_MASS}, {1, 1, 1}}}; //K*(1680)
	    


pPiPlus_nPiMinus={
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
            {3, 3, 3, 3}}} ;// RES_1950 (Delta(1950) F37)
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

            pPiMZ_nPiPZ ={
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
            };

}


