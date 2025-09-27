#include "SHERPA/SoftPhysics/Hadron_Init.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour_Tags.H"
#include "ATOOLS/Phys/KF_Table.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

void AddHadron(const kf_code &kfc, const double &mass, const double &radius,
               const double &width, const int icharge, const int spin,
               const bool majorana, const bool on, const int stable,
               const std::string &idname, const std::string &texname) {
  AddParticle(kfc, mass, radius, width, icharge, spin, on, stable, idname,
              texname);
  if (!majorana)
    s_kftable[kfc]->m_majorana = -1;
}

void AddHadron(const kf_code &kfc, const double &mass, const double &radius,
               const double &width, const int icharge, const int strong,
               const int spin, const bool majorana, const bool on,
               const int stable, const bool massive, const std::string &idname,
               const std::string &antiname, const std::string &texname,
               const std::string &antitexname) {
  AddParticle(kfc, mass, radius, width, icharge, strong, spin, majorana, on,
              stable, massive, idname, antiname, texname, antitexname);
}

void Hadron_Init::Init() {
  msg_Info()<<"Initializing hadron particle information ...\n";

  // ##########################################################################
  // ##########################################################################
  // Particle being used in SHERPA listed below
  // ##########################################################################
  //
  // ##########################################################################
  // Diquarks - not exactly particles but we need them ########################
  // ##########################################################################
  // ##########################################################################
  AddHadron(kf_dd_1,0.77133,0.,0,-2,-3,2,0,0,1,1, "dd_1","dd_1b","dd_1b","dd_1b");
  AddHadron(kf_ud_0,0.57933,0.,0,1,-3,0,0,0,1,1,  "ud_0","ud_0b","ud_0b","ud_0b");
  AddHadron(kf_ud_1,0.77133,0.,0,1,-3,2,0,0,1,1,  "ud_1","ud_1b","ud_1b","ud_1b");
  AddHadron(kf_uu_1,0.77133,0.,0,4,-3,2,0,0,1,1,  "uu_1","uu_1b","uu_1b","uu_1b");
  AddHadron(kf_sd_0,0.80473,0.,0,-2,-3,0,0,0,1,1, "sd_0","sd_0b","sd_0b","sd_0b");
  AddHadron(kf_sd_1,0.92953,0.,0,-2,-3,2,0,0,1,1, "sd_1","sd_1b","sd_1b","sd_1b");
  AddHadron(kf_su_0,0.80473,0.,0,1,-3,0,0,0,1,1,  "su_0","su_0b","su_0b","su_0b");
  AddHadron(kf_su_1,0.92953,0.,0,1,-3,2,0,0,1,1,  "su_1","su_1b","su_1b","su_1b");
  AddHadron(kf_ss_1,1.09361,0.,0,-2,-3,2,0,0,1,1, "ss_1","ss_1b","ss_1b","ss_1b");
  AddHadron(kf_cd_0,1.96908,0.,0,1,-3,0,0,0,1,1,  "cd_0","cd_0b","cd_0b","cd_0b");
  AddHadron(kf_cd_1,2.00808,0.,0,1,-3,2,0,0,1,1,  "cd_1","cd_1b","cd_1b","cd_1b");
  AddHadron(kf_cu_0,1.96908,0.,0,4,-3,0,0,0,1,1,  "cu_0","cu_0b","cu_0b","cu_0b");
  AddHadron(kf_cu_1,2.00808,0.,0,4,-3,2,0,0,1,1,  "cu_1","cu_1b","cu_1b","cu_1b");
  AddHadron(kf_cs_0,2.15432,0.,0,1,-3,0,0,0,1,1,  "cs_0","cs_0b","cs_0b","cs_0b");
  AddHadron(kf_cs_1,2.17967,0.,0,1,-3,2,0,0,1,1,  "cs_1","cs_1b","cs_1b","cs_1b");
  AddHadron(kf_cc_1,3.27531,0.,0,4,-3,2,0,0,1,1,  "cc_1","cc_1b","cc_1b","cc_1b");
  AddHadron(kf_bd_0,5.38897,0.,0,-2,-3,0,0,0,1,1, "bd_0","bd_0b","bd_0b","bd_0b");
  AddHadron(kf_bd_1,5.40145,0.,0,-2,-3,2,0,0,1,1, "bd_1","bd_1b","bd_1b","bd_1b");
  AddHadron(kf_bu_0,5.38897,0.,0,1,-3,0,0,0,1,1,  "bu_0","bu_0b","bu_0b","bu_0b");
  AddHadron(kf_bu_1,5.40145,0.,0,1,-3,2,0,0,1,1,  "bu_1","bu_1b","bu_1b","bu_1b");
  AddHadron(kf_bs_0,5.56725,0.,0,-2,-3,0,0,0,1,1, "bs_0","bs_0b","bs_0b","bs_0b");
  AddHadron(kf_bs_1,5.57536,0.,0,-2,-3,2,0,0,1,1, "bs_1","bs_1b","bs_1b","bs_1b");
  AddHadron(kf_bc_0,6.67143,0.,0,1,-3,0,0,0,1,1,  "bc_0","bc_0b","bc_0b","bc_0b");
  AddHadron(kf_bc_1,6.67397,0.,0,1,-3,2,0,0,1,1,  "bc_1","bc_1b","bc_1b","bc_1b");
  AddHadron(kf_bb_1,10.07354,0.,0,-2,-3,2,0,0,1,1,"bb_1","bb_1b","bb_1b","bb_1b");
  // ##################################################################################################
  // MESON MULTIPLETS:
  //   - will uniformly assume a radius of 0.65 fm
  // ##################################################################################################
  // Pseudoscalars   ##################################################################################
  AddHadron(kf_pi,           0.134976,0.65,7.8486e-09,0,0,false,1,0,"pi","pi");
  AddHadron(kf_pi_plus,      0.13957,0.65,2.5242e-17,3,0,true,1,1,  "pi+","pi^{+}");
  AddHadron(kf_eta,          0.5473,0.65,1.18e-06,0,0,false,1,0,    "eta","eta");
  AddHadron(kf_K,            0.49767,0.65,1.e-16,0,0,true,1,0,      "K","K");
  AddHadron(kf_K_L,          0.49767,0.65,1.273e-17,0,0,false,1,1,  "K(L)","K_{L}");
  AddHadron(kf_K_S,          0.49767,0.65,7.373e-15,0,0,false,1,0,  "K(S)","K_{S}");
  AddHadron(kf_K_plus,       0.493677,0.65,5.314e-17,3,0,true,1,1,  "K+","K^{+}");
  AddHadron(kf_eta_prime_958,0.95778,0.65,0.000203,0,0,false,1,0,   "eta'(958)","eta'(958)");
  AddHadron(kf_D_plus,       1.8693,0.65,6.23e-13,3,0,true,1,0,     "D+","D^{+}");
  AddHadron(kf_D,            1.8646,0.65,1.586e-12,0,0,true,1,0,    "D","D");
  AddHadron(kf_D_s_plus,     1.9685,0.65,1.41e-12,3,0,true,1,0,     "D(s)+","D(s)^{+}");
  AddHadron(kf_eta_c_1S,     2.9798,0.65,0.0132,0,0,false,1,0,      "eta(c)(1S)","eta_{c}(1S)");
  AddHadron(kf_B,            5.2792,0.65,4.22e-13,0,0,true,1,0,     "B","B");
  AddHadron(kf_B_plus,       5.2789,0.65,3.99e-13,3,0,true,1,0,     "B+","B^{+}");
  AddHadron(kf_B_s,          5.3693,0.65,4.27e-13,0,0,true,1,0,     "B(s)","B_{s}");
  AddHadron(kf_B_c,          6.4,0.65,1.43e-12,3,0,true,1,0,        "B(c)+","B_{c}^{+}");
  AddHadron(kf_eta_b,        9.4,0.65,0.050,0,0,false,1,0,          "eta(b)(1S)","eta_{b}(1S)");
  // Vectors         ##################################################################################
  AddHadron(kf_rho_770,         0.77,0.65,0.1507,0,2,false,1,0,     "rho(770)","rho(770)");
  AddHadron(kf_rho_770_plus,    0.77,0.65,0.1507,3,2,true,1,0,      "rho(770)+","rho^{+}(770)");
  AddHadron(kf_omega_782,       0.78194,0.65,0.00841,0,2,false,1,0, "omega(782)","omega(782)");
  AddHadron(kf_K_star_892,      0.8961,0.65,0.0505,0,2,true,1,0,    "K*(892)","K*(892)");
  AddHadron(kf_K_star_892_plus, 0.89166,0.65,0.0508,3,2,true,1,0,   "K*(892)+","K*^{+}(892)");
  AddHadron(kf_phi_1020,        1.01941,0.65,0.00443,0,2,false,1,0, "phi(1020)","phi(1020)");
  AddHadron(kf_D_star_2010_plus,2.01022,0.65,0.000083,3,2,true,1,0, "D*(2010)+","D*^{+}(2010)");
  AddHadron(kf_D_star_2007,     2.0067,0.65,0.001,0,2,true,1,0,     "D*(2007)","D*(2007)");
  AddHadron(kf_D_s_star_plus,   2.1124,0.65,0.001,3,2,true,1,0,     "D(s)*+","D_{s}*^{+}");
  AddHadron(kf_J_psi_1S,        3.09688,0.65,8.7e-05,0,2,false,1,0, "J/psi(1S)","J/psi(1S)");
  AddHadron(kf_B_star,          5.3249,0.65,0.0001,0,2,true,1,0,    "B*","B*");
  AddHadron(kf_B_star_plus,     5.3249,0.65,0.0001,3,2,true,1,0,    "B*+","B*^{+}");
  AddHadron(kf_B_s_star,        5.41630,0.65,0.0001,0,2,true,1,0,   "B(s)*","B_{s}*");
  AddHadron(kf_B_c_star,        6.602,0.65,0.0001,3,2,true,1,0,     "B(c)*+","B_{c}*^{+}");
  AddHadron(kf_Upsilon_1S,      9.46037,0.65,5.25e-05,0,2,false,1,0,"Upsilon(1S)","Upsilon(1S)");
  // Tensors 2       ##################################################################################
  AddHadron(kf_a_2_1320,          1.3181,0.65,0.107,0,4,false,1,0, "a(2)(1320)","a_{2}(1320)");
  AddHadron(kf_a_2_1320_plus,     1.3181,0.65,0.107,3,4,true,1,0,  "a(2)(1320)+","a_{2}^{+}(1320)");
  AddHadron(kf_f_2_1270,          1.275,0.65,0.1855,0,4,false,1,0, "f(2)(1270)","f_{2}(1270)");
  AddHadron(kf_K_2_star_1430,     1.4324,0.65,0.109,0,4,true,1,0,  "K(2)*(1430)","K_{2}*(1430)");
  AddHadron(kf_K_2_star_1430_plus,1.4256,0.65,0.0985,3,4,true,1,0, "K(2)*(1430)+","K_{2}*^{+}(1430)");
  AddHadron(kf_f_2_prime_1525,    1.525,0.65,0.076,0,4,false,1,0,  "f(2)'(1525)","f_{2}'(1525)");
  AddHadron(kf_D_2_star_2460_plus,2.459,0.65,0.025,3,4,true,1,0,   "D(2)*(2460)+","D_{2}*^{+}(2460)");
  AddHadron(kf_D_2_star_2460,     2.4589,0.65,0.023,0,4,true,1,0,  "D(2)*(2460)","D_{2}*(2460)");
  AddHadron(kf_D_s2_star_2573,    2.5735,0.65,0.015,3,4,true,1,0,  "D(s2)*(2573)+","D_{s2}*^{+}(2573)");
  AddHadron(kf_chi_c2_1P,         3.55617,0.65,0.002,0,4,false,1,0,"chi(c2)(1P)","chi_{c2}(1P)");
  AddHadron(kf_B_2_star,          5.83,0.65,0.02,0,4,true,1,0,     "B(2)*","B_{2}*");
  AddHadron(kf_B_2_star_plus,     5.83,0.65,0.02,3,4,true,1,0,     "B(2)*+","B_{2}*^{+}");
  AddHadron(kf_B_s2_star,         5.8397,0.65,0.02,0,4,true,1,0,   "B(s2)*","B_{s2}*");
  AddHadron(kf_B_c2_star,         7.35,0.65,0.02,3,4,true,1,0,     "B(c2)*+","B_{c2}*^{+}");
  AddHadron(kf_chi_b2_1P,         9.9132,0.65,0.001,0,4,false,1,0, "chi(b2)(1P)","chi_{b2}(1P)");
  // Scalars         ##################################################################################
  AddHadron(kf_a_0_1450,          1.474,0.65,0.265,0,0,false,1,0, "a(0)(1450)","a_{0}(1450)");
  AddHadron(kf_a_0_1450_plus,     1.474,0.65,0.265,3,0,true,1,0,  "a(0)(1450)+","a_{0}^{+}(1450)");
  AddHadron(kf_f_0_1370,          1.4,0.65,0.5,0,0,false,1,0,     "f(0)(1370)","f_{0}(1370)");
  AddHadron(kf_K_0_star_1430,     1.429,0.65,0.287,0,0,true,1,0,  "K(0)*(1430)","K_{0}*(1430)");
  AddHadron(kf_K_0_star_1430_plus,1.429,0.65,0.287,3,0,true,1,0,  "K(0)*(1430)+","K_{0}*^{+}(1430)");
  AddHadron(kf_f_0_1710,          1.7,0.65,0.5,0,0,false,1,0,     "f(0)(1710)","f_{0}(1710)");
  AddHadron(kf_D_0_star_plus,     2.351,0.65,0.23,3,0,true,1,0,   "D(0)*(2400)+","D_{0}*^{+}(2400)");
  AddHadron(kf_D_0_star,          2.318,0.65,0.267,0,0,true,1,0,  "D(0)*(2400)","D_{0}*(2400)");
  AddHadron(kf_D_s0_star,         2.3177,0.65,0.03,3,0,true,1,0,  "D(s0)*(2317)+","D_{s0}*^{+}(2317)");
  AddHadron(kf_chi_c0_1P,         3.4173,0.65,0.014,0,0,false,1,0,"chi(c0)(1P)","chi_{c0}(1P)");
  AddHadron(kf_B_0_star,          5.68,0.65,0.05,0,0,true,1,0,    "B(0)*","B_{0}*");
  AddHadron(kf_B_0_star_plus,     5.68,0.65,0.05,3,0,true,1,0,    "B(0)*+","B_{0}*^{+}");
  AddHadron(kf_B_s0_star,         5.92,0.65,0.05,0,0,true,1,0,    "B(s0)*","B_{s0}*");
  AddHadron(kf_B_c0_star,         7.25,0.65,0.05,3,0,true,1,0,    "B(c0)*+","B_{c0}*^{+}");
  AddHadron(kf_chi_b0_1P,         9.8598,0.65,0.050,0,0,false,1,0,"chi(b0)(1P)","chi_{b0}(1P)");
  // Axial vectors   ##################################################################################
  AddHadron(kf_b_1_1235,      1.2295,0.65,0.142,0,2,false,1,0, "b(1)(1235)","b_{1}(1235)");
  AddHadron(kf_b_1_1235_plus, 1.2295,0.65,0.142,3,2,true,1,0,  "b(1)(1235)+","b_{1}^{+}(1235)");
  AddHadron(kf_h_1_1170,      1.17,0.65,0.36,0,2,false,1,0,    "h(1)(1170)","h_{1}(1170)");
  AddHadron(kf_K_1_1270,      1.272,0.65,0.09,0,2,true,1,0,    "K(1)(1270)","K_{1}(1270)");
  AddHadron(kf_K_1_1270_plus, 1.272,0.65,0.09,3,2,true,1,0,    "K(1)(1270)+","K_{1}^{+}(1270)");
  AddHadron(kf_h_1_1380,      1.386,0.65,0.091,0,2,false,1,0,  "h(1)(1380)","h_{1}(1380)");
  AddHadron(kf_D_1_2420_plus, 2.4232,0.65,0.025,3,2,true,1,0,  "D(1)(2420)+","D_{1}^{+}(2420)");
  AddHadron(kf_D_1_2420,      2.4208,0.65,0.0317,0,2,true,1,0, "D(1)(2420)","D_{1}(2420)");
  AddHadron(kf_D_s1_2536_plus,2.5351,0.65,0.00092,3,2,true,1,0,"D(s1)(2536)+","D_{s1}^{+}(2536)");
  AddHadron(kf_h_c1,          3.46,0.65,0.01,0,2,false,1,0,    "h(c)(1P)","h_{c}(1P)");
  AddHadron(kf_B_1,           5.73,0.65,0.05,0,2,true,1,0,     "B(1)(L)","B_{1}(L)");
  AddHadron(kf_B_1_plus,      5.73,0.65,0.05,3,2,true,1,0,     "B(1)(L)+","B_{1}^{+}(L)");
  AddHadron(kf_B_s1,          5.97,0.65,0.05,0,2,true,1,0,     "B(s1)(L)","B_{s1}(L)");
  AddHadron(kf_B_c1,          7.3,0.65,0.05,3,2,true,1,0,      "B(c1)(L)+","B_{c1}^{+}(L)");
  AddHadron(kf_h_b1,          9.875,0.65,0.01,0,2,false,1,0,   "h(b)(1P)","h_{b}(1P)");
  // Vectors         ##################################################################################
  AddHadron(kf_a_1_1260,     1.23,0.65,0.400,0,2,false,1,0,     "a(1)(1260)","a_{1}(1260)");
  AddHadron(kf_a_1_1260_plus,1.23,0.65,0.400,3,2,true,1,0,      "a(1)(1260)+","a_{1}^{+}(1260)");
  AddHadron(kf_f_1_1285,     1.2819,0.65,0.024,0,2,false,1,0,   "f(1)(1285)","f_{1}(1285)");
  AddHadron(kf_K_1_1400,     1.402,0.65,0.174,0,2,true,1,0,     "K(1)(1400)","K_{1}(1400)");
  AddHadron(kf_K_1_1400_plus,1.402,0.65,0.174,3,2,true,1,0,     "K(1)(1400)+","K_{1}^{+}(1400)");
  AddHadron(kf_f_1_1420,     1.4262,0.65,0.055,0,2,false,1,0,   "f(1)(1420)","f_{1}(1420)");
  AddHadron(kf_D_1_H_plus,   2.427,0.65,0.384,3,2,true,1,0,     "D(1)(H)+","D_{1}^{+}(H)");
  AddHadron(kf_D_1_H,        2.427,0.65,0.384,0,2,true,1,0,     "D(1)(2430)","D_{1}(2430)");
  AddHadron(kf_D_s1_H,       2.4595,0.65,0.003,3,2,true,1,0,    "D(s1)(2460)+","D_{s1}^{+}(2460)");
  AddHadron(kf_chi_c1_1P,    3.51053,0.65,0.00088,0,2,false,1,0,"chi(c1)(1P)","chi_{c1}(1P)");
  AddHadron(kf_B_1_H,        5.78,0.65,0.05,0,2,true,1,0,       "B(1)(H)","B_{1}(H)");
  AddHadron(kf_B_1_H_plus,   5.78,0.65,0.05,3,2,true,1,0,       "B(1)(H)+","B_{1}^{+}(H)");
  AddHadron(kf_B_s1_H,       6.02,0.65,0.05,0,2,true,1,0,       "B(s1)(H)","B_{s1}(H)");
  AddHadron(kf_B_c1_H,       7.3,0.65,0.05,3,2,true,1,0,        "B(c1)(H)+","B_{c1}^{+}(H)");
  AddHadron(kf_chi_b1_1P,    9.8919,0.65,0.001,0,2,false,1,0,   "chi(b1)(1P)","chi_{b1}(1P)");
  // ##################################################################################################
  // BARYON MULTIPLETS
  //   - will uniformly assume a radius of 0.8783 fm
  // ##################################################################################################
  // Nucleons (octet) #################################################################################
  AddHadron(kf_n,                     0.939566,0.8783,7.424e-28,0,1,true,1,1,"n","n");
  AddHadron(kf_p_plus,                0.938272,0.8783,0,3,1,true,1,1,        "P+","P^{+}");
  AddHadron(kf_Sigma_minus,           1.19745,0.8783,4.45e-15,-3,1,true,1,0, "Sigma-","\\Sigma^{-}");
  AddHadron(kf_Sigma,                 1.19264,0.8783,8.9e-06,0,1,true,1,0,   "Sigma","\\Sigma");
  AddHadron(kf_Sigma_plus,            1.18937,0.8783,8.24e-15,3,1,true,1,0,  "Sigma+","\\Sigma^{+}");
  AddHadron(kf_Lambda,                1.11568,0.8783,2.501e-15,0,1,true,1,0, "Lambda","\\Lambda");
  AddHadron(kf_Xi_minus,              1.32132,0.8783,4.02e-15,-3,1,true,1,0, "Xi-","\\Xi^{-}");
  AddHadron(kf_Xi,                    1.3149,0.8783,2.27e-15,0,1,true,1,0,   "Xi","\\Xi");
  AddHadron(kf_Sigma_c_2455,          2.4522,0.8783,0.0022,0,1,true,1,0,     "Sigma(c)(2455)","\\Sigma_{c}(2455)");
  AddHadron(kf_Sigma_c_2455_plus,     2.4536,0.8783,0.0023,3,1,true,1,0,     "Sigma(c)(2455)+","\\Sigma_{c}^{+}(2455)");
  AddHadron(kf_Sigma_c_2455_plus_plus,2.4528,0.8783,0.0023,6,1,true,1,0,     "Sigma(c)(2455)++","\\Sigma_{c}^{++}(2455)");
  AddHadron(kf_Lambda_c_plus,         2.2849,0.8783,3.19e-12,3,1,true,1,0,   "Lambda(c)+","\\Lambda_{c}^{+}");
  AddHadron(kf_Xi_c_2466,             2.4703,0.8783,5.875e-12,0,1,true,1,0,  "Xi(c)","\\Xi_{c}");
  AddHadron(kf_Xi_c_2466_plus,        2.4656,0.8783,1.489e-12,3,1,true,1,0,  "Xi(c)+","\\Xi_{c}^{+}");
  AddHadron(kf_Xi_c_2574,             2.575,0.8783,0.001,0,1,true,1,0,       "Xi(c)'","\\Xi'_{c}");
  AddHadron(kf_Xi_c_2574_plus,        2.578,0.8783,0.001,3,1,true,1,0,       "Xi(c)'+","\\Xi'_{c}^{+}");
  AddHadron(kf_Omega_c_0,             2.704,0.8783,1.02e-11,0,1,true,1,0,    "Omega(c)","\\Omega_{c}");
  AddHadron(kf_Sigma_b_5820_minus,    5.810,0.8783,0.001,-3,1,true,1,0,      "Sigma(b)-","\\Sigma_{b}^{-}");
  AddHadron(kf_Sigma_b_5820,          5.810,0.8783,0.001,0,1,true,1,0,       "Sigma(b)","\\Sigma_{b}");
  AddHadron(kf_Sigma_b_5820_plus,     5.810,0.8783,0.001,3,1,true,1,0,       "Sigma(b)+","\\Sigma_{b}^{+}");
  AddHadron(kf_Lambda_b,              5.624,0.8783,5.31e-13,0,1,true,1,0,    "Lambda(b)","\\Lambda_{b}");
  AddHadron(kf_Xi_b_5840,             5.790,0.8783,1.e-13,-3,1,true,1,0,     "Xi(b)-","\\Xi_{b}^{-}");
  AddHadron(kf_Xi_b_5840_minus,       5.790,0.8783,1.e-13,0,1,true,1,0,      "Xi(b)","\\Xi_{b}");
  AddHadron(kf_Xi_b_5960,             5.890,0.8783,1.e-13,-3,1,true,1,0,     "Xi(b)'-_fict","\\Xi'_{b}_fict^{-}");
  AddHadron(kf_Xi_b_5960_minus,       5.890,0.8783,1.e-13,0,1,true,1,0,      "Xi(b)'_fict","\\Xi'_{b}_fict");
  AddHadron(kf_Omega_b_0,             6.071,0.8783,1.e-13,-3,1,true,1,0,     "Omega(b)-","\\Omega_{b}^{-}");
  // Deltas (decuplet) ################################################################################
  AddHadron(kf_Delta_1232_minus,      1.232,0.8783,0.12,-3,3,true,1,0,      "Delta(1232)-","\\Delta-(1232)");
  AddHadron(kf_Delta_1232,            1.232,0.8783,0.12,0,3,true,1,0,       "Delta(1232)","\\Delta(1232)");
  AddHadron(kf_Delta_1232_plus,       1.232,0.8783,0.12,3,3,true,1,0,       "Delta(1232)+","\\Delta^{+}(1232)");
  AddHadron(kf_Delta_1232_plus_plus,  1.232,0.8783,0.12,6,3,true,1,0,       "Delta(1232)++","\\Delta^{++}(1232)");
  AddHadron(kf_Sigma_1385_minus,      1.3872,0.8783,0.0394,-3,3,true,1,0,   "Sigma(1385)-","\\Sigma-(1385)");
  AddHadron(kf_Sigma_1385,            1.3837,0.8783,0.036,0,3,true,1,0,     "Sigma(1385)","\\Sigma(1385)");
  AddHadron(kf_Sigma_1385_plus,       1.3828,0.8783,0.0358,3,3,true,1,0,    "Sigma(1385)+","\\Sigma^{+}(1385)");
  AddHadron(kf_Xi_1530_minus,         1.535,0.8783,0.0099,-3,3,true,1,0,    "Xi(1530)-","\\Xi-(1530)");
  AddHadron(kf_Xi_1530,               1.5318,0.8783,0.0091,0,3,true,1,0,    "Xi(1530)","\\Xi(1530)");
  AddHadron(kf_Omega_minus,           1.67245,0.8783,8.01e-15,-3,3,true,1,0,"Omega-","\\Omega^{-}");
  AddHadron(kf_Sigma_c_2520,          2.5175,0.8783,0.0150,0,3,true,1,0,    "Sigma(c)(2520)","\\Sigma_{c}(2520)");
  AddHadron(kf_Sigma_c_2520_plus,     2.5159,0.8783,0.0150,3,3,true,1,0,    "Sigma(c)(2520)+","\\Sigma_{c}^{+}(2520)");
  AddHadron(kf_Sigma_c_2520_plus_plus,2.5194,0.8783,0.0150,6,3,true,1,0,    "Sigma(c)(2520)++","\\Sigma_{c}^{++}(2520)");
  AddHadron(kf_Xi_c_2645,             2.645,0.8783,0.003,0,3,true,1,0,      "Xi(c)*","\\Xi*_{c}");
  AddHadron(kf_Xi_c_2645_plus,        2.645,0.8783,0.003,3,3,true,1,0,      "Xi(c)*+","\\Xi*_{c}^{+}");
  AddHadron(kf_Omega_c_star,          2.8,0.8783,0.001,0,3,true,1,0,        "Omega(c)*","\\Omega*_{c}");
  AddHadron(kf_Sigma_b_5840_minus,    5.829,0.8783,0.01,-3,3,true,1,0,      "Sigma(b)*-","\\Sigma*_{b}^{-}");
  AddHadron(kf_Sigma_b_5840,          5.829,0.8783,0.01,0,3,true,1,0,       "Sigma(b)*","\\Sigma*_{b}");
  AddHadron(kf_Sigma_b_5840_plus,     5.829,0.8783,0.01,3,3,true,1,0,       "Sigma(b)*+","\\Sigma*_{b}^{+}");
  AddHadron(kf_Xi_b_5940,             5.930,0.8783,0.001,-3,3,true,1,0,     "Xi(b)*-_fict","\\Xi*_{b}_fict^{-}");
  AddHadron(kf_Xi_b_5940_minus,       5.930,0.8783,0.001,0,3,true,1,0,      "Xi(b)*_fict","\\Xi*_{b}_fict");
  AddHadron(kf_Omega_b_star,          6.090,0.8783,0.0003,-3,3,true,1,0,    "Omega(b)*-_fict","\\Omega*_{b}_fict^{-}");
  // Nucleons (octet) - L_N = 1 ##############################################
  // careful - we will have to add some heavy states here!
  AddHadron(kf_N_1535,          1.535,0.8783,0.15,0,1,true,1,0,   "N(1535)","N(1535)");
  AddHadron(kf_N_1535_plus,     1.535,0.8783,0.15,3,1,true,1,0,   "N(1535)+","N^{+}(1535)");
  AddHadron(kf_Sigma_1620_minus,1.62,0.8783,0.09,-3,1,true,1,0,   "Sigma(1620)-","\\Sigma-(1620)");
  AddHadron(kf_Sigma_1620,      1.62,0.8783,0.09,0,1,true,1,0,    "Sigma(1620)","\\Sigma(1620)");
  AddHadron(kf_Sigma_1620_plus, 1.62,0.8783,0.09,3,1,true,1,0,    "Sigma(1620)+","\\Sigma^{+}(1620)");
  AddHadron(kf_Lambda_1670,     1.67,0.8783,0.06,0,1,true,1,0,    "Lambda(1670)","\\Lambda(1670)");
  AddHadron(kf_Xi_1750_minus,   1.75,0.8783,0.09,-3,1,true,1,0,   "Xi(1750)-","\\Xi^{-}(1750)");
  AddHadron(kf_Xi_1750,         1.75,0.8783,0.09,0,1,true,1,0,    "Xi(1750)","\\Xi(1750)");
  AddHadron(kf_Lambda_1405,     1.407,0.8783,0.05,0,1,true,1,0,   "Lambda(1405)","\\Lambda(1405)");
  AddHadron(102142,             2.5954,0.8783,0.0036,3,1,true,1,0,"Lambda(c)(2595)+","\\Lambda_{c}(2595)^{+}");
  // Nucleons (octet) - the Roper resonance - L_N = 2 #########################
  // careful - we will have to add some heavy states here!
  AddHadron(kf_N_1440,          1.44,0.8783,0.35,0,1,true,1,0,   "N(1440)","N(1440)");
  AddHadron(kf_N_1440_plus,     1.44,0.8783,0.35,3,1,true,1,0,   "N(1440)+","N^{+}(1440)");
  AddHadron(kf_Sigma_1660_minus,1.66,0.8783,0.1,-3,1,true,1,0,   "Sigma(1660)-","\\Sigma-(1660)");
  AddHadron(kf_Sigma_1660,      1.66,0.8783,0.1,0,1,true,1,0,    "Sigma(1660)","\\Sigma(1660)");
  AddHadron(kf_Sigma_1660_plus, 1.66,0.8783,0.1,3,1,true,1,0,    "Sigma(1660)+","\\Sigma^{+}(1660)");
  AddHadron(kf_Lambda_1600,     1.6,0.8783,0.15,0,1,true,1,0,    "Lambda(1600)","\\Lambda(1600)");
  AddHadron(kf_Xi_1690_minus,   1.696,0.8783,0.010,-3,1,true,1,0,"Xi(1690)-","\\Xi-(1690)");
  AddHadron(kf_Xi_1690,         1.696,0.8783,0.010,0,1,true,1,0, "Xi(1690)","\\Xi(1690)");
  // Nucleons (octet) - L_N = 1_1 -- plus "singlet heavies" ###########################################
  AddHadron(kf_N_1520,          1.52,0.8783,0.12,0,3,true,1,0,    "N(1520)","N(1520)");
  AddHadron(kf_N_1520_plus,     1.52,0.8783,0.12,3,3,true,1,0,    "N(1520)+","N^{+}(1520)");
  AddHadron(kf_Sigma_1670_minus,1.67,0.8783,0.06,-3,3,true,1,0,   "Sigma(1670)-","\\Sigma-(1670)");
  AddHadron(kf_Sigma_1670,      1.67,0.8783,0.06,0,3,true,1,0,    "Sigma(1670)","\\Sigma(1670)");
  AddHadron(kf_Sigma_1670_plus, 1.67,0.8783,0.06,3,3,true,1,0,    "Sigma(1670)+","\\Sigma^{+}(1670)");
  AddHadron(kf_Lambda_1690,     1.69,0.8783,0.06,0,3,true,1,0,    "Lambda(1690)","\\Lambda(1690)");
  AddHadron(kf_Xi_1820_minus,   1.823,0.8783,0.024,-3,3,true,1,0, "Xi(1820)-","\\Xi-(1820)");
  AddHadron(kf_Xi_1820,         1.823,0.8783,0.024,0,3,true,1,0,  "Xi(1820)","\\Xi(1820)");
  AddHadron(kf_Lambda_1520,     1.5195,0.8783,0.0156,0,3,true,1,0,"Lambda(1520)","\\Lambda(1520)");
  AddHadron(102144,             2.625,0.8783,0.002,3,3,true,1,0,  "Lambda(c)(2625)+","\\Lambda_{c}^{+}(2625)");
  AddHadron(104314,             2.815,0.8783,0.002,0,3,true,1,0,  "Xi(c)(2815)","\\Xi_{c}(2815)");
  AddHadron(104324,             2.815,0.8783,0.002,3,3,true,1,0,  "Xi(c)(2815)+","\\Xi_{c}^{+}(2815)");
  AddHadron(102154,             5.91,0.8783,0.002,0,3,true,1,0,   "Lambda(b)(5910)","\\Lambda_{b}(5910)");
  AddHadron(kf_N_1710,          1.71,0.8783,0.12,0,1,true,1,0,    "N(1710)","N(1710)");
  AddHadron(kf_N_1710_plus,     1.71,0.8783,0.12,3,1,true,1,0,    "N(1710)+","N(1710)+");
  // #########################################################################
  // Obsolete multiple heavy baryons #########################################
  // they will not be produced in our code (we have no heavy di-quarks
  // that we can produce yet)
  // #########################################################################
  AddHadron(4412,3.59798,0.8783,0.,3,1,true,1,0,"Xi(cc)+","\\Xi(cc)^{+}");
  AddHadron(4414,3.65648,0.8783,0.,3,3,true,1,0,"Xi(cc)*+","\\Xi(cc)*^{+}");
  AddHadron(4422,3.59798,0.8783,0.,6,1,true,1,0,"Xi(cc)++","\\Xi(cc)^{++}");
  AddHadron(4424,3.65648,0.8783,0.,6,3,true,1,0,"Xi(cc)*++","\\Xi(cc)*^{++}");
  AddHadron(4432,3.78663,0.8783,0.,3,1,true,1,0,"Omega(cc)+","\\Omega(cc)^{+}");
  AddHadron(4434,3.82466,0.8783,0.,3,3,true,1,0,"Omega(cc)*+","\\Omega(cc)*^{+}");
  AddHadron(4444,4.91594,0.8783,0.,6,3,true,1,0,"Omega(ccc)*++","\\Omega(ccc)*^{++}");
  AddHadron(5142,7.00575,0.8783,0.,0,1,true,1,0,"Xi(bc)","\\Xi(bc)");
  AddHadron(5242,7.00575,0.8783,0.,3,1,true,1,0,"Xi(bc)+","\\Xi(bc)^{+}");
  AddHadron(5342,7.19099,0.8783,0.,0,1,true,1,0,"Omega(bc)","\\Omega(bc)");
  AddHadron(5412,7.03724,0.8783,0.,0,1,true,1,0,"Xi(bc)'","\\Xi(bc)'");
  AddHadron(5414,7.0485,0.8783,0.,0,3,true,1,0,"Xi(bc)*","\\Xi(bc)*");
  AddHadron(5422,7.03724,0.8783,0.,3,1,true,1,0,"Xi(bc)'+","\\Xi(bc)'^{+}");
  AddHadron(5424,7.0485,0.8783,0.,3,3,true,1,0,"Xi(bc)*+","\\Xi(bc)*^{+}");
  AddHadron(5432,7.21101,0.8783,0.,0,1,true,1,0,"Omega(bc)'","\\Omega(bc)'");
  AddHadron(5434,7.219,0.8783,0.,0,3,true,1,0,"Omega(bc)*","\\Omega(bc)*");
  AddHadron(5442,8.30945,0.8783,0.,3,1,true,1,0,"Omega(bcc)+","\\Omega(bcc)^{+}");
  AddHadron(5444,8.31325,0.8783,0.,3,3,true,1,0,"Omega(bcc)*+","\\Omega(bcc)*^{+}");
  AddHadron(5512,10.42272,0.8783,0.,-3,1,true,1,0,"Xi(bb)-","\\Xi(bb)^{-}");
  AddHadron(5514,10.44144,0.8783,0.,-3,3,true,1,0,"Xi(bb)*-","\\Xi(bb)*^{-}");
  AddHadron(5522,10.42272,0.8783,0.,0,1,true,1,0,"Xi(bb)","\\Xi(bb)");
  AddHadron(5524,10.44144,0.8783,0.,0,3,true,1,0,"Xi(bb)*","\\Xi(bb)*");
  AddHadron(5532,10.60209,0.8783,0.,-3,1,true,1,0,"Omega(bb)-","\\Omega(bb)^{-}");
  AddHadron(5534,10.61426,0.8783,0.,-3,3,true,1,0,"Omega(bb)*-","\\Omega(bb)*^{-}");
  AddHadron(5542,11.70767,0.8783,0.,0,1,true,1,0,"Omega(bbc)","\\Omega(bbc)");
  AddHadron(5544,11.71147,0.8783,0.,0,3,true,1,0,"Omega(bbc)*","\\Omega(bbc)*");
  AddHadron(5554,15.11061,0.8783,0.,-3,3,true,1,0,"Omega(bbb)*-","\\Omega(bbb)*^{-}");

  // ##########################################################################
  // ##########################################################################
  // Particles NOT being used in SHERPA listed below
  // ##########################################################################
  // ##########################################################################
  // The following hadron multiplets are incomplete and sometimes dodgy.
  // To play it safe we will not encode any multiplet weights, which means
  // that even if they are switched on, they will be ignored in the
  // Multiplet_Constructor
  //
  // My guess is that we will have to tidy up ... and kick out everything
  // we do not need in the decays.
  //
  // Meson states first
  // ##########################################################################
  // Tensors 3       ##########################################################
  // heavy ones missing - the rho(1690) could be useful for
  // tau/D/B decays with hadrons.
  AddHadron(kf_rho_3_1690,1.691,0.65,0.16,0,6,false,1,0,"rho(3)(1690)","rho_{3}(1690)");
  AddHadron(217,1.691,0.65,0.16,3,6,true,1,0,"rho(3)(1690)+","rho_{3}^{+}(1690)");
  AddHadron(kf_omega_3_1670,1.667,0.65,0.168,0,6,false,1,0,"omega(3)(1670)","omega_{3}(1670)");
  AddHadron(317,1.776,0.65,0.159,0,6,true,1,0,"K(3)*(1780)","K_{3}*(1780)");
  AddHadron(327,1.776,0.65,0.159,3,6,true,1,0,"K(3)*(1780)+","K_{3}*^{+}(1780)");
  AddHadron(kf_phi_3_1850,1.854,0.65,0.087,0,6,false,1,0,"phi(3)(1850)","phi_{3}(1850)");
  AddHadron(557,10.1599,0.65,0.0,0,6,true,1,0,"Upsilon(3)(1D)","Upsilon_{3}(1D)");
  // Tensors 4     ############################################################
  // heavy ones missing.
  AddHadron(kf_a_4_2040,2.014,0.65,0.361,0,8,false,1,0,"a(4)(2040)","a_{4}(2040)");
  AddHadron(219,2.014,0.65,0.361,3,8,true,1,0,"a(4)(2040)+","a_{4}^{+}(2040)");
  AddHadron(kf_f_4_2050,2.044,0.65,0.208,0,8,false,1,0,"f(4)(2050)","f_{4}(2050)");
  AddHadron(319,2.045,0.65,0.198,0,8,true,1,0,"K(4)*(2045)","K_{4}*(2045)");
  AddHadron(329,2.045,0.65,0.198,3,8,true,1,0,"K(4)*(2045)+","K_{4}*^{+}(2045)");
  // Tensors 2       ##################################################################################
  // heavy ones missing.
  AddHadron(kf_pi_2_1670,1.67,0.65,0.258,0,4,false,1,0,"pi(2)(1670)","pi_{2}(1670)");
  AddHadron(10215,1.67,0.65,0.258,3,4,true,1,0,"pi(2)(1670)+","pi_{2}^{+}(1670)");
  AddHadron(kf_eta_2_1645,1.617,0.65,0.181,0,4,false,1,0,"eta(2)(1645)","eta_{2}(1645)");
  AddHadron(10315,1.773,0.65,0.186,0,4,true,1,0,"K(2)(1770)","K_{2}(1770)");
  AddHadron(10325,1.773,0.65,0.186,3,4,true,1,0,"K(2)(1770)+","K_{2}^{+}(1770)");
  AddHadron(kf_eta_2_1870,1.842,0.65,0.225,0,4,false,1,0,"eta(2)(1870)","eta_{2}(1870)");
  AddHadron(10555,10.157,0.65,0.0,0,4,true,1,0,"eta(b2)(1D)","eta_{b2}(1D)");
  // Tensors 2       ##################################################################################
  // lots missing
  AddHadron(20315,1.816,0.65,0.276,0,4,true,1,0,"K(2)(1820)","K_{2}(1820)");
  AddHadron(20325,1.816,0.65,0.276,3,4,true,1,0,"K(2)(1820)+","K_{2}^{+}(1820)");
  AddHadron(20555,10.1562,0.65,0.0,0,4,true,1,0,"Upsilon(2)(1D)","Upsilon_{2}(1D)");
  AddHadron(30411,2.58,0.65,0.0,3,0,true,1,0,"D(2S)+","D(2S)^{+}");
  AddHadron(30421,2.58,0.65,0.0,0,0,true,1,0,"D(2S)","D(2S)");
  // Vectors 2       ##################################################################################
  // some heavy ones missing - we may have to include this because of the
  // psi(3770)
  AddHadron(kf_rho_1700,1.7,0.65,0.24,0,2,false,1,0,"rho(1700)","rho(1700)");
  AddHadron(30213,1.7,0.65,0.24,3,2,true,1,0,"rho(1700)+","rho^{+}(1700)");
  AddHadron(kf_omega_1600,1.670,0.65,0.31,0,2,false,1,0,"omega(1650)","omega(1650)");
  AddHadron(30313,1.717,0.65,0.32,0,2,true,1,0,"K*(1680)","K*(1680)");
  AddHadron(30323,1.717,0.65,0.32,3,2,true,1,0,"K*(1680)+","K*^{+}(1680)");
  AddHadron(kf_f1_1900_fict,1.900,0.65,0.32,0,2,false,1,0,"f1(1900)_fict","f1(1900)_fict");
  AddHadron(30413,2.64,0.65,0.0,3,2,true,1,0,"D*(2S)+","D*^{+}(2S)");
  AddHadron(30423,2.64,0.65,0.0,0,2,true,1,0,"D*(2S)","D*(2S)");
  AddHadron(kf_psi_3770,3.7699,0.65,0.0236,0,2,false,1,0,"psi(3770)","psi(3770)");
  AddHadron(30553,10.161,0.65,0.0,0,2,true,1,0,"Upsilon(1)(1D)","Upsilon_{1}(1D)");
  // Pseudoscalars   ##################################################################################
  // heavy ones missing.
  AddHadron(kf_pi_1300,1.3,0.65,0.400,0,0,false,1,0,"pi(1300)","pi(1300)");
  AddHadron(100211,1.3,0.65,0.400,3,0,true,1,0,"pi(1300)+","pi^{+}(1300)");
  AddHadron(kf_eta_1295,1.297,0.65,0.053,0,0,false,1,0,"eta(1295)","eta(1295)");
  AddHadron(100311,1.46,0.65,0.26,0,0,true,1,0,"K(1460)","K(1460)");
  AddHadron(100321,1.46,0.65,0.26,3,0,true,1,0,"K(1460)+","K^{+}(1460)");
  AddHadron(kf_eta_1475,1.476,0.65,0.08,0,0,false,1,0,"eta(1475)","eta(1475)");
  AddHadron(100441,3.638,0.65,0.014,0,0,true,1,0,"eta(c)(2S)","eta_{c}(2S)");
  AddHadron(100551,9.997,0.65,0.0,0,0,true,1,0,"eta(b)(2S)","eta_{b}(2S)");
  // Vectors         ##################################################################################
  // heavy ones missing.
  // the rho's may be important for tau/D/B decays
  AddHadron(kf_rho_1450,1.465,0.65,0.31,0,2,false,1,0,"rho(1450)","rho(1450)");
  AddHadron(100213,1.465,0.65,0.31,3,2,true,1,0,"rho(1450)+","rho^{+}(1450)");
  AddHadron(kf_omega_1420,1.419,0.65,0.17,0,2,false,1,0,"omega(1420)","omega(1420)");
  AddHadron(100313,1.414,0.65,0.232,0,2,true,1,0,"K*(1410)","K*(1410)");
  AddHadron(100323,1.414,0.65,0.232,3,2,true,1,0,"K*(1410)+","K*^{+}(1410)");
  AddHadron(kf_phi_1680,1.68,0.65,0.15,0,2,false,1,0,"phi(1680)","phi(1680)");
  AddHadron(kf_psi_2S,3.686,0.65,0.000277,0,2,false,1,0,"psi(2S)","psi(2S)");
  AddHadron(kf_Upsilon_2S,10.0233,0.65,4.4e-05,0,2,false,1,0,"Upsilon(2S)","Upsilon(2S)");
  // Tensors 2       ##################################################################################
  // heavy ones missing.
  AddHadron(kf_f_2_2010,2.01,0.65,0.2,0,4,false,1,0,"f(2)(2010)","f_{2}(2010)");
  AddHadron(100445,3.929,0.65,0.029,0,4,true,1,0,"chi(c2)(2P)","chi_{c2}(2P)");
  // More states without full multiplets ########################################################
  // light ones missing - we know A LOT about heavy-heavy states ....
  AddHadron(kf_chi_b2_2P,10.2685,0.65,0.001,0,4,false,1,0,"chi(b2)(2P)","chi_{b2}(2P)");
  AddHadron(100557,10.4443,0.65,0.0,0,6,true,1,0,"Upsilon(3)(2D)","Upsilon_{3}(2D)");
  AddHadron(kf_chi_b0_2P,10.2321,0.65,0.001,0,0,false,1,0,"chi(b0)(2P)","chi_{b0}(2P)");
  AddHadron(110553,10.255,0.65,0.0,0,2,true,1,0,"h(b)(2P)","h_{b}(2P)");
  AddHadron(110555,10.441,0.65,0.0,0,4,true,1,0,"eta(b2)(2D)","eta_{b2}(2D)");
  AddHadron(kf_chi_b1_2P,10.2552,0.65,0.001,0,2,false,1,0,"chi(b1)(2P)","chi_{b1}(2P)");
  AddHadron(120555,10.4406,0.65,0.0,0,4,true,1,0,"Upsilon(2)(2D)","Upsilon_{2}(2D)");
  AddHadron(130553,10.4349,0.65,0.0,0,2,true,1,0,"Upsilon(1)(2D)","Upsilon_{1}(2D)");
  AddHadron(200551,10.335,0.65,0.0,0,0,true,1,0,"eta(b)(3S)","eta_{b}(3S)");
  AddHadron(kf_Upsilon_3S,10.3553,0.65,2.63e-05,0,2,false,1,0,"Upsilon(3S)","Upsilon(3S)");
  AddHadron(200555,10.5264,0.65,0.0,0,4,true,1,0,"chi(b2)(3P)","chi_{b2}(3P)");
  AddHadron(210551,10.5007,0.65,0.0,0,0,true,1,0,"chi(b0)(3P)","chi_{b0}(3P)");
  AddHadron(210553,10.516,0.65,0.0,0,2,true,1,0,"h(b)(3P)","h_{b}(3P)");
  AddHadron(220553,10.516,0.65,0.0,0,2,true,1,0,"chi(b1)(3P)","chi_{b1}(3P)");
  AddHadron(kf_Upsilon_4S,10.58,0.65,0.01,0,2,false,1,0,"Upsilon(4S)","Upsilon(4S)");
  // a0(980) and friends #####################################################
  // These are the "funny" state with some iso number.  We will have to
  // figure out what to do with them.
  // #########################################################################
  AddHadron(kf_a_0_980,0.996,0.65,0.075,0,0,false,1,0,"a(0)(980)","a_{0}(980)");
  AddHadron(kf_f_0_600,0.600,0.65,0.600,0,0,false,1,0,"f(0)(600)","f_{0}(600)");
  AddHadron(kf_f_0_980,0.98,0.65,0.070,0,0,false,1,0,"f(0)(980)","f_{0}(980)");
  AddHadron(9000211,0.996,0.65,0.075,3,0,true,1,0,"a(0)(980)+","a_{0}^{+}(980)");
  AddHadron(9000311,0.841,0.65,0.618,0,0,true,1,0,"K(0)*(800)","K_{0}*(800)");
  AddHadron(9000321,0.841,0.65,0.618,3,0,true,1,0,"K(0)*(800)+","K_{0}*^{+}(800)");
  AddHadron(9000113,1.376,0.65,0.3,0,2,true,1,0,"pi(1)(1400)","pi_{1}(1400)");
  AddHadron(9000213,1.376,0.65,0.3,3,2,true,1,0,"pi(1)(1400)+","pi_{1}^{+}(1400)");
  AddHadron(9000223,1.518,0.65,0.073,0,2,true,1,0,"f(1)(1510)","f_{1}(1510)");
  AddHadron(9000313,1.65,0.65,0.15,0,2,true,1,0,"K(1)(1650)","K_{1}(1650)");
  AddHadron(9000323,1.65,0.65,0.15,3,2,true,1,0,"K(1)(1650)+","K_{1}^{+}(1650)");
  AddHadron(kf_psi_4040,4.04,0.65,0.052,0,2,false,1,0,"psi(4040)","psi(4040)");
  AddHadron(kf_Upsilon_10860,10.865,0.65,0.11,0,2,false,1,0,"Upsilon(10860)","Upsilon(10860)");
  AddHadron(9000115,1.732,0.65,0.194,0,4,true,1,0,"a(2)(1700)","a_{2}(1700)");
  AddHadron(9000215,1.732,0.65,0.194,3,4,true,1,0,"a(2)(1700)+","a_{2}^{+}(1700)");
  AddHadron(9000225,1.43,0.65,0.02,0,4,true,1,0,"f(2)(1430)","f_{2}(1430)");
  AddHadron(9000315,1.58,0.65,0.11,0,4,true,1,0,"K(2)(1580)","K_{2}(1580)");
  AddHadron(9000325,1.58,0.65,0.11,3,4,true,1,0,"K(2)(1580)+","K_{2}^{+}(1580)");
  AddHadron(9000117,1.982,0.65,0.188,0,6,true,1,0,"rho(3)(1990)","rho_{3}(1990)");
  AddHadron(9000217,1.982,0.65,0.188,3,6,true,1,0,"rho(3)(1990)+","rho_{3}^{+}(1990)");
  AddHadron(9000229,2.2311,0.65,0.023,0,8,true,1,0,"f(J)(2220)","f(J)(2220)");
  AddHadron(9000319,2.045,0.65,0.198,0,8,true,1,0,"K(4)(2500)","K_{4}(2500)");
  AddHadron(9000329,2.045,0.65,0.198,3,8,true,1,0,"K(4)(2500)+","K_{4}^{+}(2500)");
  AddHadron(9010111,1.812,0.65,0.207,0,0,true,1,0,"pi(1800)","pi(1800)");
  AddHadron(9010211,1.812,0.65,0.207,3,0,true,1,0,"pi(1800)+","pi^{+}(1800)");
  AddHadron(9010311,1.83,0.65,0.25,0,0,true,1,0,"K(1830)","K(1830)");
  AddHadron(9010321,1.83,0.65,0.25,3,0,true,1,0,"K(1830)+","K^{+}(1830)");
  AddHadron(9010113,1.653,0.65,0.225,0,2,true,1,0,"pi(1)(1600)","pi_{1}(1600)");
  AddHadron(9010213,1.653,0.65,0.225,3,2,true,1,0,"pi(1)(1600)+","pi_{1}^{+}(1600)");
  AddHadron(9010223,1.594,0.65,0.384,0,2,true,1,0,"h(1)(1595)","h_{1}(1595)");
  AddHadron(kf_Upsilon_11020,11.019,0.65,0.079,0,2,false,1,0,"Upsilon(11020)","Upsilon(11020)");
  AddHadron(kf_psi_4160,4.159,0.65,0.078,0,2,false,1,0,"psi(4160)","psi(4160)");
  AddHadron(9010115,2.09,0.65,0.625,0,4,true,1,0,"pi(2)(2100)","pi_{2}(2100)");
  AddHadron(9010215,2.09,0.65,0.625,3,4,true,1,0,"pi(2)(2100)+","pi_{2}^{+}(2100)");
  AddHadron(9010225,1.546,0.65,0.126,0,4,true,1,0,"f(2)(1565)","f_{2}(1565)");
  AddHadron(9010315,1.973,0.65,0.373,0,4,true,1,0,"K(2)*(1980)","K_{2}*(1980)");
  AddHadron(9010325,1.973,0.65,0.373,3,4,true,1,0,"K(2)*(1980)+","K_{2}*^{+}(1980)");
  AddHadron(9010117,2.25,0.65,0.2,0,6,true,1,0,"rho(3)(2250)","rho_{3}(2250)");
  AddHadron(9010217,2.25,0.65,0.2,3,6,true,1,0,"rho(3)(2250)+","rho_{3}^{+}(2250)");
  AddHadron(9010317,2.324,0.65,0.15,0,6,true,1,0,"K(3)(2320)","K_{3}(2320)");
  AddHadron(9010327,2.324,0.65,0.15,3,6,true,1,0,"K(3)(2320)+","K_{3}^{+}(2320)");
  AddHadron(9010229,2.332,0.65,0.26,0,8,true,1,0,"f(4)(2300)","f_{4}(2300)");
  AddHadron(kf_f_0_1500,1.4098,0.65,0.0511,0,0,false,1,0,"eta(1405)","eta(1405)");
  AddHadron(9020311,1.945,0.65,0.201,0,0,true,1,0,"K(0)*(1950)","K_{0}*(1950)");
  AddHadron(9020321,1.945,0.65,0.201,3,0,true,1,0,"K(0)*(1950)+","K_{0}*^{+}(1950)");
  AddHadron(9020113,1.647,0.65,0.254,0,2,true,1,0,"a(1)(1640)","a_{1}(1640)");
  AddHadron(9020213,1.647,0.65,0.254,3,2,true,1,0,"a(1)(1600)+","a_{1}^{+}(1600)");
  AddHadron(kf_psi_4415,4.415,0.65,0.043,0,2,false,1,0,"psi(4415)","psi(4415)");
  AddHadron(9020225,1.638,0.65,0.099,0,4,true,1,0,"f(2)(1640)","f_{2}(1640)");
  AddHadron(9020315,2.247,0.65,0.18,0,4,true,1,0,"K(2)(2250)","K_{2}(2250)");
  AddHadron(9020325,2.247,0.65,0.18,3,4,true,1,0,"K(2)(2250)+","K_{2}^{+}(2250)");
  AddHadron(kf_f_J_1710,1.5,0.65,0.112,0,0,false,1,0,"f(0)(1500)","f_{0}(1500)");
  AddHadron(9030113,1.9,0.65,0.16,0,2,true,1,0,"rho(1900)","rho(1900)");
  AddHadron(9030213,1.9,0.65,0.16,3,2,true,1,0,"rho(1900)+","rho^{+}(1900)");
  AddHadron(9030225,1.815,0.65,0.197,0,4,true,1,0,"f(2)(1810)","f_{2}(1810)");
  AddHadron(9040221,1.76,0.65,0.06,0,0,true,1,0,"eta(1760)","eta(1760)");
  AddHadron(9040113,2.149,0.65,0.363,0,2,true,1,0,"rho(2150)","rho(2150)");
  AddHadron(9040213,2.149,0.65,0.363,3,2,true,1,0,"rho(2150)+","rho^{+}(2150)");
  AddHadron(9040225,1.915,0.65,0.163,0,4,true,1,0,"f(2)(1910)","f_{2}(1910)");
  AddHadron(9050221,1.992,0.65,0.442,0,0,true,1,0,"f(0)(2020)","f_{0}(2020)");
  AddHadron(kf_f_2_2300,1.944,0.65,0.472,0,4,false,1,0,"f(2)(1950)","f_{2}(1950)");
  AddHadron(9060221,2.103,0.65,0.206,0,0,true,1,0,"f(0)(2100)","f_{0}(2100)");
  AddHadron(kf_f_2_2340,2.011,0.65,0.202,0,4,false,1,0,"f(2)(2011)","f_{2}(2011)");
  AddHadron(9070221,2.189,0.65,0.238,0,0,true,1,0,"f(0)(2200)","f_{0}(2200)");
  AddHadron(9070225,2.15,0.65,0.2,0,4,true,1,0,"f(2)(2150)","f_{2}(2150)");
  AddHadron(9080221,2.220,0.65,0.15,0,0,true,1,0,"eta(2225)","eta(2225)");
  AddHadron(9080225,2.297,0.65,0.15,0,4,true,1,0,"f(2)(2300)","f_{2}(2300)");
  AddHadron(9090225,2.34,0.65,0.32,0,4,true,1,0,"f(2)(2340)","f_{2}(2340)");

  // ##########################################################################
  // ##########################################################################
  // Former members of the Sherpa team - made immortal as particles here ######
  // ##########################################################################
  // ##########################################################################
  AddHadron(5505,1000000.0,1.e6,1000,0,0,true,0,0,"ralf","ralf");
  AddHadron(5506,1000000.0,1.e6,1000,0,0,true,0,0,"ande","ande");
  AddHadron(5507,1000000.0,1.e6,1000,0,0,true,0,0,"thomas","thomas");
  AddHadron(5508,1000000.0,1.e6,1000,0,0,true,0,0,"tanju","tanju");
  AddHadron(5509,1000000.0,1.e6,1000,0,0,true,0,0,"jennifer","jennifer");
  AddHadron(6505,1000000.0,1.e6,1000,0,0,true,0,0,"hendrik","hendrik");
  AddHadron(6506,1000000.0,1.e6,1000,0,0,true,0,0,"jan","jan");
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################
}
