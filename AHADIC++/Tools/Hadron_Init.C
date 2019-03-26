#include "AHADIC++/Tools/Hadron_Init.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

void Hadron_Init::Init() {
  msg_Info()<<METHOD<<"(): Initializing kf table for hadrons.\n";

  // ##########################################################################
  // ##########################################################################
  // Particle being used in AHADIC listed below
  // ##########################################################################
  //
  // ##########################################################################
  // Diquarks - not exactly particles but we need them ########################
  // ##########################################################################
  // ##########################################################################
  if(s_kftable.find(1103)==s_kftable.end()) { // if not initialized in Beam
    InitHadron(1103,0.77133,0,-2,-3,2,0,0,1,1,"dd_1",  "dd_1b","dd_1b","dd_1b");
    InitHadron(2101,0.57933,0,1,-3,0,0,0,1,1,"ud_0",   "ud_0b","ud_0b","ud_0b");
    InitHadron(2103,0.77133,0,1,-3,2,0,0,1,1,"ud_1",   "ud_1b","ud_1b","ud_1b");
    InitHadron(2203,0.77133,0,4,-3,2,0,0,1,1,"uu_1",   "uu_1b","uu_1b","uu_1b");
    InitHadron(3101,0.80473,0,-2,-3,0,0,0,1,1,"sd_0",  "sd_0b","sd_0b","sd_0b");
    InitHadron(3103,0.92953,0,-2,-3,2,0,0,1,1,"sd_1",  "sd_1b","sd_1b","sd_1b");
    InitHadron(3201,0.80473,0,1,-3,0,0,0,1,1,"su_0",   "su_0b","su_0b","su_0b");
    InitHadron(3203,0.92953,0,1,-3,2,0,0,1,1,"su_1",   "su_1b","su_1b","su_1b");
    InitHadron(3303,1.09361,0,-2,-3,2,0,0,1,1,"ss_1",  "ss_1b","ss_1b","ss_1b");
    InitHadron(4101,1.96908,0,1,-3,0,0,0,1,1,"cd_0",   "cd_0b","cd_0b","cd_0b");
    InitHadron(4103,2.00808,0,1,-3,2,0,0,1,1,"cd_1",   "cd_1b","cd_1b","cd_1b");
    InitHadron(4201,1.96908,0,4,-3,0,0,0,1,1,"cu_0",   "cu_0b","cu_0b","cu_0b");
    InitHadron(4203,2.00808,0,4,-3,2,0,0,1,1,"cu_1",   "cu_1b","cu_1b","cu_1b");
    InitHadron(4301,2.15432,0,1,-3,0,0,0,1,1,"cs_0",   "cs_0b","cs_0b","cs_0b");
    InitHadron(4303,2.17967,0,1,-3,2,0,0,1,1,"cs_1",   "cs_1b","cs_1b","cs_1b");
    InitHadron(4403,3.27531,0,4,-3,2,0,0,1,1,"cc_1",   "cc_1b","cc_1b","cc_1b");
    InitHadron(5101,5.38897,0,-2,-3,0,0,0,1,1,"bd_0",  "bd_0b","bd_0b","bd_0b");
    InitHadron(5103,5.40145,0,-2,-3,2,0,0,1,1,"bd_1",  "bd_1b","bd_1b","bd_1b");
    InitHadron(5201,5.38897,0,1,-3,0,0,0,1,1,"bu_0",   "bu_0b","bu_0b","bu_0b");
    InitHadron(5203,5.40145,0,1,-3,2,0,0,1,1,"bu_1",   "bu_1b","bu_1b","bu_1b");
    InitHadron(5301,5.56725,0,-2,-3,0,0,0,1,1,"bs_0",  "bs_0b","bs_0b","bs_0b");
    InitHadron(5303,5.57536,0,-2,-3,2,0,0,1,1,"bs_1",  "bs_1b","bs_1b","bs_1b");
    InitHadron(5401,6.67143,0,1,-3,0,0,0,1,1,"bc_0",   "bc_0b","bc_0b","bc_0b");
    InitHadron(5403,6.67397,0,1,-3,2,0,0,1,1,"bc_1",   "bc_1b","bc_1b","bc_1b");
    InitHadron(5503,10.07354,0,-2,-3,2,0,0,1,1,"bb_1", "bb_1b","bb_1b","bb_1b");
  }
  // ##########################################################################
  // ##########################################################################
  // Other objects that may show up in an event record ########################
  InitHadron(kf_cluster,0.,0.,0,0,0,0,1,1,0,"cluster","cluster","cluster","cluster");
  InitHadron(kf_string,0.,0.,0,0,0,0,1,1,0,"string","string","string","string");
  // ##################################################################################################
  // MESON MULTIPLETS
  // ##################################################################################################
  // Pseudoscalars   ##################################################################################
  InitHadron(kf_pi,true,0.134976,7.8486e-09,0,0,false,1,0,"pi","pi");
  InitHadron(211,true,0.13957,2.5242e-17,3,0,true,1,1,"pi+","pi^{+}");
  InitHadron(kf_eta,false,0.5473,1.18e-06,0,0,false,1,0,"eta","eta");
  InitHadron(311,false,0.49767,1.e-16,0,0,true,1,0,"K","K");
  InitHadron(kf_K_L,false,0.49767,1.273e-17,0,0,false,1,1,"K(L)","K_{L}");
  InitHadron(kf_K_S,false,0.49767,7.373e-15,0,0,false,1,0,"K(S)","K_{S}");
  InitHadron(321,false,0.493677,5.314e-17,3,0,true,1,1,"K+","K^{+}");
  InitHadron(kf_eta_prime_958,false,0.95778,0.000203,0,0,false,1,0,"eta'(958)","eta'(958)");
  InitHadron(411,false,1.8693,6.23e-13,3,0,true,1,0,"D+","D^{+}");
  InitHadron(421,false,1.8646,1.586e-12,0,0,true,1,0,"D","D");
  InitHadron(431,false,1.9685,1.41e-12,3,0,true,1,0,"D(s)+","D(s)^{+}");
  InitHadron(kf_eta_c_1S,false,2.9798,0.0132,0,0,false,1,0,"eta(c)(1S)","eta_{c}(1S)");
  InitHadron(511,false,5.2792,4.22e-13,0,0,true,1,0,"B","B");
  InitHadron(521,false,5.2789,3.99e-13,3,0,true,1,0,"B+","B^{+}");
  InitHadron(531,false,5.3693,4.27e-13,0,0,true,1,0,"B(s)","B_{s}");
  InitHadron(541,false,6.4,1.43e-12,3,0,true,1,0,"B(c)+","B_{c}^{+}");
  InitHadron(kf_eta_b,false,9.4,0.050,0,0,false,1,0,"eta(b)(1S)","eta_{b}(1S)");
  // Vectors         ##################################################################################
  InitHadron(kf_rho_770,false,0.77,0.1507,0,2,false,1,0,"rho(770)","rho(770)");
  InitHadron(213,false,0.77,0.1507,3,2,true,1,0,"rho(770)+","rho^{+}(770)");
  InitHadron(kf_omega_782,false,0.78194,0.00841,0,2,false,1,0,"omega(782)","omega(782)");
  InitHadron(313,false,0.8961,0.0505,0,2,true,1,0,"K*(892)","K*(892)");
  InitHadron(323,false,0.89166,0.0508,3,2,true,1,0,"K*(892)+","K*^{+}(892)");
  InitHadron(kf_phi_1020,false,1.01941,0.00443,0,2,false,1,0,"phi(1020)","phi(1020)");
  InitHadron(413,false,2.01,0.001,3,2,true,1,0,"D*(2010)+","D*^{+}(2010)");
  InitHadron(423,false,2.0067,0.001,0,2,true,1,0,"D*(2007)","D*(2007)");
  InitHadron(433,false,2.1124,0.001,3,2,true,1,0,"D(s)*+","D_{s}*^{+}");
  InitHadron(kf_J_psi_1S,false,3.09688,8.7e-05,0,2,false,1,0,"J/psi(1S)","J/psi(1S)");
  InitHadron(513,false,5.3249,0.0001,0,2,true,1,0,"B*","B*");
  InitHadron(523,false,5.3249,0.0001,3,2,true,1,0,"B*+","B*^{+}");
  InitHadron(533,false,5.41630,0.0001,0,2,true,1,0,"B(s)*","B_{s}*");
  InitHadron(543,false,6.602,0.0001,3,2,true,1,0,"B(c)*+","B_{c}*^{+}");
  InitHadron(kf_Upsilon_1S,false,9.46037,5.25e-05,0,2,false,1,0,"Upsilon(1S)","Upsilon(1S)");
  // Tensors 2       ##################################################################################
  InitHadron(kf_a_2_1320,false,1.3181,0.107,0,4,false,1,0,"a(2)(1320)","a_{2}(1320)");
  InitHadron(215,false,1.3181,0.107,3,4,true,1,0,"a(2)(1320)+","a_{2}^{+}(1320)");
  InitHadron(kf_f_2_1270,false,1.275,0.1855,0,4,false,1,0,"f(2)(1270)","f_{2}(1270)");
  InitHadron(315,false,1.4324,0.109,0,4,true,1,0,"K(2)*(1430)","K_{2}*(1430)");
  InitHadron(325,false,1.4256,0.0985,3,4,true,1,0,"K(2)*(1430)+","K_{2}*^{+}(1430)");
  InitHadron(kf_f_2_prime_1525,false,1.525,0.076,0,4,false,1,0,"f(2)'(1525)","f_{2}'(1525)");
  InitHadron(415,false,2.459,0.025,3,4,true,1,0,"D(2)*(2460)+","D_{2}*^{+}(2460)");
  InitHadron(425,false,2.4589,0.023,0,4,true,1,0,"D(2)*(2460)","D_{2}*(2460)");
  InitHadron(435,false,2.5735,0.015,3,4,true,1,0,"D(s2)*(2573)+","D_{s2}*^{+}(2573)");
  InitHadron(kf_chi_c2_1P,false,3.55617,0.002,0,4,false,1,0,"chi(c2)(1P)","chi_{c2}(1P)");
  InitHadron(515,false,5.83,0.02,0,4,true,1,0,"B(2)*","B_{2}*");
  InitHadron(525,false,5.83,0.02,3,4,true,1,0,"B(2)*+","B_{2}*^{+}");
  InitHadron(535,false,5.8397,0.02,0,4,true,1,0,"B(s2)*","B_{s2}*");
  InitHadron(545,false,7.35,0.02,3,4,true,1,0,"B(c2)*+","B_{c2}*^{+}");
  InitHadron(kf_chi_b2_1P,false,9.9132,0.001,0,4,false,1,0,"chi(b2)(1P)","chi_{b2}(1P)");
  // Scalars         ##################################################################################
  InitHadron(kf_a_0_1450,false,1.474,0.265,0,0,false,1,0,"a(0)(1450)","a_{0}(1450)");
  InitHadron(kf_a_0_1450_plus,false,1.474,0.265,3,0,true,1,0,"a(0)(1450)+","a_{0}^{+}(1450)");
  InitHadron(kf_f_0_1370,false,1.4,0.5,0,0,false,1,0,"f(0)(1370)","f_{0}(1370)");
  InitHadron(10311,false,1.429,0.287,0,0,true,1,0,"K(0)*(1430)","K_{0}*(1430)");
  InitHadron(10321,false,1.429,0.287,3,0,true,1,0,"K(0)*(1430)+","K_{0}*^{+}(1430)");
  InitHadron(kf_f_0_1710,false,1.7,0.5,0,0,false,1,0,"f(0)(1710)","f_{0}(1710)");
  InitHadron(10411,false,2.351,0.23,3,0,true,1,0,"D(0)*(2400)+","D_{0}*^{+}(2400)");
  InitHadron(10421,false,2.318,0.267,0,0,true,1,0,"D(0)*(2400)","D_{0}*(2400)");
  InitHadron(10431,false,2.3177,0.03,3,0,true,1,0,"D(s0)*(2317)+","D_{s0}*^{+}(2317)");
  InitHadron(kf_chi_c0_1P,false,3.4173,0.014,0,0,false,1,0,"chi(c0)(1P)","chi_{c0}(1P)");
  InitHadron(10511,false,5.68,0.05,0,0,true,1,0,"B(0)*","B_{0}*");
  InitHadron(10521,false,5.68,0.05,3,0,true,1,0,"B(0)*+","B_{0}*^{+}");
  InitHadron(10531,false,5.92,0.05,0,0,true,1,0,"B(s0)*","B_{s0}*");
  InitHadron(10541,false,7.25,0.05,3,0,true,1,0,"B(c0)*+","B_{c0}*^{+}");
  InitHadron(kf_chi_b0_1P,false,9.8598,0.050,0,0,false,1,0,"chi(b0)(1P)","chi_{b0}(1P)");
  // Axial vectors   ##################################################################################
  InitHadron(kf_b_1_1235,false,1.2295,0.142,0,2,false,1,0,"b(1)(1235)","b_{1}(1235)");
  InitHadron(10213,false,1.2295,0.142,3,2,true,1,0,"b(1)(1235)+","b_{1}^{+}(1235)");
  InitHadron(kf_h_1_1170,false,1.17,0.36,0,2,false,1,0,"h(1)(1170)","h_{1}(1170)");
  InitHadron(10313,false,1.272,0.09,0,2,true,1,0,"K(1)(1270)","K_{1}(1270)");
  InitHadron(10323,false,1.272,0.09,3,2,true,1,0,"K(1)(1270)+","K_{1}^{+}(1270)");
  InitHadron(kf_h_1_1380,false,1.386,0.091,0,2,false,1,0,"h(1)(1380)","h_{1}(1380)");
  InitHadron(10413,false,2.4232,0.025,3,2,true,1,0,"D(1)(2420)+","D_{1}^{+}(2420)");
  InitHadron(10423,false,2.4208,0.0317,0,2,true,1,0,"D(1)(2420)","D_{1}(2420)");
  InitHadron(10433,false,2.5351,0.00092,3,2,true,1,0,"D(s1)(2536)+","D_{s1}^{+}(2536)");
  InitHadron(kf_h_c1,false,3.46,0.01,0,2,false,1,0,"h(c)(1P)","h_{c}(1P)");
  InitHadron(10513,false,5.73,0.05,0,2,true,1,0,"B(1)(L)","B_{1}(L)");
  InitHadron(10523,false,5.73,0.05,3,2,true,1,0,"B(1)(L)+","B_{1}^{+}(L)");
  InitHadron(10533,false,5.97,0.05,0,2,true,1,0,"B(s1)(L)","B_{s1}(L)");
  InitHadron(10543,false,7.3,0.05,3,2,true,1,0,"B(c1)(L)+","B_{c1}^{+}(L)");
  InitHadron(kf_h_b1,false,9.875,0.01,0,2,false,1,0,"h(b)(1P)","h_{b}(1P)");
  // Vectors         ##################################################################################
  InitHadron(kf_a_1_1260,false,1.23,0.400,0,2,false,1,0,"a(1)(1260)","a_{1}(1260)");
  InitHadron(20213,false,1.23,0.400,3,2,true,1,0,"a(1)(1260)+","a_{1}^{+}(1260)");
  InitHadron(kf_f_1_1285,false,1.2819,0.024,0,2,false,1,0,"f(1)(1285)","f_{1}(1285)");
  InitHadron(20313,false,1.402,0.174,0,2,true,1,0,"K(1)(1400)","K_{1}(1400)");
  InitHadron(20323,false,1.402,0.174,3,2,true,1,0,"K(1)(1400)+","K_{1}^{+}(1400)");
  InitHadron(kf_f_1_1420,false,1.4262,0.055,0,2,false,1,0,"f(1)(1420)","f_{1}(1420)");
  InitHadron(20413,false,2.427,0.384,3,2,true,1,0,"D(1)(H)+","D_{1}^{+}(H)");
  InitHadron(20423,false,2.427,0.384,0,2,true,1,0,"D(1)(2430)","D_{1}(2430)");
  InitHadron(20433,false,2.4595,0.003,3,2,true,1,0,"D(s1)(2460)+","D_{s1}^{+}(2460)");
  InitHadron(kf_chi_c1_1P,false,3.51053,0.00088,0,2,false,1,0,"chi(c1)(1P)","chi_{c1}(1P)");
  InitHadron(20513,false,5.78,0.05,0,2,true,1,0,"B(1)(H)","B_{1}(H)");
  InitHadron(20523,false,5.78,0.05,3,2,true,1,0,"B(1)(H)+","B_{1}^{+}(H)");
  InitHadron(20533,false,6.02,0.05,0,2,true,1,0,"B(s1)(H)","B_{s1}(H)");
  InitHadron(20543,false,7.3,0.05,3,2,true,1,0,"B(c1)(H)+","B_{c1}^{+}(H)");
  InitHadron(kf_chi_b1_1P,false,9.8919,0.001,0,2,false,1,0,"chi(b1)(1P)","chi_{b1}(1P)");
  // ##################################################################################################
  // BARYON MULTIPLETS
  // ##################################################################################################
  // Nucleons (octet) #################################################################################
  InitHadron(2112,true,0.939566,7.424e-28,0,1,true,1,1,"n","n");
  InitHadron(2212,true,0.938272,0,3,1,true,1,1,"P+","P^{+}");
  InitHadron(3112,false,1.19745,4.45e-15,-3,1,true,1,0,"Sigma-","\\Sigma^{-}");
  InitHadron(3212,false,1.19264,8.9e-06,0,1,true,1,0,"Sigma","\\Sigma");
  InitHadron(3222,false,1.18937,8.24e-15,3,1,true,1,0,"Sigma+","\\Sigma^{+}");
  InitHadron(3122,false,1.11568,2.501e-15,0,1,true,1,0,"Lambda","\\Lambda");
  InitHadron(3312,false,1.32132,4.02e-15,-3,1,true,1,0,"Xi-","\\Xi^{-}");
  InitHadron(3322,false,1.3149,2.27e-15,0,1,true,1,0,"Xi","\\Xi");
  InitHadron(4112,false,2.4522,0.0022,0,1,true,1,0,"Sigma(c)(2455)","\\Sigma_{c}(2455)");
  InitHadron(4212,false,2.4536,0.0023,3,1,true,1,0,"Sigma(c)(2455)+","\\Sigma_{c}^{+}(2455)");
  InitHadron(4222,false,2.4528,0.0023,6,1,true,1,0,"Sigma(c)(2455)++","\\Sigma_{c}^{++}(2455)");
  InitHadron(4122,false,2.2849,3.19e-12,3,1,true,1,0,"Lambda(c)+","\\Lambda_{c}^{+}");
  InitHadron(4132,false,2.4703,5.875e-12,0,1,true,1,0,"Xi(c)","\\Xi_{c}");
  InitHadron(4232,false,2.4656,1.489e-12,3,1,true,1,0,"Xi(c)+","\\Xi_{c}^{+}");
  InitHadron(4312,false,2.575,0.001,0,1,true,1,0,"Xi(c)'","\\Xi'_{c}");
  InitHadron(4322,false,2.578,0.001,3,1,true,1,0,"Xi(c)'+","\\Xi'_{c}^{+}");
  InitHadron(4332,false,2.704,1.02e-11,0,1,true,1,0,"Omega(c)","\\Omega_{c}");
  InitHadron(5112,false,5.810,0.001,-3,1,true,1,0,"Sigma(b)-","\\Sigma_{b}^{-}");
  InitHadron(5212,false,5.810,0.001,0,1,true,1,0,"Sigma(b)","\\Sigma_{b}");
  InitHadron(5222,false,5.810,0.001,3,1,true,1,0,"Sigma(b)+","\\Sigma_{b}^{+}");
  InitHadron(5122,false,5.624,5.31e-13,0,1,true,1,0,"Lambda(b)","\\Lambda_{b}");
  InitHadron(5132,false,5.790,1.e-13,-3,1,true,1,0,"Xi(b)-","\\Xi_{b}^{-}");
  InitHadron(5232,false,5.790,1.e-13,0,1,true,1,0,"Xi(b)","\\Xi_{b}");
  InitHadron(5312,false,5.890,1.e-13,-3,1,true,1,0,"Xi(b)'-_fict","\\Xi'_{b}_fict^{-}");
  InitHadron(5322,false,5.890,1.e-13,0,1,true,1,0,"Xi(b)'_fict","\\Xi'_{b}_fict");
  InitHadron(5332,false,6.071,1.e-13,-3,1,true,1,0,"Omega(b)-","\\Omega_{b}^{-}");
  // Deltas (decuplet) ################################################################################
  InitHadron(1114,false,1.232,0.12,-3,3,true,1,0,"Delta(1232)-","\\Delta-(1232)");
  InitHadron(2114,false,1.232,0.12,0,3,true,1,0,"Delta(1232)","\\Delta(1232)");
  InitHadron(2214,false,1.232,0.12,3,3,true,1,0,"Delta(1232)+","\\Delta^{+}(1232)");
  InitHadron(2224,false,1.232,0.12,6,3,true,1,0,"Delta(1232)++","\\Delta^{++}(1232)");
  InitHadron(3114,false,1.3872,0.0394,-3,3,true,1,0,"Sigma(1385)-","\\Sigma-(1385)");
  InitHadron(3214,false,1.3837,0.036,0,3,true,1,0,"Sigma(1385)","\\Sigma(1385)");
  InitHadron(3224,false,1.3828,0.0358,3,3,true,1,0,"Sigma(1385)+","\\Sigma^{+}(1385)");
  InitHadron(3314,false,1.535,0.0099,-3,3,true,1,0,"Xi(1530)-","\\Xi-(1530)");
  InitHadron(3324,false,1.5318,0.0091,0,3,true,1,0,"Xi(1530)","\\Xi(1530)");
  InitHadron(3334,false,1.67245,8.01e-15,-3,3,true,1,0,"Omega-","\\Omega^{-}");
  InitHadron(4114,false,2.5175,0.0150,0,3,true,1,0,"Sigma(c)(2520)","\\Sigma_{c}(2520)");
  InitHadron(4214,false,2.5159,0.0150,3,3,true,1,0,"Sigma(c)(2520)+","\\Sigma_{c}^{+}(2520)");
  InitHadron(4224,false,2.5194,0.0150,6,3,true,1,0,"Sigma(c)(2520)++","\\Sigma_{c}^{++}(2520)");
  InitHadron(4314,false,2.645,0.003,0,3,true,1,0,"Xi(c)*","\\Xi*_{c}");
  InitHadron(4324,false,2.645,0.003,3,3,true,1,0,"Xi(c)*+","\\Xi*_{c}^{+}");
  InitHadron(4334,false,2.8,0.001,0,3,true,1,0,"Omega(c)*","\\Omega*_{c}");
  InitHadron(5114,false,5.829,0.01,-3,3,true,1,0,"Sigma(b)*-","\\Sigma*_{b}^{-}");
  InitHadron(5214,false,5.829,0.01,0,3,true,1,0,"Sigma(b)*","\\Sigma*_{b}");
  InitHadron(5224,false,5.829,0.01,3,3,true,1,0,"Sigma(b)*+","\\Sigma*_{b}^{+}");
  InitHadron(5314,false,5.930,0.001,-3,3,true,1,0,"Xi(b)*-_fict","\\Xi*_{b}_fict^{-}");
  InitHadron(5324,false,5.930,0.001,0,3,true,1,0,"Xi(b)*_fict","\\Xi*_{b}_fict");
  InitHadron(5334,false,6.090,0.0003,-3,3,true,1,0,"Omega(b)*-_fict","\\Omega*_{b}_fict^{-}");
  // Nucleons (octet) - L_N = 1 ##############################################
  // careful - we will have to add some heavy states here!
  InitHadron(102112,false,1.535,0.15,0,1,true,1,0,"N(1535)","N(1535)");
  InitHadron(102212,false,1.535,0.15,3,1,true,1,0,"N(1535)+","N^{+}(1535)");
  InitHadron(102132,false,1.407,0.05,0,1,true,1,0,"Lambda(1405)","\\Lambda(1405)");
  InitHadron(103122,false,1.67,0.06,0,1,true,1,0,"Lambda(1670)","\\Lambda(1670)");
  InitHadron(103112,false,1.62,0.09,-3,1,true,1,0,"Sigma(1620)-","\\Sigma-(1620)");
  InitHadron(103212,false,1.62,0.09,0,1,true,1,0,"Sigma(1620)","\\Sigma(1620)");
  InitHadron(103222,false,1.62,0.09,3,1,true,1,0,"Sigma(1620)+","\\Sigma^{+}(1620)");
  InitHadron(103312,false,1.75,0.09,-3,1,true,1,0,"Xi(1750)","\\Xi^{-}(1750)");
  InitHadron(103322,false,1.75,0.09,0,1,true,1,0,"Xi(1750)","\\Xi(1750)");
  InitHadron(102142,false,2.5954,0.0036,3,1,true,1,0,"Lambda(c)(2595)+","\\Lambda_{c}(2595)^{+}");
  // Nucleons (octet) - the Roper resonance - L_N = 2 #########################
  // careful - we will have to add some heavy states here!
  InitHadron(202112,true,1.44,0.35,0,1,true,1,0,"N(1440)","N(1440)");
  InitHadron(202212,true,1.44,0.35,3,1,true,1,0,"N(1440)+","N^{+}(1440)");
  InitHadron(203112,false,1.66,0.1,-3,1,true,1,0,"Sigma(1660)-","\\Sigma-(1660)");
  InitHadron(203212,false,1.66,0.1,0,1,true,1,0,"Sigma(1660)","\\Sigma(1660)");
  InitHadron(203222,false,1.66,0.1,3,1,true,1,0,"Sigma(1660)+","\\Sigma^{+}(1660)");
  InitHadron(203122,false,1.6,0.15,0,1,true,1,0,"Lambda(1600)","\\Lambda(1600)");
  InitHadron(203312,false,1.696,0.010,-3,1,true,1,0,"Xi(1690)-","\\Xi-(1690)");
  InitHadron(203322,false,1.696,0.010,0,1,true,1,0,"Xi(1690)","\\Xi(1690)");
  // Nucleons (octet) - L_N = 1_1 -- plus "singlet heavies" ###########################################
  InitHadron(102114,false,1.52,0.12,0,3,true,1,0,"N(1520)","N(1520)");
  InitHadron(102214,false,1.52,0.12,3,3,true,1,0,"N(1520)+","N^{+}(1520)");
  InitHadron(103114,false,1.67,0.06,-3,3,true,1,0,"Sigma(1670)-","\\Sigma-(1670)");
  InitHadron(103214,false,1.67,0.06,0,3,true,1,0,"Sigma(1670)","\\Sigma(1670)");
  InitHadron(103224,false,1.67,0.06,3,3,true,1,0,"Sigma(1670)+","\\Sigma^{+}(1670)");
  InitHadron(103124,false,1.69,0.06,0,3,true,1,0,"Lambda(1690)","\\Lambda(1690)");
  InitHadron(103314,false,1.823,0.024,-3,3,true,1,0,"Xi(1820)-","\\Xi-(1820)");
  InitHadron(103324,false,1.823,0.024,0,3,true,1,0,"Xi(1820)","\\Xi(1820)");
  InitHadron(102134,false,1.5195,0.0156,0,3,true,1,0,"Lambda(1520)","\\Lambda(1520)");
  InitHadron(102144,false,2.625,0.002,3,3,true,1,0,"Lambda(c)(2625)+","\\Lambda_{c}^{+}(2625)");
  InitHadron(104314,false,2.815,0.002,0,3,true,1,0,"Xi(c)(2815)","\\Xi_{c}(2815)");
  InitHadron(104324,false,2.815,0.002,3,3,true,1,0,"Xi(c)(2815)+","\\Xi_{c}^{+}(2815)");
  InitHadron(102154,false,5.91,0.002,0,3,true,1,0,"Lambda(b)(5910)","\\Lambda_{b}(5910)");
  // #########################################################################
  // Obsolete multiple heavy baryons #########################################
  // they will not be produced in our code (we have no heavy di-quarks
  // that we can produce yet)
  // #########################################################################
  InitHadron(4412,false,3.59798,0.,3,1,true,1,0,"Xi(cc)+","\\Xi(cc)^{+}");
  InitHadron(4414,false,3.65648,0.,3,3,true,1,0,"Xi(cc)*+","\\Xi(cc)*^{+}");
  InitHadron(4422,false,3.59798,0.,6,1,true,1,0,"Xi(cc)++","\\Xi(cc)^{++}");
  InitHadron(4424,false,3.65648,0.,6,3,true,1,0,"Xi(cc)*++","\\Xi(cc)*^{++}");
  InitHadron(4432,false,3.78663,0.,3,1,true,1,0,"Omega(cc)+","\\Omega(cc)^{+}");
  InitHadron(4434,false,3.82466,0.,3,3,true,1,0,"Omega(cc)*+","\\Omega(cc)*^{+}");
  InitHadron(4444,false,4.91594,0.,6,3,true,1,0,"Omega(ccc)*++","\\Omega(ccc)*^{++}");
  InitHadron(5142,false,7.00575,0.,0,1,true,1,0,"Xi(bc)","\\Xi(bc)");
  InitHadron(5242,false,7.00575,0.,3,1,true,1,0,"Xi(bc)+","\\Xi(bc)^{+}");
  InitHadron(5342,false,7.19099,0.,0,1,true,1,0,"Omega(bc)","\\Omega(bc)");
  InitHadron(5412,false,7.03724,0.,0,1,true,1,0,"Xi(bc)'","\\Xi(bc)'");
  InitHadron(5414,false,7.0485,0.,0,3,true,1,0,"Xi(bc)*","\\Xi(bc)*");
  InitHadron(5422,false,7.03724,0.,3,1,true,1,0,"Xi(bc)'+","\\Xi(bc)'^{+}");
  InitHadron(5424,false,7.0485,0.,3,3,true,1,0,"Xi(bc)*+","\\Xi(bc)*^{+}");
  InitHadron(5432,false,7.21101,0.,0,1,true,1,0,"Omega(bc)'","\\Omega(bc)'");
  InitHadron(5434,false,7.219,0.,0,3,true,1,0,"Omega(bc)*","\\Omega(bc)*");
  InitHadron(5442,false,8.30945,0.,3,1,true,1,0,"Omega(bcc)+","\\Omega(bcc)^{+}");
  InitHadron(5444,false,8.31325,0.,3,3,true,1,0,"Omega(bcc)*+","\\Omega(bcc)*^{+}");
  InitHadron(5512,false,10.42272,0.,-3,1,true,1,0,"Xi(bb)-","\\Xi(bb)^{-}");
  InitHadron(5514,false,10.44144,0.,-3,3,true,1,0,"Xi(bb)*-","\\Xi(bb)*^{-}");
  InitHadron(5522,false,10.42272,0.,0,1,true,1,0,"Xi(bb)","\\Xi(bb)");
  InitHadron(5524,false,10.44144,0.,0,3,true,1,0,"Xi(bb)*","\\Xi(bb)*");
  InitHadron(5532,false,10.60209,0.,-3,1,true,1,0,"Omega(bb)-","\\Omega(bb)^{-}");
  InitHadron(5534,false,10.61426,0.,-3,3,true,1,0,"Omega(bb)*-","\\Omega(bb)*^{-}");
  InitHadron(5542,false,11.70767,0.,0,1,true,1,0,"Omega(bbc)","\\Omega(bbc)");
  InitHadron(5544,false,11.71147,0.,0,3,true,1,0,"Omega(bbc)*","\\Omega(bbc)*");
  InitHadron(5554,false,15.11061,0.,-3,3,true,1,0,"Omega(bbb)*-","\\Omega(bbb)*^{-}");

  // ##########################################################################
  // ##########################################################################
  // Particles NOT being used in AHADIC listed below
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
  InitHadron(kf_rho_3_1690,false,1.691,0.16,0,6,false,1,0,"rho(3)(1690)","rho_{3}(1690)");
  InitHadron(217,false,1.691,0.16,3,6,true,1,0,"rho(3)(1690)+","rho_{3}^{+}(1690)");
  InitHadron(kf_omega_3_1670,false,1.667,0.168,0,6,false,1,0,"omega(3)(1670)","omega_{3}(1670)");
  InitHadron(317,false,1.776,0.159,0,6,true,1,0,"K(3)*(1780)","K_{3}*(1780)");
  InitHadron(327,false,1.776,0.159,3,6,true,1,0,"K(3)*(1780)+","K_{3}*^{+}(1780)");
  InitHadron(kf_phi_3_1850,false,1.854,0.087,0,6,false,1,0,"phi(3)(1850)","phi_{3}(1850)");
  InitHadron(557,false,10.1599,0.0,0,6,true,1,0,"Upsilon(3)(1D)","Upsilon_{3}(1D)");
  // Tensors 4     ############################################################
  // heavy ones missing.
  InitHadron(kf_a_4_2040,false,2.014,0.361,0,8,false,1,0,"a(4)(2040)","a_{4}(2040)");
  InitHadron(219,false,2.014,0.361,3,8,true,1,0,"a(4)(2040)+","a_{4}^{+}(2040)");
  InitHadron(kf_f_4_2050,false,2.044,0.208,0,8,false,1,0,"f(4)(2050)","f_{4}(2050)");
  InitHadron(319,false,2.045,0.198,0,8,true,1,0,"K(4)*(2045)","K_{4}*(2045)");
  InitHadron(329,false,2.045,0.198,3,8,true,1,0,"K(4)*(2045)+","K_{4}*^{+}(2045)");
  // Tensors 2       ##################################################################################
  // heavy ones missing.
  InitHadron(kf_pi_2_1670,false,1.67,0.258,0,4,false,1,0,"pi(2)(1670)","pi_{2}(1670)");
  InitHadron(10215,false,1.67,0.258,3,4,true,1,0,"pi(2)(1670)+","pi_{2}^{+}(1670)");
  InitHadron(kf_eta_2_1645,false,1.617,0.181,0,4,false,1,0,"eta(2)(1645)","eta_{2}(1645)");
  InitHadron(10315,false,1.773,0.186,0,4,true,1,0,"K(2)(1770)","K_{2}(1770)");
  InitHadron(10325,false,1.773,0.186,3,4,true,1,0,"K(2)(1770)+","K_{2}^{+}(1770)");
  InitHadron(kf_eta_2_1870,false,1.842,0.225,0,4,false,1,0,"eta(2)(1870)","eta_{2}(1870)");
  InitHadron(10555,false,10.157,0.0,0,4,true,1,0,"eta(b2)(1D)","eta_{b2}(1D)");
  // Tensors 2       ##################################################################################
  // lots missing
  InitHadron(20315,false,1.816,0.276,0,4,true,1,0,"K(2)(1820)","K_{2}(1820)");
  InitHadron(20325,false,1.816,0.276,3,4,true,1,0,"K(2)(1820)+","K_{2}^{+}(1820)");
  InitHadron(20555,false,10.1562,0.0,0,4,true,1,0,"Upsilon(2)(1D)","Upsilon_{2}(1D)");
  InitHadron(30411,false,2.58,0.0,3,0,true,1,0,"D(2S)+","D(2S)^{+}");
  InitHadron(30421,false,2.58,0.0,0,0,true,1,0,"D(2S)","D(2S)");
  // Vectors 2       ##################################################################################
  // some heavy ones missing - we may have to include this because of the
  // psi(3770)
  InitHadron(kf_rho_1700,false,1.7,0.24,0,2,false,1,0,"rho(1700)","rho(1700)");
  InitHadron(30213,false,1.7,0.24,3,2,true,1,0,"rho(1700)+","rho^{+}(1700)");
  InitHadron(kf_omega_1600,false,1.670,0.31,0,2,false,1,0,"omega(1650)","omega(1650)");
  InitHadron(30313,false,1.717,0.32,0,2,true,1,0,"K*(1680)","K*(1680)");
  InitHadron(30323,false,1.717,0.32,3,2,true,1,0,"K*(1680)+","K*^{+}(1680)");
  InitHadron(kf_f1_1900_fict,false,1.900,0.32,0,2,false,1,0,"f1(1900)_fict","f1(1900)_fict");
  InitHadron(30413,false,2.64,0.0,3,2,true,1,0,"D*(2S)+","D*^{+}(2S)");
  InitHadron(30423,false,2.64,0.0,0,2,true,1,0,"D*(2S)","D*(2S)");
  InitHadron(kf_psi_3770,false,3.7699,0.0236,0,2,false,1,0,"psi(3770)","psi(3770)");
  InitHadron(30553,false,10.161,0.0,0,2,true,1,0,"Upsilon(1)(1D)","Upsilon_{1}(1D)");
  // Pseudoscalars   ##################################################################################
  // heavy ones missing.
  InitHadron(kf_pi_1300,false,1.3,0.400,0,0,false,1,0,"pi(1300)","pi(1300)");
  InitHadron(100211,false,1.3,0.400,3,0,true,1,0,"pi(1300)+","pi^{+}(1300)");
  InitHadron(kf_eta_1295,false,1.297,0.053,0,0,false,1,0,"eta(1295)","eta(1295)");
  InitHadron(100311,false,1.46,0.26,0,0,true,1,0,"K(1460)","K(1460)");
  InitHadron(100321,false,1.46,0.26,3,0,true,1,0,"K(1460)+","K^{+}(1460)");
  InitHadron(kf_eta_1475,false,1.476,0.08,0,0,false,1,0,"eta(1475)","eta(1475)");
  InitHadron(100441,false,3.638,0.014,0,0,true,1,0,"eta(c)(2S)","eta_{c}(2S)");
  InitHadron(100551,false,9.997,0.0,0,0,true,1,0,"eta(b)(2S)","eta_{b}(2S)");
  // Vectors         ##################################################################################
  // heavy ones missing.
  // the rho's may be important for tau/D/B decays
  InitHadron(kf_rho_1450,false,1.465,0.31,0,2,false,1,0,"rho(1450)","rho(1450)");
  InitHadron(100213,false,1.465,0.31,3,2,true,1,0,"rho(1450)+","rho^{+}(1450)");
  InitHadron(kf_omega_1420,false,1.419,0.17,0,2,false,1,0,"omega(1420)","omega(1420)");
  InitHadron(100313,false,1.414,0.232,0,2,true,1,0,"K*(1410)","K*(1410)");
  InitHadron(100323,false,1.414,0.232,3,2,true,1,0,"K*(1410)+","K*^{+}(1410)");
  InitHadron(kf_phi_1680,false,1.68,0.15,0,2,false,1,0,"phi(1680)","phi(1680)");
  InitHadron(kf_psi_2S,false,3.686,0.000277,0,2,false,1,0,"psi(2S)","psi(2S)");
  InitHadron(kf_Upsilon_2S,false,10.0233,4.4e-05,0,2,false,1,0,"Upsilon(2S)","Upsilon(2S)");
  // Tensors 2       ##################################################################################
  // heavy ones missing.
  InitHadron(kf_f_2_2010,false,2.01,0.2,0,4,false,1,0,"f(2)(2010)","f_{2}(2010)");
  InitHadron(100445,false,3.929,0.029,0,4,true,1,0,"chi(c2)(2P)","chi_{c2}(2P)");
  // More states without full multiplets ########################################################
  // light ones missing - we know A LOT about heavy-heavy states ....
  InitHadron(kf_chi_b2_2P,false,10.2685,0.001,0,4,false,1,0,"chi(b2)(2P)","chi_{b2}(2P)");
  InitHadron(100557,false,10.4443,0.0,0,6,true,1,0,"Upsilon(3)(2D)","Upsilon_{3}(2D)");
  InitHadron(kf_chi_b0_2P,false,10.2321,0.001,0,0,false,1,0,"chi(b0)(2P)","chi_{b0}(2P)");
  InitHadron(110553,false,10.255,0.0,0,2,true,1,0,"h(b)(2P)","h_{b}(2P)");
  InitHadron(110555,false,10.441,0.0,0,4,true,1,0,"eta(b2)(2D)","eta_{b2}(2D)");
  InitHadron(kf_chi_b1_2P,false,10.2552,0.001,0,2,false,1,0,"chi(b1)(2P)","chi_{b1}(2P)");
  InitHadron(120555,false,10.4406,0.0,0,4,true,1,0,"Upsilon(2)(2D)","Upsilon_{2}(2D)");
  InitHadron(130553,false,10.4349,0.0,0,2,true,1,0,"Upsilon(1)(2D)","Upsilon_{1}(2D)");
  InitHadron(200551,false,10.335,0.0,0,0,true,1,0,"eta(b)(3S)","eta_{b}(3S)");
  InitHadron(kf_Upsilon_3S,false,10.3553,2.63e-05,0,2,false,1,0,"Upsilon(3S)","Upsilon(3S)");
  InitHadron(200555,false,10.5264,0.0,0,4,true,1,0,"chi(b2)(3P)","chi_{b2}(3P)");
  InitHadron(210551,false,10.5007,0.0,0,0,true,1,0,"chi(b0)(3P)","chi_{b0}(3P)");
  InitHadron(210553,false,10.516,0.0,0,2,true,1,0,"h(b)(3P)","h_{b}(3P)");
  InitHadron(220553,false,10.516,0.0,0,2,true,1,0,"chi(b1)(3P)","chi_{b1}(3P)");
  InitHadron(kf_Upsilon_4S,false,10.58,0.01,0,2,false,1,0,"Upsilon(4S)","Upsilon(4S)");
  // a0(980) and friends #####################################################
  // These are the "funny" state with some iso number.  We will have to
  // figure out what to do with them.
  // #########################################################################
  InitHadron(kf_a_0_980,false,0.996,0.075,0,0,false,1,0,"a(0)(980)","a_{0}(980)");
  InitHadron(kf_f_0_600,false,0.600,0.600,0,0,false,1,0,"f(0)(600)","f_{0}(600)");
  InitHadron(kf_f_0_980,false,0.98,0.070,0,0,false,1,0,"f(0)(980)","f_{0}(980)");
  InitHadron(9000211,false,0.996,0.075,3,0,true,1,0,"a(0)(980)+","a_{0}^{+}(980)");
  InitHadron(9000311,false,0.841,0.618,0,0,true,1,0,"K(0)*(800)","K_{0}*(800)");
  InitHadron(9000321,false,0.841,0.618,3,0,true,1,0,"K(0)*(800)+","K_{0}*^{+}(800)");
  InitHadron(9000113,false,1.376,0.3,0,2,true,1,0,"pi(1)(1400)","pi_{1}(1400)");
  InitHadron(9000213,false,1.376,0.3,3,2,true,1,0,"pi(1)(1400)+","pi_{1}^{+}(1400)");
  InitHadron(9000223,false,1.518,0.073,0,2,true,1,0,"f(1)(1510)","f_{1}(1510)");
  InitHadron(9000313,false,1.65,0.15,0,2,true,1,0,"K(1)(1650)","K_{1}(1650)");
  InitHadron(9000323,false,1.65,0.15,3,2,true,1,0,"K(1)(1650)+","K_{1}^{+}(1650)");
  InitHadron(kf_psi_4040,false,4.04,0.052,0,2,false,1,0,"psi(4040)","psi(4040)");
  InitHadron(kf_Upsilon_10860,false,10.865,0.11,0,2,false,1,0,"Upsilon(10860)","Upsilon(10860)");
  InitHadron(9000115,false,1.732,0.194,0,4,true,1,0,"a(2)(1700)","a_{2}(1700)");
  InitHadron(9000215,false,1.732,0.194,3,4,true,1,0,"a(2)(1700)+","a_{2}^{+}(1700)");
  InitHadron(9000225,false,1.43,0.02,0,4,true,1,0,"f(2)(1430)","f_{2}(1430)");
  InitHadron(9000315,false,1.58,0.11,0,4,true,1,0,"K(2)(1580)","K_{2}(1580)");
  InitHadron(9000325,false,1.58,0.11,3,4,true,1,0,"K(2)(1580)+","K_{2}^{+}(1580)");
  InitHadron(9000117,false,1.982,0.188,0,6,true,1,0,"rho(3)(1990)","rho_{3}(1990)");
  InitHadron(9000217,false,1.982,0.188,3,6,true,1,0,"rho(3)(1990)+","rho_{3}^{+}(1990)");
  InitHadron(9000229,false,2.2311,0.023,0,8,true,1,0,"f(J)(2220)","f(J)(2220)");
  InitHadron(9000319,false,2.045,0.198,0,8,true,1,0,"K(4)(2500)","K_{4}(2500)");
  InitHadron(9000329,false,2.045,0.198,3,8,true,1,0,"K(4)(2500)+","K_{4}^{+}(2500)");
  InitHadron(9010111,false,1.812,0.207,0,0,true,1,0,"pi(1800)","pi(1800)");
  InitHadron(9010211,false,1.812,0.207,3,0,true,1,0,"pi(1800)+","pi^{+}(1800)");
  InitHadron(9010311,false,1.83,0.25,0,0,true,1,0,"K(1830)","K(1830)");
  InitHadron(9010321,false,1.83,0.25,3,0,true,1,0,"K(1830)+","K^{+}(1830)");
  InitHadron(9010113,false,1.653,0.225,0,2,true,1,0,"pi(1)(1600)","pi_{1}(1600)");
  InitHadron(9010213,false,1.653,0.225,3,2,true,1,0,"pi(1)(1600)+","pi_{1}^{+}(1600)");
  InitHadron(9010223,false,1.594,0.384,0,2,true,1,0,"h(1)(1595)","h_{1}(1595)");
  InitHadron(kf_Upsilon_11020,false,11.019,0.079,0,2,false,1,0,"Upsilon(11020)","Upsilon(11020)");
  InitHadron(kf_psi_4160,false,4.159,0.078,0,2,false,1,0,"psi(4160)","psi(4160)");
  InitHadron(9010115,false,2.09,0.625,0,4,true,1,0,"pi(2)(2100)","pi_{2}(2100)");
  InitHadron(9010215,false,2.09,0.625,3,4,true,1,0,"pi(2)(2100)+","pi_{2}^{+}(2100)");
  InitHadron(9010225,false,1.546,0.126,0,4,true,1,0,"f(2)(1565)","f_{2}(1565)");
  InitHadron(9010315,false,1.973,0.373,0,4,true,1,0,"K(2)*(1980)","K_{2}*(1980)");
  InitHadron(9010325,false,1.973,0.373,3,4,true,1,0,"K(2)*(1980)+","K_{2}*^{+}(1980)");
  InitHadron(9010117,false,2.25,0.2,0,6,true,1,0,"rho(3)(2250)","rho_{3}(2250)");
  InitHadron(9010217,false,2.25,0.2,3,6,true,1,0,"rho(3)(2250)+","rho_{3}^{+}(2250)");
  InitHadron(9010317,false,2.324,0.15,0,6,true,1,0,"K(3)(2320)","K_{3}(2320)");
  InitHadron(9010327,false,2.324,0.15,3,6,true,1,0,"K(3)(2320)+","K_{3}^{+}(2320)");
  InitHadron(9010229,false,2.332,0.26,0,8,true,1,0,"f(4)(2300)","f_{4}(2300)");
  InitHadron(kf_f_0_1500,false,1.4098,0.0511,0,0,false,1,0,"eta(1405)","eta(1405)");
  InitHadron(9020311,false,1.945,0.201,0,0,true,1,0,"K(0)*(1950)","K_{0}*(1950)");
  InitHadron(9020321,false,1.945,0.201,3,0,true,1,0,"K(0)*(1950)+","K_{0}*^{+}(1950)");
  InitHadron(9020113,false,1.647,0.254,0,2,true,1,0,"a(1)(1640)","a_{1}(1640)");
  InitHadron(9020213,false,1.647,0.254,3,2,true,1,0,"a(1)(1600)+","a_{1}^{+}(1600)");
  InitHadron(kf_psi_4415,false,4.415,0.043,0,2,false,1,0,"psi(4415)","psi(4415)");
  InitHadron(9020225,false,1.638,0.099,0,4,true,1,0,"f(2)(1640)","f_{2}(1640)");
  InitHadron(9020315,false,2.247,0.18,0,4,true,1,0,"K(2)(2250)","K_{2}(2250)");
  InitHadron(9020325,false,2.247,0.18,3,4,true,1,0,"K(2)(2250)+","K_{2}^{+}(2250)");
  InitHadron(kf_f_J_1710,false,1.5,0.112,0,0,false,1,0,"f(0)(1500)","f_{0}(1500)");
  InitHadron(9030113,false,1.9,0.16,0,2,true,1,0,"rho(1900)","rho(1900)");
  InitHadron(9030213,false,1.9,0.16,3,2,true,1,0,"rho(1900)+","rho^{+}(1900)");
  InitHadron(9030225,false,1.815,0.197,0,4,true,1,0,"f(2)(1810)","f_{2}(1810)");
  InitHadron(9040221,false,1.76,0.06,0,0,true,1,0,"eta(1760)","eta(1760)");
  InitHadron(9040113,false,2.149,0.363,0,2,true,1,0,"rho(2150)","rho(2150)");
  InitHadron(9040213,false,2.149,0.363,3,2,true,1,0,"rho(2150)+","rho^{+}(2150)");
  InitHadron(9040225,false,1.915,0.163,0,4,true,1,0,"f(2)(1910)","f_{2}(1910)");
  InitHadron(9050221,false,1.992,0.442,0,0,true,1,0,"f(0)(2020)","f_{0}(2020)");
  InitHadron(kf_f_2_2300,false,1.944,0.472,0,4,false,1,0,"f(2)(1950)","f_{2}(1950)");
  InitHadron(9060221,false,2.103,0.206,0,0,true,1,0,"f(0)(2100)","f_{0}(2100)");
  InitHadron(kf_f_2_2340,false,2.011,0.202,0,4,false,1,0,"f(2)(2010)","f_{2}(2010)");
  InitHadron(9070221,false,2.189,0.238,0,0,true,1,0,"f(0)(2200)","f_{0}(2200)");
  InitHadron(9070225,false,2.15,0.2,0,4,true,1,0,"f(2)(2150)","f_{2}(2150)");
  InitHadron(9080221,false,2.220,0.15,0,0,true,1,0,"eta(2225)","eta(2225)");
  InitHadron(9080225,false,2.297,0.15,0,4,true,1,0,"f(2)(2300)","f_{2}(2300)");
  InitHadron(9090225,false,2.34,0.32,0,4,true,1,0,"f(2)(2340)","f_{2}(2340)");

  // ##########################################################################
  // ##########################################################################
  // Former members of the Sherpa team - made immortal as particles here ######
  // ##########################################################################
  // ##########################################################################
  InitHadron(5505,false,1000000.0,1000,0,0,true,0,0,"ralf","ralf");
  InitHadron(5506,false,1000000.0,1000,0,0,true,0,0,"ande","ande");
  InitHadron(5507,false,1000000.0,1000,0,0,true,0,0,"thomas","thomas");
  InitHadron(5508,false,1000000.0,1000,0,0,true,0,0,"tanju","tanju");
  InitHadron(5509,false,1000000.0,1000,0,0,true,0,0,"jennifer","jennifer");
  InitHadron(6505,false,1000000.0,1000,0,0,true,0,0,"hendrik","hendrik");
  InitHadron(6506,false,1000000.0,1000,0,0,true,0,0,"jan","jan");
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################

  OverrideProperties();
}

void Hadron_Init::InitHadron(const kf_code& kfc,
                             const bool checkinitialised,
                             const double &mass,
                             const double &width,
                             const int icharge,
                             const int spin,
                             const bool majorana,
                             const bool on,
                             const int stable,
                             const std::string& idname,
                             const std::string& texname)
{
  if (checkinitialised) {
    // if not initialised in SHRIMPS, BEAM, AMISIC, ...
    if(s_kftable.find(kfc) != s_kftable.end())
      return;
  }
  s_kftable[kfc]=new Particle_Info(kfc, mass, width, icharge, spin, on, stable,
                                   idname, texname);
  if (!majorana)
    s_kftable[kfc]->m_majorana=-1;

  m_addedhadrons.insert(kfc);
}

void Hadron_Init::InitHadron(const kf_code& kfc,
                             const double &mass,
                             const double &width,
                             const int icharge,
                             const int strong,
                             const int spin,
                             const bool majorana,
                             const bool on,
                             const int stable,
                             const bool massive,
                             const std::string& idname,
                             const std::string& antiname,
                             const std::string& texname,
                             const std::string& antitexname)
{
  s_kftable[kfc]=new Particle_Info(
      kfc, mass, width, icharge, strong, spin, majorana, on, stable, massive,
      idname, antiname, texname, antitexname);
  m_addedhadrons.insert(kfc);
}

void Hadron_Init::OverrideProperties()
{
  auto pdata = Settings::GetMainSettings()["PARTICLE_DATA"];
  for (const auto& ptclname : pdata.GetKeys()) {
    kf_code kf = ToType<kf_code>(ptclname);
    if (m_addedhadrons.find(kf) == m_addedhadrons.end())
      continue;
    const auto it = s_kftable.find(kf);
    if (it != s_kftable.end()) {
      for (const auto& propertyname : pdata[ptclname].GetKeys()) {
        auto s = pdata[ptclname][propertyname];
        if (propertyname == "Mass") {
          it->second->m_mass = s.SetDefault(it->second->m_mass).Get<double>();
        } else if (propertyname == "Width") {
          it->second->m_width = s.SetDefault(it->second->m_width).Get<double>();
        } else if (propertyname == "Active") {
          it->second->m_on = s.SetDefault(it->second->m_on).Get<bool>();
        } else if (propertyname == "Stable") {
          it->second->m_stable = s.SetDefault(it->second->m_stable).Get<int>();
        } else if (propertyname == "Massive") {
          it->second->m_massive = s.SetDefault(it->second->m_massive).Get<bool>();
        }
      }
    }
  }
}
