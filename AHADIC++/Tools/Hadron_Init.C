#include "AHADIC++/Tools/Hadron_Init.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

void Hadron_Init::Init() {
  msg_Info()<<METHOD<<"(): Initializing kf table for hadrons.\n";
  // ##########################################################################
  // ##########################################################################
  // Particle being used in AHADIC listed below
  // ##########################################################################
  // ##########################################################################
  // Diquarks - not exactly particles but we need them ########################
  // ##########################################################################
  // ##########################################################################
  if(s_kftable.find(1103)==s_kftable.end()) { // if not initialized in Beam
    s_kftable[1103]=new Particle_Info(1103,0.77133,0,-2,-3,2,0,0,1,1,
				      "dd_1","dd_1b","dd_1b","dd_1b");
    s_kftable[2101]=new Particle_Info(2101,0.57933,0,1,-3,0,0,0,1,1,
				      "ud_0","ud_0b","ud_0b","ud_0b");
    s_kftable[2103]=new Particle_Info(2103,0.77133,0,1,-3,2,0,0,1,1,
				      "ud_1","ud_1b","ud_1b","ud_1b");
    s_kftable[2203]=new Particle_Info(2203,0.77133,0,4,-3,2,0,0,1,1,
				      "uu_1","uu_1b","uu_1b","uu_1b");
    s_kftable[3101]=new Particle_Info(3101,0.80473,0,-2,-3,0,0,0,1,1,
				      "sd_0","sd_0b","sd_0b","sd_0b");
    s_kftable[3103]=new Particle_Info(3103,0.92953,0,-2,-3,2,0,0,1,1,
				      "sd_1","sd_1b","sd_1b","sd_1b");
    s_kftable[3201]=new Particle_Info(3201,0.80473,0,1,-3,0,0,0,1,1,
				      "su_0","su_0b","su_0b","su_0b");
    s_kftable[3203]=new Particle_Info(3203,0.92953,0,1,-3,2,0,0,1,1,
				      "su_1","su_1b","su_1b","su_1b");
    s_kftable[3303]=new Particle_Info(3303,1.09361,0,-2,-3,2,0,0,1,1,
				      "ss_1","ss_1b","ss_1b","ss_1b");
    s_kftable[4101]=new Particle_Info(4101,1.96908,0,1,-3,0,0,0,1,1,
				      "cd_0","cd_0b","cd_0b","cd_0b");
    s_kftable[4103]=new Particle_Info(4103,2.00808,0,1,-3,2,0,0,1,1,
				      "cd_1","cd_1b","cd_1b","cd_1b");
    s_kftable[4201]=new Particle_Info(4201,1.96908,0,4,-3,0,0,0,1,1,
				      "cu_0","cu_0b","cu_0b","cu_0b");
    s_kftable[4203]=new Particle_Info(4203,2.00808,0,4,-3,2,0,0,1,1,
				      "cu_1","cu_1b","cu_1b","cu_1b");
    s_kftable[4301]=new Particle_Info(4301,2.15432,0,1,-3,0,0,0,1,1,
				      "cs_0","cs_0b","cs_0b","cs_0b");
    s_kftable[4303]=new Particle_Info(4303,2.17967,0,1,-3,2,0,0,1,1,
				      "cs_1","cs_1b","cs_1b","cs_1b");
    s_kftable[4403]=new Particle_Info(4403,3.27531,0,4,-3,2,0,0,1,1,
				      "cc_1","cc_1b","cc_1b","cc_1b");
    s_kftable[5101]=new Particle_Info(5101,5.38897,0,-2,-3,0,0,0,1,1,
				      "bd_0","bd_0b","bd_0b","bd_0b");
    s_kftable[5103]=new Particle_Info(5103,5.40145,0,-2,-3,2,0,0,1,1,
				      "bd_1","bd_1b","bd_1b","bd_1b");
    s_kftable[5201]=new Particle_Info(5201,5.38897,0,1,-3,0,0,0,1,1,
				      "bu_0","bu_0b","bu_0b","bu_0b");
    s_kftable[5203]=new Particle_Info(5203,5.40145,0,1,-3,2,0,0,1,1,
				      "bu_1","bu_1b","bu_1b","bu_1b");
    s_kftable[5301]=new Particle_Info(5301,5.56725,0,-2,-3,0,0,0,1,1,
				      "bs_0","bs_0b","bs_0b","bs_0b");
    s_kftable[5303]=new Particle_Info(5303,5.57536,0,-2,-3,2,0,0,1,1,
				      "bs_1","bs_1b","bs_1b","bs_1b");
    s_kftable[5401]=new Particle_Info(5401,6.67143,0,1,-3,0,0,0,1,1,
				      "bc_0","bc_0b","bc_0b","bc_0b");
    s_kftable[5403]=new Particle_Info(5403,6.67397,0,1,-3,2,0,0,1,1,
				      "bc_1","bc_1b","bc_1b","bc_1b");
    s_kftable[5503]=new Particle_Info(5503,10.07354,0,-2,-3,2,0,0,1,1,
				      "bb_1","bb_1b","bb_1b","bb_1b");
  }
  // ##########################################################################
  // ##########################################################################
  // Other objects that may show up in an event record ########################
  s_kftable[kf_cluster]=new Particle_Info(kf_cluster,0.,0.,0,0,0,0,1,1,0,
					  "cluster","cluster",
					  "cluster","cluster");
  s_kftable[kf_string] =new Particle_Info (kf_string,0.,0.,0,0,0,0,1,1,0,
					   "string","string",
					   "string","string");

  // ##########################################################################
  // ##########################################################################
  // MESON MULTIPLETS
  // ##########################################################################
  // Pseudoscalars ############################################################
  // kfc, mass, width, icharge, spin2, on, stable, idname, texname
  if(s_kftable.find(111)==s_kftable.end()) // if not initialized in amisic
  s_kftable[111]=new Particle_Info(111,0.134976,7.8486e-09,0,0,1,0,
				   "pi","pi");
  if(s_kftable.find(211)==s_kftable.end()) // if not initialized yet
  s_kftable[211]=new Particle_Info(211,0.13957,2.5242e-17,3,0,1,1,
				   "pi+","pi^{+}");
  s_kftable[221]=new Particle_Info(221,0.5473,1.18e-06,0,0,1,0,
				   "eta","eta");
  s_kftable[311]=new Particle_Info(311,0.49767,1.e-16,0,0,1,0,
				   "K","K");
  s_kftable[130]=new Particle_Info(130,0.49767,1.273e-17,0,0,1,1,
				   "K(L)","K_{L}");
  s_kftable[310]=new Particle_Info(310,0.49767,7.373e-15,0,0,1,0,
				   "K(S)","K_{S}");
  s_kftable[321]=new Particle_Info(321,0.493677,5.314e-17,3,0,1,1,
				   "K+","K^{+}");
  s_kftable[331]=new Particle_Info(331,0.95778,0.000203,0,0,1,0,
				   "eta'(958)","eta'(958)");
  s_kftable[411]=new Particle_Info(411,1.8693,6.23e-13,3,0,1,0,
				   "D+","D^{+}");
  s_kftable[421]=new Particle_Info(421,1.8646,1.586e-12,0,0,1,0,
				   "D","D");
  s_kftable[431]=new Particle_Info(431,1.9685,1.41e-12,3,0,1,0,
				   "D(s)+","D(s)^{+}");
  s_kftable[441]=new Particle_Info(441,2.9798,0.0132,0,0,1,0,
				   "eta(c)(1S)","eta_{c}(1S)");
  s_kftable[511]=new Particle_Info(511,5.2792,4.22e-13,0,0,1,0,
				   "B","B");
  s_kftable[521]=new Particle_Info(521,5.2789,3.99e-13,3,0,1,0,
				   "B+","B^{+}");
  s_kftable[531]=new Particle_Info(531,5.3693,4.27e-13,0,0,1,0,
				   "B(s)","B_{s}");
  s_kftable[541]=new Particle_Info(541,6.4,1.43e-12,3,0,1,0,
				   "B(c)+","B_{c}^{+}");
  s_kftable[551]=new Particle_Info(551,9.4,0.050,0,0,1,0,
				   "eta(b)(1S)","eta_{b}(1S)");
  // Vectors       ############################################################
  s_kftable[113]=new Particle_Info(113,0.77,0.1507,0,2,1,0,
				   "rho(770)","rho(770)");
  s_kftable[213]=new Particle_Info(213,0.77,0.1507,3,2,1,0,
				   "rho(770)+","rho^{+}(770)");
  s_kftable[223]=new Particle_Info(223,0.78194,0.00841,0,2,1,0,
				   "omega(782)","omega(782)");
  s_kftable[313]=new Particle_Info(313,0.8961,0.0505,0,2,1,0,
				   "K*(892)","K*(892)");
  s_kftable[323]=new Particle_Info(323,0.89166,0.0508,3,2,1,0,
				   "K*(892)+","K*^{+}(892)");
  s_kftable[333]=new Particle_Info(333,1.01941,0.00443,0,2,1,0,
				   "phi(1020)","phi(1020)");
  s_kftable[413]=new Particle_Info(413,2.01,0.001,3,2,1,0,
				   "D*(2010)+","D*^{+}(2010)");
  s_kftable[423]=new Particle_Info(423,2.0067,0.001,0,2,1,0,
				   "D*(2007)","D*(2007)");
  s_kftable[433]=new Particle_Info(433,2.1124,0.001,3,2,1,0,
				   "D(s)*+","D_{s}*^{+}");
  s_kftable[443]=new Particle_Info(443,3.09688,8.7e-05,0,2,1,0,
				   "J/psi(1S)","J/psi(1S)");
  s_kftable[513]=new Particle_Info(513,5.3249,0.0001,0,2,1,0,
				   "B*","B*");
  s_kftable[523]=new Particle_Info(523,5.3249,0.0001,3,2,1,0,
				   "B*+","B*^{+}");
  s_kftable[533]=new Particle_Info(533,5.41630,0.0001,0,2,1,0,
				   "B(s)*","B_{s}*");
  s_kftable[543]=new Particle_Info(543,6.602,0.0001,3,2,1,0,
				   "B(c)*+","B_{c}*^{+}");
  s_kftable[553]=new Particle_Info(553,9.46037,5.25e-05,0,2,1,0,
				   "Upsilon(1S)","Upsilon(1S)");
  // Tensors 2     ############################################################
  s_kftable[115]=new Particle_Info(115,1.3181,0.107,0,4,1,0,
				   "a(2)(1320)","a_{2}(1320)");
  s_kftable[215]=new Particle_Info(215,1.3181,0.107,3,4,1,0,
				   "a(2)(1320)+","a_{2}^{+}(1320)");
  s_kftable[225]=new Particle_Info(225,1.275,0.1855,0,4,1,0,
				   "f(2)(1270)","f_{2}(1270)");
  s_kftable[315]=new Particle_Info(315,1.4324,0.109,0,4,1,0,
				   "K(2)*(1430)","K_{2}*(1430)");
  s_kftable[325]=new Particle_Info(325,1.4256,0.0985,3,4,1,0,
				   "K(2)*(1430)+","K_{2}*^{+}(1430)");
  s_kftable[335]=new Particle_Info(335,1.525,0.076,0,4,1,0,
				   "f(2)'(1525)","f_{2}'(1525)");
  s_kftable[415]=new Particle_Info(415,2.459,0.025,3,4,1,0,
				   "D(2)*(2460)+","D_{2}*^{+}(2460)");
  s_kftable[425]=new Particle_Info(425,2.4589,0.023,0,4,1,0,
				   "D(2)*(2460)","D_{2}*(2460)");
  s_kftable[435]=new Particle_Info(435,2.5735,0.015,3,4,1,0,
				   "D(s2)*(2573)+","D_{s2}*^{+}(2573)");
  s_kftable[445]=new Particle_Info(445,3.55617,0.002,0,4,1,0,
				   "chi(c2)(1P)","chi_{c2}(1P)");
  s_kftable[515]=new Particle_Info(515,5.83,0.02,0,4,1,0,
				   "B(2)*","B_{2}*");
  s_kftable[525]=new Particle_Info(525,5.83,0.02,3,4,1,0,
				   "B(2)*+","B_{2}*^{+}");
  s_kftable[535]=new Particle_Info(535,5.8397,0.02,0,4,1,0,
				   "B(s2)*","B_{s2}*");
  s_kftable[545]=new Particle_Info(545,7.35,0.02,3,4,1,0,
				   "B(c2)*+","B_{c2}*^{+}");
  s_kftable[555]=new Particle_Info(555,9.9132,0.001,0,4,1,0,
				   "chi(b2)(1P)","chi_{b2}(1P)");
  // Scalars       ############################################################
  s_kftable[10111]=new Particle_Info(10111,1.474,0.265,0,0,1,0,
				     "a(0)(1450)","a_{0}(1450)");
  s_kftable[10211]=new Particle_Info(10211,1.474,0.265,3,0,1,0,
				     "a(0)(1450)+","a_{0}^{+}(1450)");
  s_kftable[10221]=new Particle_Info(10221,1.4,0.5,0,0,1,0,
				     "f(0)(1370)","f_{0}(1370)");
  s_kftable[10311]=new Particle_Info(10311,1.429,0.287,0,0,1,0,
				     "K(0)*(1430)","K_{0}*(1430)");
  s_kftable[10321]=new Particle_Info(10321,1.429,0.287,3,0,1,0,
				     "K(0)*(1430)+","K_{0}*^{+}(1430)");
  s_kftable[10331]=new Particle_Info(10331,1.7,0.5,0,0,1,0,
				     "f(0)(1710)","f_{0}(1710)");
  s_kftable[10411]=new Particle_Info(10411,2.26511,0.05,3,0,1,0,
				     "D(0)*(2400)+","D_{0}*^{+}(2400)");
  s_kftable[10421]=new Particle_Info(10421,2.272,0.05,0,0,1,0,
				     "D(0)*(2400)","D_{0}*(2400)");
  s_kftable[10431]=new Particle_Info(10431,2.3178,0.05,3,0,1,0,
				     "D(s0)*(2317)+","D_{s0}*^{+}(2317)");
  s_kftable[10441]=new Particle_Info(10441,3.4173,0.014,0,0,1,0,
				     "chi(c0)(1P)","chi_{c0}(1P)");
  s_kftable[10511]=new Particle_Info(10511,5.68,0.05,0,0,1,0,
				     "B(0)*","B_{0}*");
  s_kftable[10521]=new Particle_Info(10521,5.68,0.05,3,0,1,0,
				     "B(0)*+","B_{0}*^{+}");
  s_kftable[10531]=new Particle_Info(10531,5.92,0.05,0,0,1,0,
				     "B(s0)*","B_{s0}*");
  s_kftable[10541]=new Particle_Info(10541,7.25,0.05,3,0,1,0,
				     "B(c0)*+","B_{c0}*^{+}");
  s_kftable[10551]=new Particle_Info(10551,9.8598,0.05,0,0,1,0,
				     "chi(b0)(1P)","chi_{b0}(1P)");
  // Axial vectors x############################################################
  s_kftable[10113]=new Particle_Info(10113,1.2295,0.142,0,2,1,0,
				     "b(1)(1235)","b_{1}(1235)");
  s_kftable[10213]=new Particle_Info(10213,1.2295,0.142,3,2,1,0,
				     "b(1)(1235)+","b_{1}^{+}(1235)");
  s_kftable[10223]=new Particle_Info(10223,1.17,0.36,0,2,1,0,
				     "h(1)(1170)","h_{1}(1170)");
  s_kftable[10313]=new Particle_Info(10313,1.272,0.09,0,2,1,0,
				     "K(1)(1270)","K_{1}(1270)");
  s_kftable[10323]=new Particle_Info(10323,1.272,0.09,3,2,1,0,
				     "K(1)(1270)+","K_{1}^{+}(1270)");
  s_kftable[10333]=new Particle_Info(10333,1.386,0.091,0,2,1,0,
				     "h(1)(1380)","h_{1}(1380)");
  s_kftable[10413]=new Particle_Info(10413,2.424,0.02,3,2,1,0,
				     "D(1)(2420)+","D_{1}^{+}(2420)");
  s_kftable[10423]=new Particle_Info(10423,2.4222,0.0189,0,2,1,0,
				     "D(1)(2420)","D_{1}(2420)");
  s_kftable[10433]=new Particle_Info(10433,2.5353,0.001,3,2,1,0,
				     "D(s1)(2536)+","D_{s1}^{+}(2536)");
  s_kftable[10443]=new Particle_Info(10443,3.46,0.01,0,2,1,0,
				     "h(c)(1P)","h_{c}(1P)");
  s_kftable[10513]=new Particle_Info(10513,5.73,0.05,0,2,1,0,
				     "B(1)(L)","B_{1}(L)");
  s_kftable[10523]=new Particle_Info(10523,5.73,0.05,3,2,1,0,
				     "B(1)(L)+","B_{1}^{+}(L)");
  s_kftable[10533]=new Particle_Info(10533,5.97,0.05,0,2,1,0,
				     "B(s1)(L)","B_{s1}(L)");
  s_kftable[10543]=new Particle_Info(10543,7.3,0.05,3,2,1,0,
				     "B(c1)(L)+","B_{c1}^{+}(L)");
  s_kftable[10553]=new Particle_Info(10553,9.875,0.01,0,2,1,0,
				     "h(b)(1P)","h_{b}(1P)");
  // Vectors       ############################################################
  s_kftable[20113]=new Particle_Info(20113,1.23,0.400,0,2,1,0,
				     "a(1)(1260)","a_{1}(1260)");
  s_kftable[20213]=new Particle_Info(20213,1.23,0.400,3,2,1,0,
				     "a(1)(1260)+","a_{1}^{+}(1260)");
  s_kftable[20223]=new Particle_Info(20223,1.2819,0.024,0,2,1,0,
				     "f(1)(1285)","f_{1}(1285)");
  s_kftable[20313]=new Particle_Info(20313,1.402,0.174,0,2,1,0,
				     "K(1)(1400)","K_{1}(1400)");
  s_kftable[20323]=new Particle_Info(20323,1.402,0.174,3,2,1,0,
				     "K(1)(1400)+","K_{1}^{+}(1400)");
  s_kftable[20333]=new Particle_Info(20333,1.4262,0.055,0,2,1,0,
				     "f(1)(1420)","f_{1}(1420)");
  s_kftable[20413]=new Particle_Info(20413,2.39155,0.003,3,2,1,0,
				     "D(1)(H)+","D_{1}^{+}(H)");
  s_kftable[20423]=new Particle_Info(20423,2.43267,0.003,0,2,1,0,
				     "D(1)(2430)","D_{1}(2430)");
  s_kftable[20433]=new Particle_Info(20433,2.4596,0.003,3,2,1,0,
				     "D(s1)(2460)+","D_{s1}^{+}(2460)");
  s_kftable[20443]=new Particle_Info(20443,3.51053,0.00088,0,2,1,0,
				     "chi(c1)(1P)","chi_{c1}(1P)");
  s_kftable[20513]=new Particle_Info(20513,5.78,0.05,0,2,1,0,
				     "B(1)(H)","B_{1}(H)");
  s_kftable[20523]=new Particle_Info(20523,5.78,0.05,3,2,1,0,
				     "B(1)(H)+","B_{1}^{+}(H)");
  s_kftable[20533]=new Particle_Info(20533,6.02,0.05,0,2,1,0,
				     "B(s1)(H)","B_{s1}(H)");
  s_kftable[20543]=new Particle_Info(20543,7.3,0.05,3,2,1,0,
				     "B(c1)(H)+","B_{c1}^{+}(H)");
  s_kftable[20553]=new Particle_Info(20553,9.8919,0.001,0,2,1,0,
				     "chi(b1)(1P)","chi_{b1}(1P)");

  // #########################################################################
  // BARYON MULTIPLETS
  // #########################################################################
  // Nucleons (octet) ########################################################
  if(s_kftable.find(2112)==s_kftable.end()) // if not initialized in Beam
    s_kftable[2112]=new Particle_Info(2112,0.939566,7.424e-28,0,1,1,1,
				      "n","n");
  if(s_kftable.find(2212)==s_kftable.end()) // if not initialized in Beam
    s_kftable[2212]=new Particle_Info(2212,0.938272,0,3,1,1,1,
				      "P+","P^{+}");
  s_kftable[3112]=new Particle_Info(3112,1.19745,4.45e-15,-3,1,1,0,
				    "Sigma-","\\Sigma^{-}");
  s_kftable[3212]=new Particle_Info(3212,1.19264,8.9e-06,0,1,1,0,
				    "Sigma","\\Sigma");
  s_kftable[3222]=new Particle_Info(3222,1.18937,8.24e-15,3,1,1,0,
				    "Sigma+","\\Sigma^{+}");
  s_kftable[3122]=new Particle_Info(3122,1.11568,2.501e-15,0,1,1,0,
				    "Lambda","\\Lambda");
  s_kftable[3312]=new Particle_Info(3312,1.32132,4.02e-15,-3,1,1,0,
				    "Xi-","\\Xi^{-}");
  s_kftable[3322]=new Particle_Info(3322,1.3149,2.27e-15,0,1,1,0,
				    "Xi","\\Xi");
  s_kftable[4112]=new Particle_Info(4112,2.4522,0.0022,0,1,1,0,
				    "Sigma(c)(2455)","\\Sigma_{c}(2455)");
  s_kftable[4212]=new Particle_Info(4212,2.4536,0.0023,3,1,1,0,
				    "Sigma(c)(2455)+","\\Sigma_{c}^{+}(2455)");
  s_kftable[4222]=new Particle_Info(4222,2.4528,0.0023,6,1,1,0,
				  "Sigma(c)(2455)++","\\Sigma_{c}^{++}(2455)");
  s_kftable[4122]=new Particle_Info(4122,2.2849,3.19e-12,3,1,1,0,
				    "Lambda(c)+","\\Lambda_{c}^{+}");
  s_kftable[4132]=new Particle_Info(4132,2.4703,5.875e-12,0,1,1,0,
				    "Xi(c)","\\Xi_{c}");
  s_kftable[4232]=new Particle_Info(4232,2.4656,1.489e-12,3,1,1,0,
				    "Xi(c)+","\\Xi_{c}^{+}");
  s_kftable[4312]=new Particle_Info(4312,2.575,0.001,0,1,1,0,
				    "Xi(c)'","\\Xi'_{c}");
  s_kftable[4322]=new Particle_Info(4322,2.578,0.001,3,1,1,0,
				    "Xi(c)'+","\\Xi'_{c}^{+}");
  s_kftable[4332]=new Particle_Info(4332,2.704,1.02e-11,0,1,1,0,
				    "Omega(c)","\\Omega_{c}");
  s_kftable[5112]=new Particle_Info(5112,5.810,0.001,-3,1,1,0,
				    "Sigma(b)-","\\Sigma_{b}^{-}");
  s_kftable[5212]=new Particle_Info(5212,5.810,0.001,0,1,1,0,
				    "Sigma(b)","\\Sigma_{b}");
  s_kftable[5222]=new Particle_Info(5222,5.810,0.001,3,1,1,0,
				    "Sigma(b)+","\\Sigma_{b}^{+}");
  s_kftable[5122]=new Particle_Info(5122,5.624,5.31e-13,0,1,1,0,
				    "Lambda(b)","\\Lambda_{b}");
  s_kftable[5132]=new Particle_Info(5132,5.790,1.e-13,-3,1,1,0,
				    "Xi(b)-","\\Xi_{b}^{-}");
  s_kftable[5232]=new Particle_Info(5232,5.790,1.e-13,0,1, 1,0,
				    "Xi(b)","\\Xi_{b}");
  s_kftable[5312]=new Particle_Info(5312,5.890,1.e-13,-3,1,1,0,
				    "Xi(b)'-_fict","\\Xi'_{b}_fict^{-}");
  s_kftable[5322]=new Particle_Info(5322,5.890,1.e-13,0,1,1,0,
				    "Xi(b)'_fict","\\Xi'_{b}_fict");
  s_kftable[5332]=new Particle_Info(5332,6.071,1.e-13,-3,1,1,0,
				    "Omega(b)-","\\Omega_{b}^{-}");
  // Deltas (decuplet) ########################################################
  s_kftable[1114]=new Particle_Info(1114,1.232,0.12,-3,3,1,0,
				    "Delta(1232)-","\\Delta-(1232)");
  s_kftable[2114]=new Particle_Info(2114,1.232,0.12,0,3,1,0,
				    "Delta(1232)","\\Delta(1232)");
  s_kftable[2214]=new Particle_Info(2214,1.232,0.12,3,3,1,0,
				    "Delta(1232)+","\\Delta^{+}(1232)");
  s_kftable[2224]=new Particle_Info(2224,1.232,0.12,6,3,1,0,
				    "Delta(1232)++","\\Delta^{++}(1232)");
  s_kftable[3114]=new Particle_Info(3114,1.3872,0.0394,-3,3,1,0,
				    "Sigma(1385)-","\\Sigma-(1385)");
  s_kftable[3214]=new Particle_Info(3214,1.3837,0.036,0,3,1,0,
				    "Sigma(1385)","\\Sigma(1385)");
  s_kftable[3224]=new Particle_Info(3224,1.3828,0.0358,3,3,1,0,
				    "Sigma(1385)+","\\Sigma^{+}(1385)");
  s_kftable[3314]=new Particle_Info(3314,1.535,0.0099,-3,3,1,0,
				    "Xi(1530)-","\\Xi-(1530)");
  s_kftable[3324]=new Particle_Info(3324,1.5318,0.0091,0,3,1,0,
				    "Xi(1530)","\\Xi(1530)");
  s_kftable[3334]=new Particle_Info(3334,1.67245,8.01e-15,-3,3,1,0,
				    "Omega-","\\Omega^{-}");
  s_kftable[4114]=new Particle_Info(4114,2.5175,0.0150,0,3,1,0,
				    "Sigma(c)(2520)","\\Sigma_{c}(2520)");
  s_kftable[4214]=new Particle_Info(4214,2.5159,0.0150,3,3,1,0,
				    "Sigma(c)(2520)+","\\Sigma_{c}^{+}(2520)");
  s_kftable[4224]=new Particle_Info(4224,2.5194,0.0150,6,3,1,0,
				  "Sigma(c)(2520)++","\\Sigma_{c}^{++}(2520)");
  s_kftable[4314]=new Particle_Info(4314,2.645,0.003,0,3,1,0,
				    "Xi(c)*","\\Xi*_{c}");
  s_kftable[4324]=new Particle_Info(4324,2.645,0.003,3,3,1,0,
				    "Xi(c)*+","\\Xi*_{c}^{+}");
  s_kftable[4334]=new Particle_Info(4334,2.8,0.001,0,3,1,0,
				    "Omega(c)*","\\Omega*_{c}");
  s_kftable[5114]=new Particle_Info(5114,5.829,0.01,-3,3,1,0,
				    "Sigma(b)*-","\\Sigma*_{b}^{-}");
  s_kftable[5214]=new Particle_Info(5214,5.829,0.01,0,3,1,0,
				    "Sigma(b)*","\\Sigma*_{b}");
  s_kftable[5224]=new Particle_Info(5224,5.829,0.01,3,3,1,0,
				    "Sigma(b)*+","\\Sigma*_{b}^{+}");
  s_kftable[5314]=new Particle_Info(5314,5.930,0.001,-3,3,1,0,
				    "Xi(b)*-_fict","\\Xi*_{b}_fict^{-}");
  s_kftable[5324]=new Particle_Info(5324,5.930,0.001,0,3,1,0,
				    "Xi(b)*_fict","\\Xi*_{b}_fict");
  s_kftable[5334]=new Particle_Info(5334,6.090,0.0003,-3,3,1,0,
				    "Omega(b)*-_fict","\\Omega*_{b}_fict^{-}");
  // Nucleons (octet) - L_N = 1 ##############################################
  // careful - we will have to add some heavy states here!
  s_kftable[102112]=new Particle_Info(102112,1.535,0.15,0,1,1,0,
				      "N(1535)","N(1535)");
  s_kftable[102212]=new Particle_Info(102212,1.535,0.15,3,1,1,0,
				      "N(1535)+","N^{+}(1535)");
  s_kftable[102132]=new Particle_Info(102132,1.407,0.05,0,1,1,0,
				      "Lambda(1405)","\\Lambda(1405)");
  s_kftable[103122]=new Particle_Info(103122,1.67,0.06,0,1,1,0,
				      "Lambda(1670)","\\Lambda(1670)");
  s_kftable[103112]=new Particle_Info(103112,1.75,0.09,-3,1,1,0,
				      "Sigma(1620)-","\\Sigma-(1620)");
  s_kftable[103212]=new Particle_Info(103212,1.75,0.09,0,1,1,0,
				      "Sigma(1620)","\\Sigma(1620)");
  s_kftable[103222]=new Particle_Info(103222,1.75,0.09,3,1,1,0,
				      "Sigma(1620)+","\\Sigma^{+}(1620)");
  s_kftable[103312]=new Particle_Info(103312,1.75,0.09,-3,1,1,0,
				      "Xi(1750)","\\Xi^{-}(1750)");
  s_kftable[103322]=new Particle_Info(103322,1.75,0.09,0,1,1,0,
				      "Xi(1750)","\\Xi(1750)");
  s_kftable[102142]=new Particle_Info(102142,2.5954,0.0036,3,1,1,0,
				      "Lambda(c)(2595)+","\\Lambda_{c}(2595)^{+}");
  // Nucleons (octet) - the Roper resonance - L_N = 2 #########################
  // careful - we will have to add some heavy states here!
  if(s_kftable.find(202112)==s_kftable.end()) // if not initialised in SHRIMPS
    s_kftable[202112]=new Particle_Info(202112,1.44,0.35,0,1,1,0,
					"N(1440)","N(1440)");
  if(s_kftable.find(202212)==s_kftable.end()) // if not initialised in SHRIMPS 
    s_kftable[202212]=new Particle_Info(202212,1.44,0.35,3,1,1,0,
					"N(1440)+","N^{+}(1440)");
  s_kftable[203112]=new Particle_Info(203112,1.66,0.1,-3,1,1,0,
				      "Sigma(1660)-","\\Sigma-(1660)");
  s_kftable[203212]=new Particle_Info(203212,1.66,0.1,0,1,1,0,
				      "Sigma(1660)","\\Sigma(1660)");
  s_kftable[203222]=new Particle_Info(203222,1.66,0.1,3,1,1,0,
				      "Sigma(1660)+","\\Sigma^{+}(1660)");
  s_kftable[203122]=new Particle_Info(203122,1.6,0.15,0,1,1,0,
				      "Lambda(1600)","\\Lambda(1600)");
  s_kftable[203312]=new Particle_Info(203312,1.696,0.010,-3,1,1,0,
				      "Xi(1690)-","\\Xi-(1690)");
  s_kftable[203322]=new Particle_Info(203322,1.696,0.010,0,1,1,0,
				      "Xi(1690)","\\Xi(1690)");
  // Nucleons (octet) - L_N = 1 -- plus singlet ##############################
  // careful - we will have to add some heavy states here!
  s_kftable[102114]=new Particle_Info(102114,1.52,0.12,0,3,1,0,
				      "N(1520)","N(1520)");
  s_kftable[102214]=new Particle_Info(102214,1.52,0.12,3,3,1,0,
				      "N(1520)+","N^{+}(1520)");
  s_kftable[103114]=new Particle_Info(103114,1.67,0.06,-3,3,1,0,
				      "Sigma(1670)-","\\Sigma-(1670)");
  s_kftable[103214]=new Particle_Info(103214,1.67,0.06,0,3,1,0,
				      "Sigma(1670)","\\Sigma(1670)");
  s_kftable[103224]=new Particle_Info(103224,1.67,0.06,3,3,1,0,
				      "Sigma(1670)+","\\Sigma^{+}(1670)");
  s_kftable[103124]=new Particle_Info(103124,1.69,0.06,0,3,1,0,
				      "Lambda(1690)","\\Lambda(1690)");
  s_kftable[103314]=new Particle_Info(103314,1.823,0.024,-3,3,1,0,
				      "Xi(1820)-","\\Xi-(1820)");
  s_kftable[103324]=new Particle_Info(103324,1.823,0.024,0,3,1,0,
				      "Xi(1820)","\\Xi(1820)");
  s_kftable[102134]=new Particle_Info(102134,1.5195,0.0156,0,3,1,0,
				      "Lambda(1520)","\\Lambda(1520)");
  s_kftable[102144]=new Particle_Info(104214,2.625,0.002,3,3,1,0,
				  "Lambda(c)(2625)+","\\Lambda_{c}^{+}(2625)");
  s_kftable[104312]=new Particle_Info(104312,2.815,0.002,0,3,1,0,
				      "Xi(c)(2815)","\\Xi_{c}(2815)");
  s_kftable[104322]=new Particle_Info(104322,2.815,0.002,3,3,1,0,
				      "Xi(c)(2815)+","\\Xi_{c}^{+}(2815)");
  s_kftable[102154]=new Particle_Info(102154,5.91,0.002,0,3,1,0,
				      "Lambda(b)(5910)","\\Lambda_{b}(5910)");
  // #########################################################################
  // Obsolete multiple heavy baryons #########################################
  // they will not be produced in our code (we have no heavy di-quarks
  // that we can produce yet)
  // #########################################################################
  s_kftable[4412]=new Particle_Info(4412,3.59798,0.,3,1,1,0,
				    "Xi(cc)+","\\Xi(cc)^{+}");
  s_kftable[4414]=new Particle_Info(4414,3.65648,0.,3,3,1,0,
				    "Xi(cc)*+","\\Xi(cc)*^{+}");
  s_kftable[4422]=new Particle_Info(4422,3.59798,0.,6,1,1,0,
				    "Xi(cc)++","\\Xi(cc)^{++}");
  s_kftable[4424]=new Particle_Info(4424,3.65648,0.,6,3,1,0,
				    "Xi(cc)*++","\\Xi(cc)*^{++}");
  s_kftable[4432]=new Particle_Info(4432,3.78663,0.,3,1,1,0,
				    "Omega(cc)+","\\Omega(cc)^{+}");
  s_kftable[4434]=new Particle_Info(4434,3.82466,0.,3,3,1,0,
				    "Omega(cc)*+","\\Omega(cc)*^{+}");
  s_kftable[4444]=new Particle_Info(4444,4.91594,0.,6,3,1,0,
				    "Omega(ccc)*++","\\Omega(ccc)*^{++}");
  s_kftable[5142]=new Particle_Info(5142,7.00575,0.,0,1,1,0,
				    "Xi(bc)","\\Xi(bc)");
  s_kftable[5242]=new Particle_Info(5242,7.00575,0.,3,1,1,0,
				    "Xi(bc)+","\\Xi(bc)^{+}");
  s_kftable[5342]=new Particle_Info(5342,7.19099,0.,0,1,1,0,
				    "Omega(bc)","\\Omega(bc)");
  s_kftable[5412]=new Particle_Info(5412,7.03724,0.,0,1,1,0,
				    "Xi(bc)'","\\Xi(bc)'");
  s_kftable[5414]=new Particle_Info(5414,7.0485,0.,0,3,1,0,
				    "Xi(bc)*","\\Xi(bc)*");
  s_kftable[5422]=new Particle_Info(5422,7.03724,0.,3,1,1,0,
				    "Xi(bc)'+","\\Xi(bc)'^{+}");
  s_kftable[5424]=new Particle_Info(5424,7.0485,0.,3,3,1,0,
				    "Xi(bc)*+","\\Xi(bc)*^{+}");
  s_kftable[5432]=new Particle_Info(5432,7.21101,0.,0,1,1,0,
				    "Omega(bc)'","\\Omega(bc)'");
  s_kftable[5434]=new Particle_Info(5434,7.219,0.,0,3,1,0,
				    "Omega(bc)*","\\Omega(bc)*");
  s_kftable[5442]=new Particle_Info(5442,8.30945,0.,3,1,1,0,
				    "Omega(bcc)+","\\Omega(bcc)^{+}");
  s_kftable[5444]=new Particle_Info(5444,8.31325,0.,3,3,1,0,
				    "Omega(bcc)*+","\\Omega(bcc)*^{+}");
  s_kftable[5512]=new Particle_Info(5512,10.42272,0.,-3,1,1,0,
				    "Xi(bb)-","\\Xi(bb)^{-}");
  s_kftable[5514]=new Particle_Info(5514,10.44144,0.,-3,3,1,0,
				    "Xi(bb)*-","\\Xi(bb)*^{-}");
  s_kftable[5522]=new Particle_Info(5522,10.42272,0.,0,1,1,0,
				    "Xi(bb)","\\Xi(bb)");
  s_kftable[5524]=new Particle_Info(5524,10.44144,0.,0,3,1,0,
				    "Xi(bb)*","\\Xi(bb)*");
  s_kftable[5532]=new Particle_Info(5532,10.60209,0.,-3,1,1,0,
				    "Omega(bb)-","\\Omega(bb)^{-}");
  s_kftable[5534]=new Particle_Info(5534,10.61426,0.,-3,3,1,0,
				    "Omega(bb)*-","\\Omega(bb)*^{-}");
  s_kftable[5542]=new Particle_Info(5542,11.70767,0.,0,1,1,0,
				    "Omega(bbc)","\\Omega(bbc)");
  s_kftable[5544]=new Particle_Info(5544,11.71147,0.,0,3,1,0,
				    "Omega(bbc)*","\\Omega(bbc)*");
  s_kftable[5554]=new Particle_Info(5554,15.11061,0.,-3,3,1,0,
				    "Omega(bbb)*-","\\Omega(bbb)*^{-}");

  
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
  s_kftable[117]=new Particle_Info(117,1.691,0.16,0,6,1,0,
				   "rho(3)(1690)","rho_{3}(1690)");
  s_kftable[217]=new Particle_Info(217,1.691,0.16,3,6,1,0,
				   "rho(3)(1690)+","rho_{3}^{+}(1690)");
  s_kftable[227]=new Particle_Info(227,1.667,0.168,0,6,1,0,
				   "omega(3)(1670)","omega_{3}(1670)");
  s_kftable[317]=new Particle_Info(317,1.776,0.159,0,6,1,0,
				   "K(3)*(1780)","K_{3}*(1780)");
  s_kftable[327]=new Particle_Info(327,1.776,0.159,3,6,1,0,
				   "K(3)*(1780)+","K_{3}*^{+}(1780)");
  s_kftable[337]=new Particle_Info(337,1.854,0.087,0,6,1,0,
				   "phi(3)(1850)","phi_{3}(1850)");
  s_kftable[557]=new Particle_Info(557,10.1599,0.0,0,6,1,0,
				   "Upsilon(3)(1D)","Upsilon_{3}(1D)");
  // Tensors 4     ############################################################
  // heavy ones missing.
  s_kftable[119]=new Particle_Info(119,2.014,0.361,0,8,1,0,
				   "a(4)(2040)","a_{4}(2040)");
  s_kftable[219]=new Particle_Info(219,2.014,0.361,3,8,1,0,
				   "a(4)(2040)+","a_{4}^{+}(2040)");
  s_kftable[229]=new Particle_Info(229,2.044,0.208,0,8,1,0,
				   "f(4)(2050)","f_{4}(2050)");
  s_kftable[319]=new Particle_Info(319,2.045,0.198,0,8,1,0,
				   "K(4)*(2045)","K_{4}*(2045)");
  s_kftable[329]=new Particle_Info(329,2.045,0.198,3,8,1,0,
				   "K(4)*(2045)+","K_{4}*^{+}(2045)");
  // Tensors 2     ############################################################
  // heavy ones missing.
  s_kftable[10115]=new Particle_Info(10115,1.67,0.258,0,4,1,0,
				     "pi(2)(1670)","pi_{2}(1670)");
  s_kftable[10215]=new Particle_Info(10215,1.67,0.258,3,4,1,0,
				     "pi(2)(1670)+","pi_{2}^{+}(1670)");
  s_kftable[10225]=new Particle_Info(10225,1.617,0.181,0,4,1,0,
				     "eta(2)(1645)","eta_{2}(1645)");
  s_kftable[10315]=new Particle_Info(10315,1.773,0.186,0,4,1,0,
				     "K(2)(1770)","K_{2}(1770)");
  s_kftable[10325]=new Particle_Info(10325,1.773,0.186,3,4,1,0,
				     "K(2)(1770)+","K_{2}^{+}(1770)");
  s_kftable[10335]=new Particle_Info(10335,1.842,0.225,0,4,1,0,
				     "eta(2)(1870)","eta_{2}(1870)");
  s_kftable[10555]=new Particle_Info(10555,10.157,0.0,0,4,1,0,
				     "eta(b2)(1D)","eta_{b2}(1D)");
  // Tensors 2     ############################################################
  // lots missing.
  s_kftable[20315]=new Particle_Info(20315,1.816,0.276,0,4,1,0,
				     "K(2)(1820)","K_{2}(1820)");
  s_kftable[20325]=new Particle_Info(20325,1.816,0.276,3,4,1,0,
				     "K(2)(1820)+","K_{2}^{+}(1820)");
  s_kftable[20555]=new Particle_Info(20555,10.1562,0.0,0,4,1,0,
				     "Upsilon(2)(1D)","Upsilon_{2}(1D)");
  s_kftable[30411]=new Particle_Info(30411,2.58,0.0,3,0,1,0,
				     "D(2S)+","D(2S)^{+}");
  s_kftable[30421]=new Particle_Info(30421,2.58,0.0,0,0,1,0,
				     "D(2S)","D(2S)");
  // Vectors 2     ############################################################
  // some heavy ones missing - we may have to include this because of the
  // psi(3770)
  s_kftable[30113]=new Particle_Info(30113,1.7,0.24,0,2,1,0,
				     "rho(1700)","rho(1700)");
  s_kftable[30213]=new Particle_Info(30213,1.7,0.24,3,2,1,0,
				     "rho(1700)+","rho^{+}(1700)");
  s_kftable[30223]=new Particle_Info(30223,1.670,0.31,0,2,1,0,
				     "omega(1650)","omega(1650)");
  s_kftable[30313]=new Particle_Info(30313,1.717,0.32,0,2,1,0,
				     "K*(1680)","K*(1680)");
  s_kftable[30323]=new Particle_Info(30323,1.717,0.32,3,2,1,0,
				     "K*(1680)+","K*^{+}(1680)");
  s_kftable[30333]=new Particle_Info(30333,1.900,0.32,0,2,1,0,
				     "f1(1900)_fict","f1(1900)_fict");
  s_kftable[30413]=new Particle_Info(30413,2.64,0.0,3,2,1,0,
				     "D*(2S)+","D*^{+}(2S)");
  s_kftable[30423]=new Particle_Info(30423,2.64,0.0,0,2,1,0,
				     "D*(2S)","D*(2S)");
  s_kftable[30443]=new Particle_Info(30443,3.7699,0.0236,0,2,1,0,
				     "psi(3770)","psi(3770)");
  s_kftable[30553]=new Particle_Info(30553,10.161,0.0,0,2,1,0,
				     "Upsilon(1)(1D)","Upsilon_{1}(1D)");
  // Pseudoscalars ############################################################
  // heavy ones missing.
  s_kftable[100111]=new Particle_Info(100111,1.3,0.400,0,0,1,0,
				      "pi(1300)","pi(1300)");
  s_kftable[100211]=new Particle_Info(100211,1.3,0.400,3,0,1,0,
				      "pi(1300)+","pi^{+}(1300)");
  s_kftable[100221]=new Particle_Info(100221,1.297,0.053,0,0,1,0,
				      "eta(1295)","eta(1295)");
  s_kftable[100311]=new Particle_Info(100311,1.46,0.26,0,0,1,0,
				      "K(1460)","K(1460)");
  s_kftable[100321]=new Particle_Info(100321,1.46,0.26,3,0,1,0,
				      "K(1460)+","K^{+}(1460)");
  s_kftable[100331]=new Particle_Info(100331,1.476,0.08,0,0,1,0,
				      "eta(1475)","eta(1475)");
  s_kftable[100441]=new Particle_Info(100441,3.638,0.014,0,0,1,0,
				      "eta(c)(2S)","eta_{c}(2S)");
  s_kftable[100551]=new Particle_Info(100551,9.997,0.0,0,0,1,0,
				      "eta(b)(2S)","eta_{b}(2S)");
  // Vectors       ###########################################################
  // heavy ones missing.
  // the rho's may be important for tau/D/B decays
  s_kftable[100113]=new Particle_Info(100113,1.465,0.31,0,2,1,0,
				      "rho(1450)","rho(1450)");
  s_kftable[100213]=new Particle_Info(100213,1.465,0.31,3,2,1,0,
				      "rho(1450)+","rho^{+}(1450)");
  s_kftable[100223]=new Particle_Info(100223,1.419,0.17,0,2,1,0,
				      "omega(1420)","omega(1420)");
  s_kftable[100313]=new Particle_Info(100313,1.414,0.232,0,2,1,0,
				      "K*(1410)","K*(1410)");
  s_kftable[100323]=new Particle_Info(100323,1.414,0.232,3,2,1,0,
				      "K*(1410)+","K*^{+}(1410)");
  s_kftable[100333]=new Particle_Info(100333,1.68,0.15,0,2,1,0,
				      "phi(1680)","phi(1680)");
  s_kftable[100443]=new Particle_Info(100443,3.686,0.000277,0,2,1,0,
				      "psi(2S)","psi(2S)");
  s_kftable[100553]=new Particle_Info(100553,10.0233,4.4e-05,0,2,1,0,
				      "Upsilon(2S)","Upsilon(2S)");
  // Tensors 2     ##########################################################
  // heavy ones missing.
  s_kftable[100335]=new Particle_Info(100335,2.01,0.2,0,4,1,0,
				      "f(2)(2010)","f_{2}(2010)");
  s_kftable[100445]=new Particle_Info(100445,3.929,0.029,0,4,1,0,
				      "chi(c2)(2P)","chi_{c2}(2P)");
  // More states without full multiplets ####################################
  // light ones missing - we know A LOT about heavy-heavy states ....
  s_kftable[100555]=new Particle_Info(100555,10.2685,0.001,0,4,1,0,
				      "chi(b2)(2P)","chi_{b2}(2P)");
  s_kftable[100557]=new Particle_Info(100557,10.4443,0.0,0,6,1,0,
				      "Upsilon(3)(2D)","Upsilon_{3}(2D)");
  s_kftable[110551]=new Particle_Info(110551,10.2321,0.001,0,0,1,0,
				      "chi(b0)(2P)","chi_{b0}(2P)");
  s_kftable[110553]=new Particle_Info(110553,10.255,0.0,0,2,1,0,
				      "h(b)(2P)","h_{b}(2P)");
  s_kftable[110555]=new Particle_Info(110555,10.441,0.0,0,4,1,0,
				      "eta(b2)(2D)","eta_{b2}(2D)");
  s_kftable[120553]=new Particle_Info(120553,10.2552,0.001,0,2,1,0,
				      "chi(b1)(2P)","chi_{b1}(2P)");
  s_kftable[120555]=new Particle_Info(120555,10.4406,0.0,0,4,1,0,
				      "Upsilon(2)(2D)","Upsilon_{2}(2D)");
  s_kftable[130553]=new Particle_Info(130553,10.4349,0.0,0,2,1,0,
				      "Upsilon(1)(2D)","Upsilon_{1}(2D)");
  s_kftable[200551]=new Particle_Info(200551,10.335,0.0,0,0,1,0,
				      "eta(b)(3S)","eta_{b}(3S)");
  s_kftable[200553]=new Particle_Info(200553,10.3553,2.63e-05,0,2,1,0,
				      "Upsilon(3S)","Upsilon(3S)");
  s_kftable[200555]=new Particle_Info(200555,10.5264,0.0,0,4,1,0,
				      "chi(b2)(3P)","chi_{b2}(3P)");
  s_kftable[210551]=new Particle_Info(210551,10.5007,0.0,0,0,1,0,
				      "chi(b0)(3P)","chi_{b0}(3P)");
  s_kftable[210553]=new Particle_Info(210553,10.516,0.0,0,2,1,0,
				      "h(b)(3P)","h_{b}(3P)");
  s_kftable[220553]=new Particle_Info(220553,10.516,0.0,0,2,1,0,
				      "chi(b1)(3P)","chi_{b1}(3P)");
  s_kftable[300553]=new Particle_Info(300553,10.58,0.01,0,2,1,0,
				      "Upsilon(4S)","Upsilon(4S)");
  // a0(980) and friends #####################################################
  // These are the "funny" state with some iso number.  We will have to
  // figure out what to do with them.
  // #########################################################################
  s_kftable[9000111]=new Particle_Info(9000111,0.996,0.075,0,0,1,0,
				       "a(0)(980)","a_{0}(980)");
  s_kftable[9000211]=new Particle_Info(9000211,0.996,0.075,3,0,1,0,
				       "a(0)(980)+","a_{0}^{+}(980)");
  s_kftable[9010221]=new Particle_Info(9010221,0.98,0.070,0,0,1,0,
				       "f(0)(980)","f_{0}(980)");
  s_kftable[9000221]=new Particle_Info(9000221,0.600,0.600,0,0,1,0,
				       "f(0)(600)","f_{0}(600)");
  s_kftable[9000311]=new Particle_Info(9000311,0.841,0.618,0,0,1,0,
				       "K(0)*(800)","K_{0}*(800)");
  s_kftable[9000321]=new Particle_Info(9000321,0.841,0.618,3,0,1,0,
				       "K(0)*(800)+","K_{0}*^{+}(800)");
  s_kftable[9000113]=new Particle_Info(9000113,1.376,0.3,0,2,1,0,
				       "pi(1)(1400)","pi_{1}(1400)");
  s_kftable[9000213]=new Particle_Info(9000213,1.376,0.3,3,2,1,0,
				       "pi(1)(1400)+","pi_{1}^{+}(1400)");
  s_kftable[9000223]=new Particle_Info(9000223,1.518,0.073,0,2,1,0,
				       "f(1)(1510)","f_{1}(1510)");
  s_kftable[9000313]=new Particle_Info(9000313,1.65,0.15,0,2,1,0,
				       "K(1)(1650)","K_{1}(1650)");
  s_kftable[9000323]=new Particle_Info(9000323,1.65,0.15,3,2,1,0,
				       "K(1)(1650)+","K_{1}^{+}(1650)");
  s_kftable[9000443]=new Particle_Info(9000443,4.04,0.052,0,2,1,0,
				       "psi(4040)","psi(4040)");
  s_kftable[9000553]=new Particle_Info(9000553,10.865,0.11,0,2,1,0,
				       "Upsilon(10860)","Upsilon(10860)");
  s_kftable[9000115]=new Particle_Info(9000115,1.732,0.194,0,4,1,0,
				       "a(2)(1700)","a_{2}(1700)");
  s_kftable[9000215]=new Particle_Info(9000215,1.732,0.194,3,4,1,0,
				       "a(2)(1700)+","a_{2}^{+}(1700)");
  s_kftable[9000225]=new Particle_Info(9000225,1.43,0.02,0,4,1,0,
				       "f(2)(1430)","f_{2}(1430)");
  s_kftable[9000315]=new Particle_Info(9000315,1.58,0.11,0,4,1,0,
				       "K(2)(1580)","K_{2}(1580)");
  s_kftable[9000325]=new Particle_Info(9000325,1.58,0.11,3,4,1,0,
				       "K(2)(1580)+","K_{2}^{+}(1580)");
  s_kftable[9000117]=new Particle_Info(9000117,1.982,0.188,0,6,1,0,
				       "rho(3)(1990)","rho_{3}(1990)");
  s_kftable[9000217]=new Particle_Info(9000217,1.982,0.188,3,6,1,0,
				       "rho(3)(1990)+","rho_{3}^{+}(1990)");
  s_kftable[9000229]=new Particle_Info(9000229,2.2311,0.023,0,8,1,0,
				       "f(J)(2220)","f(J)(2220)");
  s_kftable[9000319]=new Particle_Info(9000319,2.045,0.198,0,8,1,0,
				       "K(4)(2500)","K_{4}(2500)");
  s_kftable[9000329]=new Particle_Info(9000329,2.045,0.198,3,8,1,0,
				       "K(4)(2500)+","K_{4}^{+}(2500)");
  s_kftable[9010111]=new Particle_Info(9010111,1.812,0.207,0,0,1,0,
				       "pi(1800)","pi(1800)");
  s_kftable[9010211]=new Particle_Info(9010211,1.812,0.207,3,0,1,0,
				       "pi(1800)+","pi^{+}(1800)");
  s_kftable[9010311]=new Particle_Info(9010311,1.83,0.25,0,0,1,0,
				       "K(1830)","K(1830)");
  s_kftable[9010321]=new Particle_Info(9010321,1.83,0.25,3,0,1,0,
				       "K(1830)+","K^{+}(1830)");
  s_kftable[9010113]=new Particle_Info(9010113,1.653,0.225,0,2,1,0,
				       "pi(1)(1600)","pi_{1}(1600)");
  s_kftable[9010213]=new Particle_Info(9010213,1.653,0.225,3,2,1,0,
				       "pi(1)(1600)+","pi_{1}^{+}(1600)");
  s_kftable[9010223]=new Particle_Info(9010223,1.594,0.384,0,2,1,0,
				       "h(1)(1595)","h_{1}(1595)");
  s_kftable[9010553]=new Particle_Info(9010553,11.019,0.079,0,2,1,0,
				       "Upsilon(11020)","Upsilon(11020)");
  s_kftable[9010443]=new Particle_Info(9010443,4.159,0.078,0,2,1,0,
				       "psi(4160)","psi(4160)");
  s_kftable[9010115]=new Particle_Info(9010115,2.09,0.625,0,4,1,0,
				       "pi(2)(2100)","pi_{2}(2100)");
  s_kftable[9010215]=new Particle_Info(9010215,2.09,0.625,3,4,1,0,
				       "pi(2)(2100)+","pi_{2}^{+}(2100)");
  s_kftable[9010225]=new Particle_Info(9010225,1.546,0.126,0,4,1,0,
				       "f(2)(1565)","f_{2}(1565)");
  s_kftable[9010315]=new Particle_Info(9010315,1.973,0.373,0,4,1,0,
				       "K(2)*(1980)","K_{2}*(1980)");
  s_kftable[9010325]=new Particle_Info(9010325,1.973,0.373,3,4,1,0,
				       "K(2)*(1980)+","K_{2}*^{+}(1980)");
  s_kftable[9010117]=new Particle_Info(9010117,2.25,0.2,0,6,1,0,
				       "rho(3)(2250)","rho_{3}(2250)");
  s_kftable[9010217]=new Particle_Info(9010217,2.25,0.2,3,6,1,0,
				       "rho(3)(2250)+","rho_{3}^{+}(2250)");
  s_kftable[9010317]=new Particle_Info(9010317,2.324,0.15,0,6,1,0,
				       "K(3)(2320)","K_{3}(2320)");
  s_kftable[9010327]=new Particle_Info(9010327,2.324,0.15,3,6,1,0,
				       "K(3)(2320)+","K_{3}^{+}(2320)");
  s_kftable[9010229]=new Particle_Info(9010229,2.332,0.26,0,8,1,0,
				       "f(4)(2300)","f_{4}(2300)");
  s_kftable[9020221]=new Particle_Info(9020221,1.4098,0.0511,0,0,1,0,
				       "eta(1405)","eta(1405)");
  s_kftable[9020311]=new Particle_Info(9020311,1.945,0.201,0,0,1,0,
				       "K(0)*(1950)","K_{0}*(1950)");
  s_kftable[9020321]=new Particle_Info(9020321,1.945,0.201,3,0,1,0,
				       "K(0)*(1950)+","K_{0}*^{+}(1950)");
  s_kftable[9020113]=new Particle_Info(9020113,1.647,0.254,0,2,1,0,
				       "a(1)(1640)","a_{1}(1640)");
  s_kftable[9020213]=new Particle_Info(9020213,1.647,0.254,3,2,1,0,
				       "a(1)(1600)+","a_{1}^{+}(1600)");
  s_kftable[9020443]=new Particle_Info(9020443,4.415,0.043,0,2,1,0,
				       "psi(4415)","psi(4415)");
  s_kftable[9020225]=new Particle_Info(9020225,1.638,0.099,0,4,1,0,
				       "f(2)(1640)","f_{2}(1640)");
  s_kftable[9020315]=new Particle_Info(9020315,2.247,0.18,0,4,1,0,
				       "K(2)(2250)","K_{2}(2250)");
  s_kftable[9020325]=new Particle_Info(9020325,2.247,0.18,3,4,1,0,
				       "K(2)(2250)+","K_{2}^{+}(2250)");
  s_kftable[9030221]=new Particle_Info(9030221,1.5,0.112,0,0,1,0,
				       "f(0)(1500)","f_{0}(1500)");
  s_kftable[9030113]=new Particle_Info(9030113,1.9,0.16,0,2,1,0,
				       "rho(1900)","rho(1900)");
  s_kftable[9030213]=new Particle_Info(9030213,1.9,0.16,3,2,1,0,
				       "rho(1900)+","rho^{+}(1900)");
  s_kftable[9030225]=new Particle_Info(9030225,1.815,0.197,0,4,1,0,
				       "f(2)(1810)","f_{2}(1810)");
  s_kftable[9040221]=new Particle_Info(9040221,1.76,0.06,0,0,1,0,
				       "eta(1760)","eta(1760)");
  s_kftable[9040113]=new Particle_Info(9040113,2.149,0.363,0,2,1,0,
				       "rho(2150)","rho(2150)");
  s_kftable[9040213]=new Particle_Info(9040213,2.149,0.363,3,2,1,0,
				       "rho(2150)+","rho^{+}(2150)");
  s_kftable[9040225]=new Particle_Info(9040225,1.915,0.163,0,4,1,0,
				       "f(2)(1910)","f_{2}(1910)");
  s_kftable[9050221]=new Particle_Info(9050221,1.992,0.442,0,0,1,0,
				       "f(0)(2020)","f_{0}(2020)");
  s_kftable[9050225]=new Particle_Info(9050225,1.944,0.472,0,4,1,0,
				       "f(2)(1950)","f_{2}(1950)");
  s_kftable[9060221]=new Particle_Info(9060221,2.103,0.206,0,0,1,0,
				       "f(0)(2100)","f_{0}(2100)");
  s_kftable[9060225]=new Particle_Info(9060225,2.011,0.202,0,4,1,0,
				       "f(2)(2010)","f_{2}(2010)");
  s_kftable[9070221]=new Particle_Info(9070221,2.189,0.238,0,0,1,0,
				       "f(0)(2200)","f_{0}(2200)");
  s_kftable[9070225]=new Particle_Info(9070225,2.15,0.2,0,4,1,0,
				       "f(2)(2150)","f_{2}(2150)");
  s_kftable[9080221]=new Particle_Info(9080221,2.220,0.15,0,0,1,0,
				       "eta(2225)","eta(2225)");
  s_kftable[9080225]=new Particle_Info(9080225,2.297,0.15,0,4,1,0,
				       "f(2)(2300)","f_{2}(2300)");
  s_kftable[9090225]=new Particle_Info(9090225,2.34,0.32,0,4,1,0,
				       "f(2)(2340)","f_{2}(2340)");

  // ##########################################################################
  // ##########################################################################
  // Former members of the Sherpa team - made immortal as particles here ######
  // ##########################################################################
  // ##########################################################################
  s_kftable[5505]=new Particle_Info(5505,1000000.0,1000,0,0,0,0,
				    "ralf","ralf");
  s_kftable[5506]=new Particle_Info(5506,1000000.0,1000,0,0,0,0,
				    "ande","ande");
  s_kftable[5507]=new Particle_Info(5507,1000000.0,1000,0,0,0,0,
				    "thomas","thomas");
  s_kftable[5508]=new Particle_Info(5508,1000000.0,1000,0,0,0,0,
				    "tanju","tanju");
  s_kftable[5508]=new Particle_Info(5509,1000000.0,1000,0,0,0,0,
				    "jennifer","jennifer");
  s_kftable[5508]=new Particle_Info(6505,1000000.0,1000,0,0,0,0,
				    "hendrik","hendrik");
  s_kftable[5508]=new Particle_Info(6506,1000000.0,1000,0,0,0,0,
				    "jan","jan");
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################
  // ##########################################################################


  // set self-anti property
  s_kftable[kf_pi]->m_majorana=-1;
  s_kftable[kf_eta]->m_majorana=-1;
  s_kftable[kf_eta_prime_958]->m_majorana=-1;
  s_kftable[kf_K_L]->m_majorana=-1;
  s_kftable[kf_K_S]->m_majorana=-1;
  s_kftable[kf_eta_c_1S]->m_majorana=-1;
  s_kftable[kf_eta_b]->m_majorana=-1;
  s_kftable[kf_rho_770]->m_majorana=-1;
  s_kftable[kf_omega_782]->m_majorana=-1;
  s_kftable[kf_phi_1020]->m_majorana=-1;
  s_kftable[kf_J_psi_1S]->m_majorana=-1;
  s_kftable[kf_Upsilon_1S]->m_majorana=-1;
  s_kftable[kf_a_2_1320]->m_majorana=-1;
  s_kftable[kf_f_2_1270]->m_majorana=-1;
  s_kftable[kf_f_2_prime_1525]->m_majorana=-1;
  s_kftable[kf_chi_c2_1P]->m_majorana=-1;
  s_kftable[kf_chi_b2_1P]->m_majorana=-1;
  s_kftable[kf_rho_3_1690]->m_majorana=-1;
  s_kftable[kf_omega_3_1670]->m_majorana=-1;
  s_kftable[kf_phi_3_1850]->m_majorana=-1;
  s_kftable[kf_a_4_2040]->m_majorana=-1;
  s_kftable[kf_f_4_2050]->m_majorana=-1;
  //s_kftable[kf_f_J_2220]->m_majorana=-1;
  s_kftable[kf_a_0_980]->m_majorana=-1;
  s_kftable[kf_f_0_980]->m_majorana=-1;
  s_kftable[kf_f_0_1370]->m_majorana=-1;
  s_kftable[kf_chi_c0_1P]->m_majorana=-1;
  s_kftable[kf_chi_b0_1P]->m_majorana=-1;
  s_kftable[kf_b_1_1235]->m_majorana=-1;
  s_kftable[kf_h_1_1170]->m_majorana=-1;
  s_kftable[kf_h_1_1380]->m_majorana=-1;
  s_kftable[kf_h_c1]->m_majorana=-1;
  s_kftable[kf_h_b1]->m_majorana=-1;
  s_kftable[kf_pi_2_1670]->m_majorana=-1;
  s_kftable[kf_eta_2_1645]->m_majorana=-1;
  s_kftable[kf_eta_2_1870]->m_majorana=-1;
  s_kftable[kf_a_1_1260]->m_majorana=-1;
  s_kftable[kf_f_1_1285]->m_majorana=-1;
  s_kftable[kf_f_1_1420]->m_majorana=-1;
  s_kftable[kf_chi_c1_1P]->m_majorana=-1;
  s_kftable[kf_chi_b1_1P]->m_majorana=-1;
  s_kftable[kf_rho_1700]->m_majorana=-1;
  s_kftable[kf_omega_1600]->m_majorana=-1;
  s_kftable[kf_f1_1900_fict]->m_majorana=-1;
  s_kftable[kf_psi_3770]->m_majorana=-1;
  s_kftable[kf_pi_1300]->m_majorana=-1;
  s_kftable[kf_eta_1295]->m_majorana=-1;
  s_kftable[kf_eta_1475]->m_majorana=-1;
  s_kftable[kf_rho_1450]->m_majorana=-1;
  s_kftable[kf_omega_1420]->m_majorana=-1;
  s_kftable[kf_phi_1680]->m_majorana=-1;
  s_kftable[kf_psi_2S]->m_majorana=-1;
  s_kftable[kf_Upsilon_2S]->m_majorana=-1;
  s_kftable[kf_f_2_2010]->m_majorana=-1;
  s_kftable[kf_chi_b2_2P]->m_majorana=-1;
  s_kftable[kf_chi_b0_2P]->m_majorana=-1;
  s_kftable[kf_chi_b1_2P]->m_majorana=-1;
  s_kftable[kf_Upsilon_3S]->m_majorana=-1;
  s_kftable[kf_Upsilon_4S]->m_majorana=-1;
  s_kftable[kf_a_0_1450]->m_majorana=-1;
  s_kftable[kf_f_0_600]->m_majorana=-1;
  s_kftable[kf_psi_4040]->m_majorana=-1;
  s_kftable[kf_Upsilon_10860]->m_majorana=-1;
  s_kftable[kf_f_0_1710]->m_majorana=-1;
  s_kftable[kf_psi_4160]->m_majorana=-1;
  s_kftable[kf_Upsilon_11020]->m_majorana=-1;
  s_kftable[kf_f_0_1500]->m_majorana=-1;
  s_kftable[kf_psi_4415]->m_majorana=-1;
  s_kftable[kf_f_J_1710]->m_majorana=-1;
  s_kftable[kf_f_2_2300]->m_majorana=-1;
  s_kftable[kf_f_2_2340]->m_majorana=-1;
}

void Hadron_Init::OverWriteProperties(Data_Reader& dr)
{
  std::map<int,double> cdm, cdw;
  std::map<int,int> cia, cis, cim;
  std::vector<std::vector<double> > helpdvv;
  if (dr.MatrixFromFile(helpdvv,"MASS"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdm[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"WIDTH"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cdw[int(helpdvv[i][0])]=helpdvv[i][1];
  if (dr.MatrixFromFile(helpdvv,"ACTIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cia[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"STABLE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cis[int(helpdvv[i][0])]=int(helpdvv[i][1]);
  if (dr.MatrixFromFile(helpdvv,"MASSIVE"))
    for (size_t i(0);i<helpdvv.size();++i)
      if (helpdvv[i].size()==2) cim[int(helpdvv[i][0])]=int(helpdvv[i][1]);

  //set masses
  std::map<int,double>::const_iterator dit=cdm.begin();
  for (;dit!=cdm.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_mass = dit->second;
      msg_Tracking()<<" set mass of "<<Flavour(dit->first)
		    <<" to "<<dit->second<<" GeV"<<std::endl; 
    }
  }
  //set widths
  dit=cdw.begin();
  for (;dit!=cdw.end();dit++) {
    if (s_kftable.find(dit->first)!=s_kftable.end()) {
      s_kftable[dit->first]->m_width = dit->second;
      msg_Tracking()<<" set width of "<<Flavour(dit->first)<<" to "<<dit->second<<" GeV"<<std::endl; 
    }
  }
  //set (in)active
  std::map<int,int>::const_iterator iit=cia.begin();
  for (;iit!=cia.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_on = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" inactive "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" active "<<std::endl; 
      }
    }
  }
  //set (un)stable
  iit=cis.begin();
  for (;iit!=cis.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_stable = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" unstable "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" stable "<<std::endl; 
      }
    }
  }
  //set massive/massless
  iit=cim.begin();
  for (;iit!=cim.end();iit++) {
    if (s_kftable.find(iit->first)!=s_kftable.end()) {
      s_kftable[iit->first]->m_massive = iit->second;
      if (iit->second==0) {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" massless "<<std::endl; 
      }
      else {
	msg_Tracking()<<" set flavour "<<Flavour(iit->first)<<" massive "<<std::endl; 
      }
    }
  }

}
