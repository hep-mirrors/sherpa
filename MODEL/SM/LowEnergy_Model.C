#include "MODEL/SM/LowEnergy_Model.H"
#include "MODEL/SM/ThreePion_FittedCurrent.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Form_Factor.H"
#include "METOOLS/FormFactors/Line_Shapes.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(MODEL::LowEnergy_Model,"LowEnergy",MODEL::Model_Base,MODEL::Model_Arguments);

Model_Base *ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::LowEnergy_Model>::
operator()(const Model_Arguments &args) const
{
  return new LowEnergy_Model();
}

void ATOOLS::Getter<MODEL::Model_Base,MODEL::Model_Arguments,MODEL::LowEnergy_Model>::
PrintInfo(ostream &str,const size_t width) const
{
  str<<"The LowEnergy Model\n";
  str<<setw(width+4)<<" "<<"{\n"
     <<setw(width+7)<<" "<<"# possible parameters in yaml configuration [usage: \"keyword: value\"]\n"
     <<setw(width+7)<<" "<<"- 1/ALPHAQED(0) (alpha QED Thompson limit)\n"
     <<setw(width+7)<<" "<<"- SIN2THETAW (weak mixing angle)\n"
     <<setw(width+4)<<" "<<"}";
}

LowEnergy_Model::LowEnergy_Model() : Model_Base(true)
{
  msg_Out()<<METHOD<<" starts initialising.\n";
  m_name="LowEnergy";
  ParticleInit();
  InitQEDConstants();
  RegisterDefaults();
}

bool LowEnergy_Model::ModelInit()
{
  Settings& s = Settings::GetMainSettings();
  m_alpha     = 1./s["1/ALPHAQED(0)"].Get<double>();
  m_sinthetaW = sqrt(s["SIN2THETAW"].Get<double>());
  // Set THREE_PION_CONTACT: true to go back to the contact pi0 pi+ pi-
  // vertex instead of the gamma -> rho pi chain.
  m_threepion_contact =
    s["THREE_PION_CONTACT"].SetDefault(false).Get<bool>();
  msg_Out()<<METHOD<<": 1/alpha = "<<(1./m_alpha)<<"\n";
  METOOLS::LineShapes = new METOOLS::Line_Shapes();
  METOOLS::LineShapes->Init();
  return true;
}

void LowEnergy_Model::ParticleInit() {
  //          kf_code,        mass,    radius, width,  3*charge,2*spin,on,stable,idname,texname
  AddParticle(kf_n,           0.939566,0.8783, 7.424e-28, 0,    1,     true,true,  "n","n");
  AddParticle(kf_p_plus,      0.938272,0.8783, 0.0,       3,    1,     true,true,  "P+","P^{+}");
  AddParticle(kf_Sigma_minus, 1.19745, 0.8783, 4.45e-15, -3,    1,     true,false, "Sigma-","\\Sigma^{-}");
  AddParticle(kf_Sigma,       1.19264, 0.8783, 8.9e-06,   0,    1,     true,false, "Sigma","\\Sigma");
  AddParticle(kf_Sigma_plus,  1.18937, 0.8783, 8.24e-15,  3,    1,     true,false, "Sigma+","\\Sigma^{+}");
  AddParticle(kf_Lambda,      1.11568, 0.8783, 2.501e-15, 0,    1,     true,false, "Lambda","\\Lambda");
  AddParticle(kf_Xi_minus,    1.32132, 0.8783, 4.02e-15, -3,    1,     true,false, "Xi-","\\Xi^{-}");
  AddParticle(kf_Xi,          1.3149,  0.8783, 2.27e-15,  0,    1,     true,false, "Xi","\\Xi");
  AddParticle(kf_ISR_photon,  0.0,     0.0,    0.0,       0,    2,     true, true, "ISRPhoton","\\gamma");
  AddParticle(kf_FSR_photon,  0.0,     0.0,    0.0,       0,    2,     true, true, "FSRPhoton","\\gamma");

  m_Diracs.push_back(Flavour(kf_n));
  m_Diracs.push_back(Flavour(kf_p_plus));
  m_Diracs.push_back(Flavour(kf_Sigma_minus));
  m_Diracs.push_back(Flavour(kf_Sigma));
  m_Diracs.push_back(Flavour(kf_Sigma_plus));
  m_Diracs.push_back(Flavour(kf_Lambda));
  m_Diracs.push_back(Flavour(kf_Xi_minus));
  m_Diracs.push_back(Flavour(kf_Xi));

  //          kf_code,        mass,    radius, width,  3*charge,2*spin,"Majorana",on,stable,idname,texname
  AddParticle(kf_pi,           0.134976, 0.65, 7.8486e-09, 0,   0,     false,     true,false, "pi","pi");
  AddParticle(kf_pi_plus,      0.13957,  0.65, 2.5242e-17, 3,   0,     true,      true,true,  "pi+","pi^{+}");
  AddParticle(kf_eta,          0.5473,   0.65, 1.18e-06,   0,   0,     false,     true,false, "eta","eta");
  AddParticle(kf_K,            0.49767,  0.65, 1.e-16,     0,   0,     true,      true,false, "K","K");
  AddParticle(kf_K_L,          0.49767,  0.65, 1.273e-17,  0,   0,     false,     true,true,  "K(L)","K_{L}");
  AddParticle(kf_K_S,          0.49767,  0.65, 7.373e-15,  0,   0,     false,     true,false, "K(S)","K_{S}");
  AddParticle(kf_K_plus,       0.493677, 0.65, 5.314e-17,  3,   0,     true,      true,true,  "K+","K^{+}");
  AddParticle(kf_eta_prime_958,0.95778,  0.65, 0.000203,   0,   0,     false,     true,false, "eta'(958)","eta'(958)");
  AddParticle(kf_J_psi_1S,     3.0969,   0.0,  0.0000926,  0,   1,     false,     true,false, "J/psi","J/\\psi");
  
  
  m_PseudoScalars.push_back(kf_pi);
  m_PseudoScalars.push_back(kf_pi_plus);
  m_PseudoScalars.push_back(kf_eta);
  m_PseudoScalars.push_back(kf_K);
  m_PseudoScalars.push_back(kf_K_L);
  m_PseudoScalars.push_back(kf_K_S);
  m_PseudoScalars.push_back(kf_K_plus);
  m_PseudoScalars.push_back(kf_eta_prime_958);

  AddParticle(kf_f_0_600,      0.600,    0.65, 0.600,      0,   0,     false,     true,false, "f(0)(600)","f_{0}(600)");
  AddParticle(kf_rho_770,      0.77,     0.65, 0.1507,     0,   2,     false,     true,false, "rho(770)","rho(770)");
  AddParticle(kf_rho_770_plus, 0.77,     0.65, 0.1507,     3,   2,     true,      true,false, "rho(770)+","rho^{+}(770)");
  AddParticle(kf_omega_782,    0.78194,  0.65, 0.00841,    0,   2,     false,     true,false, "omega(782)","omega(782)");
  AddParticle(kf_rho_1450,     1.465,    0.65, 0.31,       0,   2,     false,     true,false, "rho(1450)","rho(1450)");
  AddParticle(100213,          1.465,    0.65, 0.31,       3,   2,     true,      true,false, "rho(1450)+","rho^{+}(1450)");
  AddParticle(kf_rho_1700,     1.7,      0.65, 0.24,       0,   2,     false,     true,false, "rho(1700)","rho(1700)");
  AddParticle(30213,           1.7,      0.65, 0.24,       3,   2,     true,      true,false, "rho(1700)+","rho^{+}(1700)");
  AddParticle(kf_rho_2150,      2.15,   0.65, 0.11327,    0,   2,     false,     true,false, "rho(2150)","rho(2150)");
  AddParticle(kf_rho_2150_plus, 2.15,   0.65, 0.11327,    3,   2,     true,      true,false, "rho(2150)+","rho^{+}(2150)");


}

void LowEnergy_Model::InitQEDConstants() {
  // This could/should become part of a data file, maybe in yaml format ... .
  (*p_constants)["Lambda2_2212_2212_22"]  =  0.71;
  (*p_constants)["Q_2212_2212_22"]        =  1.;
  (*p_constants)["Mu_2212_2212_22"]       =  2.792847;
  (*p_constants)["Lambda2_2112_2112_22"]  =  0.71;
  (*p_constants)["Q_2112_2112_22"]        =  0.;
  (*p_constants)["Mu_2112_2112_22"]       = -1.913;
  // Numbers below by taking the magnetic moments of the
  // hyperons from the PDG.
  // We should compare with cross sections listed in
  // https://arxiv.org/html/2412.07543v1
  (*p_constants)["Lambda2_3122_3122_22"]  =  0.71;
  (*p_constants)["Q_3122_3122_22"]        =  0.;
  (*p_constants)["Mu_3122_3122_22"]       = -0.613;
  (*p_constants)["Lambda2_3222_3222_22"]  =  1.04;
  (*p_constants)["Q_3222_3222_22"]        =  1.;
  (*p_constants)["Mu_3222_3222_22"]       =  2.458;
  (*p_constants)["Lambda2_3212_3212_22"]  =  1.04;
  (*p_constants)["Q_3212_3212_22"]        =  0.;
  (*p_constants)["Mu_3212_3212_22"]       =  1.61;
  (*p_constants)["Lambda2_3112_3112_22"]  =  1.04;
  (*p_constants)["Q_3112_3112_22"]        = -1.;
  (*p_constants)["Mu_3112_3112_22"]       = -1.160;
  (*p_constants)["Lambda2_3322_3322_22"]  =  1.04;
  (*p_constants)["Q_3322_3322_22"]        =  0.;
  (*p_constants)["Mu_3322_3322_22"]       = -1.250;
  (*p_constants)["Lambda2_3312_3312_22"]  =  1.04;
  (*p_constants)["Q_3312_3312_22"]        = -1.;
  (*p_constants)["Mu_3312_3312_22"]       = -0.6507;
  (*p_constants)["Lambda2_3332_3332_22"]  =  1.04;
  (*p_constants)["Q_3332_3332_22"]        = -1.;
  (*p_constants)["Mu_3332_3332_22"]       = -2.02;
  // The mass for the axial form factor can be obtained from
  // 10.1088/0954-3899/28/1/201
  // For a start we could assume: M_A^2 ~ 1.04 GeV^2
  // across the board - this needs to be implemented
}

void LowEnergy_Model::InitVertices() {
  InitQEDVertices();
  InitEWVertices();
}

void LowEnergy_Model::InitQEDVertices() {
  Kabbala g1("g_1",sqrt(4.*M_PI*m_alpha));
  Kabbala cpl=g1*Kabbala("i",Complex(0.,1.));
  Flavour flav;
  for (list<Flavour>::iterator flit=m_Diracs.begin();
       flit!=m_Diracs.end();flit++) {
    Flavour flav = *flit;
    // only create vertices for hadrons that are switched on.
    Kabbala Q("Q_{"+flav.TexName()+"}",flav.Charge());
    msg_Out()<<METHOD<<" for "<<flav<<": "<<flav.IntSpin()<<"\n";
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("FFV");
    m_v.back().FormFactor.push_back("Dirac_F1");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("FFVM");
    m_v.back().FormFactor.push_back("Dirac_F2");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;
  }
  for (list<Flavour>::iterator flit=m_PseudoScalars.begin();
       flit!=m_PseudoScalars.end();flit++) {
    Flavour flav = *flit;
    if (flav.IntCharge()==0) continue;
    // only create vertices for hadrons that are switched on.
    Kabbala Q("Q_{"+flav.TexName()+"}",flav.Charge());
    msg_Out()<<METHOD<<" for "<<flav<<": "<<flav.IntSpin()<<"\n";
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(flav.Bar());
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;

    //Seagulls  
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVSS");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(2.0*cpl*g1);
    m_v.back().order[1] = 2;
    
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_FSR_photon));
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(flav.Bar());
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;

    // seagull with one FSR photon
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_FSR_photon));
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VVSS");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(2.0*cpl*g1);
    m_v.back().order[1] = 2;
  }

  //////////////////////////////////////////////////////////////////////
  // Contact pi0 pi+ pi- vertex.  Superseded by the gamma -> rho pi chain
  // below, which carries the correct anomalous structure - the two would
  // double count, so only one of them may be active at a time.
  //////////////////////////////////////////////////////////////////////
  if (m_threepion_contact) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_pi));
    m_v.back().AddParticle(Flavour(kf_pi_plus));
    m_v.back().AddParticle(Flavour(kf_pi_plus).Bar());
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSS");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=2;
  }
  else {
    InitThreePionVertices(cpl);
  }
  InitTwoPhotonVertices();
  InitFourPionVertices(cpl);
  // add the ISR photon vertex, which is a special case of the photon vertex
  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(Flavour(kf_e));
  m_v.back().AddParticle(Flavour(kf_e).Bar());
  m_v.back().AddParticle(Flavour(kf_ISR_photon));
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("FFV");
  m_v.back().cpl.push_back(cpl);
  m_v.back().order[1]=1;
}

//////////////////////////////////////////////////////////////////////
//
// gamma* -> V -> rho pi -> pi pi pi.
//
// The anomalous gamma-rho-pi vertex carries the Wess-Zumino structure
//   L = eps_{mu,nu,al,be} A^mu rho^nu p_A^al p_rho^be
//////////////////////////////////////////////////////////////////////
void LowEnergy_Model::InitThreePionVertices(const Kabbala & cpl) {
  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  const Flavour photon(kf_photon);
  const Flavour pi0(kf_pi), pip(kf_pi_plus), pim(Flavour(kf_pi_plus).Bar());
  const Flavour rho0(kf_rho_770);
  const Flavour rhop(kf_rho_770_plus), rhom(Flavour(kf_rho_770_plus).Bar());

  ////////////////////////////////////////////////////////////////////
  // rho -> pi pi, from the rho width itself.  With the standard isospin
  // coupling L = g rho^mu.(pi x d_mu pi) the amplitude is g eps.(p1-p2),
  // whose polarisation sum is g^2 (m_rho^2 - 4 m_pi^2) = 4 g^2 p^2; folding
  // that into Gamma = p/(8 pi m^2) <|M|^2> and averaging over the three
  // polarisations gives
  //   Gamma_rho = g_rhopipi^2 p^3 / (6 pi m_rho^2),
  //   p         = sqrt(m_rho^2/4 - m_pi^2),
  // so the coupling is just the width read back.  The particle table here
  // (m = 0.77, Gamma = 0.1507) gives 6.04; the PDG rho gives 5.94, which
  // agrees to 0.1% with the KSRF relation m_rho/(sqrt(2) f_pi) = 5.93.
  //
  // SSV_LC evaluates cpl*(p1-p2)^mu, exactly the structure above, so the
  // coupling goes in as-is - the same convention under which cpl = i e on
  // the photon vertices reproduces scalar QED. (RhoPiPiCoupling, shared with
  // the fitted current's production form factors.)
  ////////////////////////////////////////////////////////////////////
  const double mrho(rho0.Mass());
  const double grpp(RhoPiPiCoupling(rho0,"g_rhopipi"));

 
  ////////////////////////////////////////////////////////////////////
  // g_gammarhopi from the CHIRAL ANOMALY, not from a fit.
  //
  // The Wess-Zumino-Witten low-energy theorem PREDICTS the gamma* -> 3 pi
  // amplitude at threshold, so this coupling is not free.  Writing the
  // isoscalar amplitude the way arXiv:1907.01556 eq. (2.13) does,
  //     alpha_A = F_3pi/3 * 1.066,   F_3pi = 1/(4 pi^2 F_pi^3),
  // and matching our two-vertex chain to it: the AVVP vertex supplies the
  // bare epsilon contraction (AVVP_LC.C), the internal rho line supplies
  // -i/(p^2-m^2+i m Gamma) = (i/m_rho^2) BW (V_C.C:243), the SSV decay
  // supplies (p_a-p_b)^mu (SSV_LC.C), and folding the SSV momentum
  // difference into the AVVP epsilon slot gives a factor 2 relative to
  // eps(e,p_+,p_-,p_0) - measured, exactly 2.  Hence
  //     2 g_gammarhopi g_rhopipi / m_rho^2 = alpha_A
  // =>  g_gammarhopi = alpha_A m_rho^2 / (2 g_rhopipi) = 0.5598 /GeV
  // with F_pi = 0.1307/sqrt(2), the same value InitTwoPhotonVertices uses.
  //
  // This replaces a hand-set 0.773 that had no derivation behind it (the
  // VMD estimate on the line below was computed and then not used).  The
  // predicted value is a factor 0.724 smaller in amplitude, 0.52 in cross
  // section.  Overriding it with g_gammarhopi: in the runcard is putting a
  // fudge factor on a theorem - if the rate is off, the physics is
  // elsewhere.
  //
  // const double ee(sqrt(4.*M_PI*m_alpha));
  // const double ggrp_vmd(ee*14.0/17.0);   // 0.25 /GeV, for reference
  //
  // (GammaRhoPiCoupling, shared with the fitted current's production form
  // factors.)
  const double ggrp(GammaRhoPiCoupling());

  msg_Out()<<METHOD<<": g_rhopipi = "<<grpp
	   <<", g_gammarhopi = "<<ggrp<<" /GeV\n";
  // e g_gammarhopi: the photon couples to the beams with i e
  const Kabbala cpl_grp(cpl*Kabbala("g_{\\gamma\\rho\\pi}",ggrp));
  const Kabbala cpl_rpp("g_{\\rho\\pi\\pi}",
			grpp*Complex(0.,1.));

  ////////////////////////////////////////////////////////////////////
  // THREE_PION_FORM_FACTOR: {Contact: true} replaces the whole
  // gamma->rho pi -> pi pi pi chain below with a single four-leg contact
  // current, gamma pi+ pi- pi0, following arXiv:1907.01556 as implemented in
  // tau-decays-in-sherpa.  The epsilon structure lives in VPPP_LC.C and the
  // dynamics - anomaly-normalised isoscalar times an all-plus rho sum over
  // the three two-pion sub-masses - in FFVMD's ee3pi mode.
  //
  // The rho->pi pi vertices are still built below, unconditionally: the rho
  // remains a real particle needed by other final states (four-pion, etc.).
  // Only the gamma->rho pi production vertices are replaced, which is what
  // removes the double route to pi+ pi- pi0.
  ////////////////////////////////////////////////////////////////////
  const bool contact3pi(s["Contact"].SetDefault(false).Get<bool>());
  // Current: Fitted - the fitted current's isoscalars on the production
  // vertices (FFVMD) and its rho line shapes on the rho -> pi pi vertices
  // (ThreePionFittedRho); the rho pi chain only
  const bool fitted3pi(METOOLS::Fitted_ThreePion_Current::Selected());
  if (fitted3pi && contact3pi)
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: Current: Fitted is built on the "
	  "rho pi chain; Contact: true is not supported with it.");
  if (fitted3pi)
    msg_Out()<<METHOD<<": THREE_PION_FORM_FACTOR: Current: Fitted, "
	     <<METOOLS::Fitted_ThreePion_Current::Instance().Describe()<<"\n";
  // RhoWidth: Fixed (Comix's fixed-width propagator) or PWave (the RhoPWave
  // form factor on the rho -> pi pi vertices turns it into the P-wave running
  // width of hep-ph/0512180 eq. (5))
  const std::string rhowidth(s["RhoWidth"].SetDefault("Fixed").Get<std::string>());
  if (rhowidth!="Fixed" && rhowidth!="PWave")
    THROW(fatal_error,"THREE_PION_FORM_FACTOR: RhoWidth must be Fixed or PWave, not '"
	  +rhowidth+"'.");
  const std::string rhopipiff(fitted3pi ? "ThreePionFittedRho" :
			      rhowidth=="PWave" ? "RhoPWave" : "");
  if (contact3pi) {
    msg_Out()<<METHOD<<": three-pion contact current (VPPP) active; "
	     <<"gamma-rho-pi production vertices disabled.\n";
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    m_v.back().AddParticle(pip);
    m_v.back().AddParticle(pim);
    m_v.back().AddParticle(pi0);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VPPP");
    m_v.back().FormFactor.push_back("VMD");
    // The anomaly normalisation sits in the form factor (its F(0) matches the
    // Wess-Zumino-Witten theorem's 9.7 GeV^-3 to 1.066, checked), so the
    // coupling is the bare electric charge: one photon leg, no free parameter.
    //
    // NB an earlier sqrt(2) here was WRONG and has been removed.  Measured
    // against BELLE_2024_I2775022 the required correction is NOT flat: ~2 in
    // sigma through the resonance region (0.9-1.4 GeV) but ~1 in the tail
    // (>1.7 GeV).  A vertex coupling is energy-independent and so lifts the
    // tail it should not; the deficit is resonant, so any correction belongs
    // on the isoscalar resonance weights, which are large at the peaks and
    // vanish at high q^2 - see the note in ConstructEE3PiFormFactor.
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;
  }

  // gamma rho pi, anomalous.  Charges have to add up to zero: all legs are
  // read as incoming.
  //
  // ISOSPIN FACTORS.  I = 0: the three channels add all-plus with equal
  // weight; the charged ones come out of the epsilon with opposite signs,
  // hence (+1, -1, +1).  rho0_isospin: weight of rho0 pi0 (formerly 0.5).
  //
  // The same statement, and the same measurement, appear in
  // tau-decays-in-sherpa (METOOLS-hadron-currents),
  // FF_0_Isoscalar3Pi::BuildChains(), which carries isospin factors
  // (-1,+1,-1) for its own channel ordering and notes the momentum epsilon
  // alone gives (-,+,-).  Without these the charged channels partially
  // cancel instead of adding, which distorts the Dalitz distribution and so
  // shows up as an ENERGY-DEPENDENT discrepancy - not something any overall
  // coupling can absorb.
  const Flavour rhos[3] = { rho0, rhop, rhom };
  const Flavour pis [3] = { pi0 , pim , pip  };
  const double  iso [3] = { s["rho0_isospin"].SetDefault(1.).Get<double>(),
			    -1. , 1. };
  for (size_t i(0);i<3 && !contact3pi;++i) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    m_v.back().AddParticle(rhos[i]);
    m_v.back().AddParticle(pis[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("AVVP");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl_grp*iso[i]);
    m_v.back().order[1]=1;
  }
  // rho -> pi pi.  Vector first, as for the photon vertices above.
  const Flavour rho2[3] = { rho0, rhop, rhom };
  const Flavour pa  [3] = { pip , pim , pip  };
  const Flavour pb  [3] = { pim , pi0 , pi0  };
  for (size_t i(0);i<3;++i) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(rho2[i]);
    m_v.back().AddParticle(pa[i]);
    m_v.back().AddParticle(pb[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    if (!rhopipiff.empty()) m_v.back().FormFactor.push_back(rhopipiff);
    m_v.back().cpl.push_back(cpl_rpp);
    m_v.back().order[1]=1;
  }

  ////////////////////////////////////////////////////////////////////
  // hep-ph/0512180 Eq. (10) terms E and F: gamma -> rho(1450) pi and
  // gamma -> rho(1700) pi, each carrying ONLY its partner isoscalar pole
  // (phi(1020) for E, omega(1631) for F) through its own VMD form factor -
  // see FFVMD::ConstructVectorPionFormFactorRhoE/F in Meson_FormFactors.C
  // and the ThreePionFF::Phi1020RhoE()/Omega1631RhoF() note in
  // FormFactorParameters.C for why these cannot be folded into the existing
  // V_pi (rho(770)) production vertex above.
  //
  // Same production coupling cpl_grp as rho(770) - Eq. (10)'s overall scale
  // A is already the anchor for v_omega782's Amplitude = 1.0, so the E/F
  // magnitude (relative to A) lives entirely in Phi1020RhoE()/
  // Omega1631RhoF()'s Amplitude/Phase, exactly as C and D do for the
  // existing omega(1375)/omega(1631) terms in V_pi - up to the chain
  // normalisation below.
  //
  // Eq. (10)'s H_rho sums R = m^2/(s-m^2+...) for each rho, so every chain is
  // -1 at s = 0 whatever its mass, while a Comix chain is
  // g_{rho pi pi}/(s-m^2+...). The rho(1450)/rho(1700) production couplings
  // therefore carry (m_n/m_rho)^2 g_rhopipi/g_{rho_n pi pi}, which puts their
  // chains on eq. (10)'s footing relative to the rho(770) one. Without it the
  // E/F chains are a factor ~3.7 too small.
  //
  // rho(1450)/rho(1700) -> pi pi couplings are read back from their
  // particle-table widths with the SAME formula as g_rhopipi above -
  // Gamma = g^2 p^3/(6 pi M^2) makes no assumption about which rho it is.
  ////////////////////////////////////////////////////////////////////
  const Flavour rho1450_0(kf_rho_1450);
  const Flavour rho1450_p(Flavour(kf_rho_1450_plus));
  const Flavour rho1450_m(Flavour(kf_rho_1450_plus).Bar());
  const Flavour rho1700_0(kf_rho_1700);
  const Flavour rho1700_p(Flavour(kf_rho_1700_plus));
  const Flavour rho1700_m(Flavour(kf_rho_1700_plus).Bar());

  const double mrho1450(rho1450_0.Mass());
  const double grpp1450(RhoPiPiCoupling(rho1450_0,"g_rho1450pipi"));
  const Kabbala cpl_rpp1450("g_{\\rho(1450)\\pi\\pi}",grpp1450*Complex(0.,1.));

  const double mrho1700(rho1700_0.Mass());
  const double grpp1700(RhoPiPiCoupling(rho1700_0,"g_rho1700pipi"));
  const Kabbala cpl_rpp1700("g_{\\rho(1700)\\pi\\pi}",grpp1700*Complex(0.,1.));

  msg_Out()<<METHOD<<": g_rho1450pipi = "<<grpp1450
	   <<", g_rho1700pipi = "<<grpp1700<<"\n";

  const Kabbala norm1450("N_{\\rho(1450)}",sqr(mrho1450/mrho)*grpp/grpp1450);
  const Kabbala norm1700("N_{\\rho(1700)}",sqr(mrho1700/mrho)*grpp/grpp1700);

  // Skipped in contact mode: these are gamma -> rho(1450)/rho(1700) pi
  // PRODUCTION vertices, so leaving them on would open a second route to
  // pi+ pi- pi0 alongside the contact current and double count.  The
  // rho(1450)/rho(1700) -> pi pi decay vertices below stay, like the
  // rho(770) ones, because those resonances are still needed elsewhere.
  // same isospin factors iso[] as the rho(770) production vertices
  const Flavour rhosE[3] = { rho1450_0, rho1450_p, rho1450_m };
  const Flavour rhosF[3] = { rho1700_0, rho1700_p, rho1700_m };
  for (size_t i(0);i<3 && !contact3pi;++i) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    m_v.back().AddParticle(rhosE[i]);
    m_v.back().AddParticle(pis[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("AVVP");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl_grp*norm1450*iso[i]);
    m_v.back().order[1]=1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    m_v.back().AddParticle(rhosF[i]);
    m_v.back().AddParticle(pis[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("AVVP");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl_grp*norm1700*iso[i]);
    m_v.back().order[1]=1;
  }

  const Flavour rho1450_2[3] = { rho1450_0, rho1450_p, rho1450_m };
  const Flavour rho1700_2[3] = { rho1700_0, rho1700_p, rho1700_m };
  for (size_t i(0);i<3;++i) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(rho1450_2[i]);
    m_v.back().AddParticle(pa[i]);
    m_v.back().AddParticle(pb[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    if (!rhopipiff.empty()) m_v.back().FormFactor.push_back(rhopipiff);
    m_v.back().cpl.push_back(cpl_rpp1450);
    m_v.back().order[1]=1;

    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(rho1700_2[i]);
    m_v.back().AddParticle(pa[i]);
    m_v.back().AddParticle(pb[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    if (!rhopipiff.empty()) m_v.back().FormFactor.push_back(rhopipiff);
    m_v.back().cpl.push_back(cpl_rpp1700);
    m_v.back().order[1]=1;
  }

  if (!contact3pi) InitRhoOmegaMixingVertices(cpl_grp);
}

//////////////////////////////////////////////////////////////////////
// g_{rho pi pi} of a rho (770, 1450, 1700) from its width,
//   Gamma = g^2 p^3/(6 pi M^2),  p = sqrt(M^2/4 - m_pi^2),
// THREE_PION_FORM_FACTOR: key overrides it.
//////////////////////////////////////////////////////////////////////
double LowEnergy_Model::RhoPiPiCoupling(const Flavour & rho,const std::string & key) {
  const double M(rho.Mass()), G(rho.Width(true)), mpi(Flavour(kf_pi_plus).Mass());
  const double p(sqrt(Max(0.25*sqr(M)-sqr(mpi),1.e-12)));
  const double def(sqrt(6.*M_PI*sqr(M)*G/pow(p,3)));
  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  return s[key].SetDefault(def).Get<double>();
}

//////////////////////////////////////////////////////////////////////
// g_{gamma rho pi} from the WZW anomaly, 2 g_{gamma rho pi} g_{rho pi pi}/m_rho^2
// = alpha_A, alpha_A = F_3pi/3 x 1.066, F_3pi = 1/(4 pi^2 F_pi^3)
// (arXiv:1907.01556 eq. (2.13)); THREE_PION_FORM_FACTOR: g_gammarhopi overrides it.
//////////////////////////////////////////////////////////////////////
double LowEnergy_Model::GammaRhoPiCoupling() {
  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"] };
  const Flavour rho0(kf_rho_770);
  const double grpp(RhoPiPiCoupling(rho0,"g_rhopipi"));
  // F_pi = 130.7 MeV (PDG 2004 f_pi+), F = F_pi/sqrt(2) in the WZW normalisation
  const double fpi_wzw(s["F_pi"].SetDefault(0.1307).Get<double>()/sqrt(2.));
  const double F3pi(1./(4.*sqr(M_PI)*pow(fpi_wzw,3.)));
  const double alphaA(F3pi/3.*1.066);
  return s["g_gammarhopi"].SetDefault(alphaA*sqr(rho0.Mass())/(2.*grpp)).Get<double>();
}

//////////////////////////////////////////////////////////////////////
// THREE_PION_FORM_FACTOR: RhoOmegaMixing: {Active: false, Amplitude, Phase,
// Br_pipi}. gamma* -> omega pi0 (AVVP, form factor RhoOmegaMix), omega ->
// pi+ pi- (SSV), built like the rho0 pi0 chain, with which it interferes.
// Amplitude and Phase are relative to the rho0 pi0 chain; g_omegapipi follows
// from Gamma(omega -> pi+ pi-) = Br_pipi Gamma_omega = g^2 p^3/(6 pi m^2).
//////////////////////////////////////////////////////////////////////
void LowEnergy_Model::InitRhoOmegaMixingVertices(const Kabbala & cpl_grp) {
  Scoped_Settings s{ Settings::GetMainSettings()["THREE_PION_FORM_FACTOR"]
		     ["RhoOmegaMixing"] };
  if (!s["Active"].SetDefault(false).Get<bool>()) return;
  // Amplitude, Phase: PHOKHARA 10's rho-omega mixing term of its pion=8 current
  // (hep-ph/0512180) in the normalisation of this chain, matched at the phi by
  // yfs-helpers/three-pion-current/rhoomega-mix.py. Br_pipi: PDG Br(omega -> pi+ pi-).
  const double amplitude(s["Amplitude"].SetDefault(3.40726).Get<double>());
  const double phase(s["Phase"].SetDefault(3.02625).Get<double>());
  const double br_pipi(s["Br_pipi"].SetDefault(0.0153).Get<double>());
  const Flavour photon(kf_photon), omega(kf_omega_782), pi0(kf_pi);
  const Flavour pip(kf_pi_plus), pim(Flavour(kf_pi_plus).Bar());
  const double p(sqrt(sqr(omega.Mass())/4.-sqr(pip.Mass())));
  const double g_ompp(sqrt(br_pipi*omega.Width()*6.*M_PI*sqr(omega.Mass())/pow(p,3)));
  const Complex mix(std::polar(amplitude,phase));
  msg_Out()<<METHOD<<": rho-omega mixing on, A e^{i phi} = "<<mix
	   <<", g_omegapipi = "<<g_ompp<<"\n";

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(photon);
  m_v.back().AddParticle(omega);
  m_v.back().AddParticle(pi0);
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("AVVP");
  m_v.back().FormFactor.push_back("RhoOmegaMix");
  m_v.back().cpl.push_back(cpl_grp*Kabbala("A_{\\rho\\omega}",mix));
  m_v.back().order[1]=1;

  m_v.push_back(Single_Vertex());
  m_v.back().AddParticle(omega);
  m_v.back().AddParticle(pip);
  m_v.back().AddParticle(pim);
  m_v.back().Color.push_back(Color_Function(cf::None));
  m_v.back().Lorentz.push_back("SSV");
  m_v.back().cpl.push_back(Kabbala("g_{\\omega\\pi\\pi}",g_ompp*Complex(0.,1.)));
  m_v.back().order[1]=1;
}

//////////////////////////////////////////////////////////////////////
// gamma* -> 4 pi (arXiv:hep-ph/0008262, arXiv:0804.0359), FourPion_Current.H:
// one contact vertex gamma* 4 pi per charge mode with the current linear in
// the pion momenta, J^mu = sum_i G_i p_i^mu (VPPPP, form factor FourPion), the
// omega pi0 term of pi+ pi- 2pi0 included in contact form. Born and initial-
// state radiation only: photon emission off the pions needs the photon
// contacts of the current, which are not built here.
//
// Every 4 pi vertex counts one order more than its photon, so that the Born
// e+e- -> 4 pi is QED order 3 (as with the omega pi0 chain of the full model).
//
// FOUR_PION_FORM_FACTOR: {Active: false, parameters in FourPion_Parameters::Read}.
//////////////////////////////////////////////////////////////////////
void LowEnergy_Model::InitFourPionVertices(const Kabbala & cpl) {
  Scoped_Settings s{ Settings::GetMainSettings()["FOUR_PION_FORM_FACTOR"] };
  if (!s["Active"].SetDefault(false).Get<bool>()) return;
  // contact vertices: Comix's own phase space (binary vertices) cannot map them
  Scoped_Settings integ{ Settings::GetMainSettings()["INTEGRATOR"] };
  bool rambo(false);
  for (const std::string &i : integ.SetDefault("Default").GetVector<std::string>())
    rambo = rambo || i=="Rambo";
  if (!rambo)
    THROW(fatal_error,"FOUR_PION_FORM_FACTOR builds five-leg contact vertices, "
	  "which Comix's phase space cannot map; set INTEGRATOR: Rambo.");
  const Flavour photon(kf_photon), pi0(kf_pi);
  const Flavour pip(kf_pi_plus), pim(Flavour(kf_pi_plus).Bar());
  const std::vector<Flavour> modes[2] = { { pip, pim, pi0, pi0 }, { pip, pip, pim, pim } };
  for (const std::vector<Flavour> &pis : modes) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    for (const Flavour &f : pis) m_v.back().AddParticle(f);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("VPPPP");
    m_v.back().FormFactor.push_back("FourPion");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=2;
  }
}


//////////////////////////////////////////////////////////////////////
//
//
//   L = (e^2/(32 pi^2 f_pi)) eps^{mu,nu,al,be} F_{mu,nu} F_{al,be} pi0,
//
//
// The two photons are generally both off shell, so the form factor here has to
// be F(q1^2, q2^2) - the single-virtuality VMD form factor used everywhere
// else in this model will not do.  Hence the separate "TFF" implementation in
// Meson_FormFactors.C.
//
// DOUBLE COUNTING against the three-pion chain.  Together with the SSV photon
// vertex on the charged pions above, this vertex opens a second route to
// pi+ pi- pi0: gamma* -> gamma* pi0, then gamma* -> pi+ pi- through the pion
// form factor.  With the rho pole in both form factors that is the same
// gamma* -> rho pi0 -> pi+ pi- pi0 transition InitThreePionVertices already
// describes, so in that final state the two add coherently and double count -
// the same relationship THREE_PION_CONTACT toggles.  Comix cannot be told to
// use a vertex in one final state and not another, so the switch is global:
// set TWO_PHOTON_FORM_FACTOR: Active: false to get the tuned three-pion
// channel back unchanged.
//////////////////////////////////////////////////////////////////////
void LowEnergy_Model::InitTwoPhotonVertices() {
  Scoped_Settings s{ Settings::GetMainSettings()["TWO_PHOTON_FORM_FACTOR"] };
  if (!s["Active"].SetDefault(true).Get<bool>()) {
    msg_Out()<<METHOD<<": switched off.\n";
    return;
  }
  ////////////////////////////////////////////////////////////////////
  // Only pi0, eta and eta' get this vertex.  gamma gamma is neutral, C even
  // and carries no strangeness, so of the eight entries in m_PseudoScalars the
  // charged ones (pi+, K+) are forbidden by charge conservation and the
  // neutral kaons (K, K_L, K_S) by strangeness - K_S -> gamma gamma exists
  // only as a weak decay, at a branching ratio of 2.6e-6, which is not what
  // this vertex models.  So this must NOT be written as a loop over
  // m_PseudoScalars filtered on IntCharge()==0, the way InitQEDVertices does
  // it above: that would pick up all three kaons.
  //
  // Each coupling is fixed by the measured two-photon width, via
  //   Gamma(P -> gamma gamma) = g_P^2 m_P^3 / (64 pi),   g_P = alpha/(pi f_P),
  // so the f_P below are back-solved from the PDG widths (7.82 eV, 0.516 keV,
  // 4.34 keV).  For eta and eta' this "f" is the effective constant defined by
  // F_P(0) = 1/(4 pi^2 f_P) - CLEO's Eqn. 6 - not the physical decay constant,
  // which involves eta-eta' mixing.  CLEO's own numbers, 97.5 and 74.4 MeV,
  // came from the 1996 widths; anchoring to the current widths instead moves
  // f_eta by 5%, so do not mix the two sources.
  //
  // NOTE the pi0 default below is 82.4 MeV, not the 91.9 MeV the PDG width
  // implies - it is a tune, and it makes pi0 inconsistent with the eta and
  // eta' entries, which are pure PDG.  Left as set; see NOTES-pi0gg.md.
  ////////////////////////////////////////////////////////////////////
  struct Anomalous_PS {
    kf_code       kf;
    const char *  tag;    // run-card sub-block, and the form factor's key
    double        f;      // GeV
    const char *  tex;
  };
  const Anomalous_PS pseudoscalars[3] = {
    { kf_pi,            "pi0" , 0.09238, "\\pi^0"   },
    { kf_eta,           "eta" , 0.0925, "\\eta"    },
    { kf_eta_prime_958, "eta'", 0.0737, "\\eta'"   }
  };
  for (const Anomalous_PS & ps : pseudoscalars) {
    Scoped_Settings ms{ s[ps.tag] };
    const double f(ms["f"].SetDefault(ps.f).Get<double>());
    const double g_def(m_alpha/(M_PI*f));
    const double g(ms["g_gammagammaP"].SetDefault(g_def).Get<double>());
    msg_Out()<<METHOD<<": "<<ps.tag<<": f = "<<f
	     <<", g_gammagammaP = "<<g<<" /GeV\n";
    const Kabbala cpl_ggp(std::string("g_{\\gamma\\gamma")+ps.tex+"}",
			  g*Complex(0.,1.));
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().AddParticle(Flavour(ps.kf));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("AVVP");
    m_v.back().FormFactor.push_back("TFF");
    m_v.back().cpl.push_back(cpl_ggp);
    m_v.back().order[1]=2;
  }
}


void LowEnergy_Model::InitEWVertices() {}
  
