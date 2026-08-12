#include "MODEL/SM/LowEnergy_Model.H"
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
  (*p_constants)["Lambda2_3122_3122_22"]  =  1.04;
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
  // the photon vertices reproduces scalar QED.
  ////////////////////////////////////////////////////////////////////
  const double mrho(rho0.Mass()), wrho(rho0.Width(true)), mpi(pip.Mass());
  const double prho(sqrt(Max(0.25*sqr(mrho)-sqr(mpi),1.e-12)));
  const double grpp_def(sqrt(6.*M_PI*sqr(mrho)*wrho/pow(prho,3)));
  const double grpp(s["g_rhopipi"].SetDefault(grpp_def).Get<double>());

 
  ////////////////////////////////////////////////////////////////////
  // const double ee(sqrt(4.*M_PI*m_alpha));
  // const double ggrp_vmd(ee*14.0/17.0);   // 0.25 /GeV, for reference
  const double ggrp_def(0.773);
  const double ggrp(s["g_gammarhopi"].SetDefault(ggrp_def).Get<double>());

  msg_Out()<<METHOD<<": g_rhopipi = "<<grpp
	   <<", g_gammarhopi = "<<ggrp<<" /GeV\n";
  const Kabbala cpl_grp("g_{\\gamma\\rho\\pi}",
			ggrp*Complex(0.,1.));
  const Kabbala cpl_rpp("g_{\\rho\\pi\\pi}",
			grpp*Complex(0.,1.));

  // gamma rho pi, anomalous.  Charges have to add up to zero: all legs are
  // read as incoming.
  const Flavour rhos[3] = { rho0, rhop, rhom };
  const Flavour pis [3] = { pi0 , pim , pip  };
  for (size_t i(0);i<3;++i) {
    m_v.push_back(Single_Vertex());
    m_v.back().AddParticle(photon);
    m_v.back().AddParticle(rhos[i]);
    m_v.back().AddParticle(pis[i]);
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("AVVP");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl_grp);
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
    m_v.back().cpl.push_back(cpl_rpp);
    m_v.back().order[1]=1;
  }
}


void LowEnergy_Model::InitEWVertices() {}
  
