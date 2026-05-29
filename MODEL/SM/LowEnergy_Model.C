#include "MODEL/SM/LowEnergy_Model.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Single_Vertex.H"
#include "METOOLS/Explicit/Form_Factor.H"
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
  msg_Out()<<METHOD<<": 1/alpha = "<<(1./m_alpha)<<"\n";
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
    m_v.back().AddParticle(flav.Bar());
    m_v.back().AddParticle(flav);
    m_v.back().AddParticle(Flavour(kf_photon));
    m_v.back().Color.push_back(Color_Function(cf::None));
    m_v.back().Lorentz.push_back("SSV");
    m_v.back().FormFactor.push_back("VMD");
    m_v.back().cpl.push_back(cpl);
    m_v.back().order[1]=1;
  }
}


void LowEnergy_Model::InitEWVertices() {}
  
