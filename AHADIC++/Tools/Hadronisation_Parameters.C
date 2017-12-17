#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Multiplet_Constructor.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;


Hadronisation_Parameters * AHADIC::hadpars = NULL;

Hadronisation_Parameters::Hadronisation_Parameters() {}

Hadronisation_Parameters::~Hadronisation_Parameters() {
  if (p_constituents!=NULL) {
    delete p_constituents;
    p_constituents=NULL;
  }
  if (p_stransitions!=NULL) {
    delete p_stransitions;
    p_stransitions=NULL;
  }
  if (p_dtransitions!=NULL) {
    delete p_dtransitions;
    p_dtransitions=NULL;
  }
}

void Hadronisation_Parameters::Init(string dir,string file)
{
  msg_Out()<<"In Hadronisation_Parameters::Init("<<dir<<file<<")"<<endl;
  ReadParameters(dir,file);

  bool test      = false;
  bool diquarks  = true;
  p_constituents = new Constituents(diquarks);
  Multiplet_Constructor multipletconstructor(false);
  Wave_Functions * wavefunctions = multipletconstructor.GetWaveFunctions();
  p_stransitions = new Single_Transitions(wavefunctions);
  p_dtransitions = new Double_Transitions(p_stransitions);
  
  if (test) {
    p_constituents->PrintConstituents();
    //multipletconstructor.PrintWaveFunctions(true); 
    multipletconstructor.PrintMultiplets();
    p_stransitions->Print();
    p_dtransitions->Print(true);
    //p_dtransitions->Print(false);
    exit(1);
  }
}


void Hadronisation_Parameters::ReadParameters(string dir,string file)
{
  Data_Reader dataread(" ",";","!","=");
  dataread.AddComment("#");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(dir);
  dataread.SetInputFile(file);
  ReadGeneralSwitches(dataread);
  ReadMassParameters(dataread);
  ReadPoppingParameters(dataread);
  ReadMesonWeights(dataread);
  ReadGluonSplittingParameters(dataread);
  ReadClusterDecayParameters(dataread);
  ReadClusterToMesonParameters(dataread);
}

double Hadronisation_Parameters::Get(string keyword) 
{
  m_piter = m_parametermap.find(keyword);
  if (m_piter!=m_parametermap.end()) return m_piter->second;
  msg_Tracking()<<"Error in Hadronisation_Parameters::Get("<<keyword<<") "
		<<"in "<<m_parametermap.size()<<".\n"
		<<"   Keyword not found. Return 0 and hope for the best.\n";
  return 0.;
}

void Hadronisation_Parameters::
ReadGluonSplittingParameters(Data_Reader & dataread) {
  m_parametermap[string("kt_o")]   =
    dataread.GetValue<double>( "KT_O",      1.00 );
  // gluon fragmentation function
  m_parametermap[string("alphaG")] =
    dataread.GetValue<double>( "ALPHA_G",   2.00 );
  m_parametermap[string("betaG")]  =
    dataread.GetValue<double>( "BETA_G",    2.00 );
  // light quark fragmentation function
  m_parametermap[string("alphaL")] =
    dataread.GetValue<double>( "ALPHA_L",   1.50 );
  m_parametermap[string("betaL")]  =
    dataread.GetValue<double>( "BETA_L",    0.25 );
  m_parametermap[string("gammaL")] =
    dataread.GetValue<double>( "GAMMA_L",   0.50 );
  // di-quark fragmentation function
  m_parametermap[string("alphaD")] =
    dataread.GetValue<double>( "ALPHA_D",   m_parametermap[string("alphaL")] );
  m_parametermap[string("betaD")]  =
    dataread.GetValue<double>( "BETA_D",    m_parametermap[string("betaL")] );
  m_parametermap[string("gammaD")] =
    dataread.GetValue<double>( "GAMMA_D",   m_parametermap[string("gammaL")] );
  // heavy quark fragmentation function
  m_parametermap[string("alphaH")] =
    dataread.GetValue<double>( "ALPHA_H",   0.25 );
  m_parametermap[string("betaH")]  =
    dataread.GetValue<double>( "BETA_H",    1.20 );
  m_parametermap[string("gammaH")] =
    dataread.GetValue<double>( "GAMMA_H",   0.25 );
}

void Hadronisation_Parameters::
ReadClusterDecayParameters(Data_Reader & dataread) {
  // Probably irrelevant as long as they are small.
  // We will probably not have to tune them.
  m_parametermap[string("decay_threshold")] =
    dataread.GetValue<double>("DECAY_THRESHOLD",     0.500);
  m_parametermap[string("piphoton_threshold")] =
    dataread.GetValue<double>("PI_PHOTON_THRESHOLD", 0.150);
  m_parametermap[string("dipion_threshold")] =
    dataread.GetValue<double>("DI_PION_THRESHOLD",   0.300);
  m_parametermap[string("open_threshold")] =
    dataread.GetValue<double>("OPEN_THRESHOLD",      0.100);
}

void Hadronisation_Parameters::
ReadDeprecatedParameters(Data_Reader & dataread) {
  m_parametermap[string("colour_reconnection_strength")] = 
    dataread.GetValue<double>("COLOUR_RECONNECTION_STRENGTH",0.23);
}

void Hadronisation_Parameters::
ReadClusterToMesonParameters(Data_Reader & dataread) {
  m_parametermap[string("mass_exponent")]   = 
    dataread.GetValue<double>("MASS_EXPONENT",4.0);
}

void Hadronisation_Parameters::ReadMesonWeights(Data_Reader & dataread) 
{
  m_parametermap[string("Singlet_Suppression")]   = 
    dataread.GetValue<double>("SINGLET_SUPPRESSION",  1.0); // 0.8
  m_parametermap[string("Mixing_Angle_0+")]    = 
    dataread.GetValue<double>("Mixing_0+",            -14.1/180.*M_PI);
  m_parametermap[string("Mixing_Angle_1-")]    = 
    dataread.GetValue<double>("Mixing_1-",             36.4/180.*M_PI);
  m_parametermap[string("Mixing_Angle_2+")]    = 
    dataread.GetValue<double>("Mixing_2+",             27.0/180.*M_PI);
  m_parametermap[string("Mixing_Angle_3-")]    = 
    dataread.GetValue<double>("Mixing_3-",             0.5411);
  m_parametermap[string("Mixing_Angle_4+")]    = 
    dataread.GetValue<double>("Mixing_4+",             0.6283);

  // Mesons currently included
  m_parametermap[string("Multiplet_Meson_R0L0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L0_PSEUDOSCALARS", 1.00);
  m_parametermap[string("Multiplet_Meson_R0L0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L0_VECTORS",       1.00);  
  m_parametermap[string("Multiplet_Meson_R0L0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L0_TENSORS2",      0.75); 
  m_parametermap[string("Multiplet_Meson_R0L1S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L1_SCALARS",       0.10);  
  m_parametermap[string("Multiplet_Meson_R0L1S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L1_AXIALVECTORS",  0.10); 
  m_parametermap[string("Multiplet_Meson_R0L2S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L2_VECTORS",       0.25); 
  // Baryons currently included
  m_parametermap[string("Multiplet_Baryon_R0L0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L0_N_1/2",         1.00); 
  m_parametermap[string("Multiplet_Baryon_R1L0S1/2")]   =  
    dataread.GetValue<double>("MULTI_WEIGHT_R1L0_N_1/2",         1.00); 
  m_parametermap[string("Multiplet_Baryon_R2L0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R2L0_N_1/2",         0.20); 
  m_parametermap[string("Multiplet_Baryon_R1_1L0S1/2")] = 
    dataread.GetValue<double>("MULTI_WEIGHT_R1_1L0_N_1/2",       0.20); 
  m_parametermap[string("Multiplet_Baryon_R0L0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_R0L0_DELTA_3/2",     0.25); 
  // Individual hadrons
  m_parametermap[string("eta_modifier")]   = 
    dataread.GetValue<double>("ETA_MODIFIER",                 0.12);
  m_parametermap[string("eta_prime_modifier")]   = 
    dataread.GetValue<double>("ETA_PRIME_MODIFIER",           1.00);
}

void Hadronisation_Parameters::ReadPoppingParameters(Data_Reader & dataread) 
{
  double sfrac; 
  m_parametermap[string("Strange_fraction")]     = sfrac = 
    dataread.GetValue<double>("STRANGE_FRACTION",                0.42);
  m_parametermap[string("Baryon_fraction")]      = 
    dataread.GetValue<double>("BARYON_FRACTION",                 1.00);
  m_parametermap[string("charm_baryon_modifier")]    = 
    dataread.GetValue<double>("CHARM_BARYON_MODIFIER",           1.00);
  m_parametermap[string("beauty_baryon_modifier")]    = 
    dataread.GetValue<double>("BEAUTY_BARYON_MODIFIER",          1.00);
  m_parametermap[string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",                   0.20);  
  m_parametermap[string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",                   0.04);
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",               0.20); 
}

void Hadronisation_Parameters::ReadMassParameters(Data_Reader & dataread) 
{
  m_parametermap[string("minmass2")] = 
    dataread.GetValue<double>("MIN_MASS2",                      0.10);
  m_parametermap[string("Mass_glue")] = 
    dataread.GetValue<double>("M_GLUE",                         0.00);
  Flavour(kf_gluon).SetHadMass(m_parametermap["Mass_glue"]);
  double mud = m_parametermap[string("Mass_updown")] = 
    dataread.GetValue<double>("M_UP_DOWN",                      0.30);
  double ms = m_parametermap[string("Mass_strange")] = 
    dataread.GetValue<double>("M_STRANGE",                      0.40);
  double mc = m_parametermap[string("Mass_charm")] = 
    dataread.GetValue<double>("M_CHARM",                        1.80);
  double mb = m_parametermap[string("Mass_bottom")] = 
    dataread.GetValue<double>("M_BOTTOM",                       5.10);
  double mdiq = m_parametermap[string("Mass_diquark")] = 
    dataread.GetValue<double>("M_DIQUARK_OFFSET",               0.30);
  double bind0 = m_parametermap[string("Mass_bind0")] = 
    dataread.GetValue<double>("M_BIND_0",                       0.12);
  double bind1 = m_parametermap[string("Mass_bind1")] = 
    dataread.GetValue<double>("M_BIND_1",                       0.50);
  Flavour(kf_d).SetHadMass(mud);
  Flavour(kf_u).SetHadMass(mud);
  Flavour(kf_s).SetHadMass(ms);
  Flavour(kf_c).SetHadMass(mc);
  Flavour(kf_b).SetHadMass(mb);
  Flavour(kf_ud_0).SetHadMass((2.*mud+mdiq)*(1.+bind0));
  Flavour(kf_uu_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_ud_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_dd_1).SetHadMass((2.*mud+mdiq)*(1.+bind1));
  Flavour(kf_su_0).SetHadMass((ms+mud+mdiq)*(1.+bind0));
  Flavour(kf_sd_0).SetHadMass((ms+mud+mdiq)*(1.+bind0));
  Flavour(kf_su_1).SetHadMass((ms+mud+mdiq)*(1.+bind1));
  Flavour(kf_sd_1).SetHadMass((ms+mud+mdiq)*(1.+bind1));
  Flavour(kf_ss_1).SetHadMass((2.*ms+mdiq)*(1.+bind1));
}

void Hadronisation_Parameters::ReadGeneralSwitches(Data_Reader & dataread) 
{
  // General switches for operational modes
  m_ana = dataread.GetValue<int>("FRAGMENTATION_ANALYSIS",0)==1;
  m_parametermap[string("colour_reconnections")] = 
    dataread.GetValue<int>("COLOUR_RECONNECTIONS",0);
}


bool Hadronisation_Parameters::
AdjustMomenta(const int n,ATOOLS::Vec4D * moms,const double * masses)
{
  Momenta_Stretcher stretcher;
  if (n==1) return false;
  bool success(true);
  if (n!=2) {
    bool  prepare=false,boost=false;
    Poincare rest;
    Vec4D cms = Vec4D(0.,0.,0.,0.);
    double mass(0.);
    for (int i=0;i<n;i++) {
      cms  += moms[i];
      mass += masses[i];
      if (dabs(moms[i].Abs2())>1.e-6) prepare = true;
    } 
    if (Vec3D(cms).Abs()>1.e-6) { 
      boost = true;
      rest  = Poincare(cms);
      for (int i=0;i<n;i++) rest.Boost(moms[i]);
    }
    if (mass>sqrt(cms.Abs2())) {
      msg_Error()<<"Error in "<<METHOD<<" : "<<"\n"
		 <<"   Total mass = "<<mass<<", "
		 <<"total E = "<<sqrt(cms.Abs2())<<":\n";
      for (int i=0;i<n;i++) {
	msg_Error()<<"   "<<i<<"th mass = "<<masses[i]<<"\n";
      }
      msg_Error()<<"   Will possibly lead to retrying the event.\n";
      return false;
    }
    if (prepare) success = success && stretcher.ZeroThem(0,n,moms,1.e-10);
    if (!success) std::cout<<METHOD<<" failed for ZeroThem(0,"<<n<<").\n";
    success = success && stretcher.MassThem(0,n,moms,masses);
    if (!success) std::cout<<METHOD<<" failed for MassThem(0,"<<n<<").\n";
    if (boost) {
      for (int i=0;i<n;i++) rest.BoostBack(moms[i]);
    }
  }
  else {
    success = stretcher.MassThem(0,n,moms,masses,1.e-10);
    if (!success) std::cout<<METHOD<<" failed for MassThem(0,"<<n<<"), 2nd.\n";
  }
  if (!success && msg->LevelIsDebugging()) {
    msg_Debugging()<<"Error in "<<METHOD<<" : "<<"\n"
		  <<"   Could not shift particles on new shells.\n";
    for (int i=0;i<n;i++) {
      msg_Debugging()<<"   "<<i<<"th mass = "<<masses[i]<<"\n";
    }
    msg_Debugging()<<"   Will lead to retrying the event.\n";
  }
  return success;
}
