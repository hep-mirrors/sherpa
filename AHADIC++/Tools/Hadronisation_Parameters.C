#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
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


Hadronisation_Parameters *AHADIC::hadpars=NULL;


Hadronisation_Parameters::Hadronisation_Parameters() :
  p_constituents(NULL),p_multiplets(NULL),
  p_singletransitions(NULL),p_doubletransitions(NULL),
  p_softclusters(NULL)
{ }

Hadronisation_Parameters::~Hadronisation_Parameters() {
  if (p_constituents!=NULL)      { delete p_constituents;      p_constituents=NULL;       }
  if (p_multiplets!=NULL)        { delete p_multiplets;        p_multiplets=NULL;         }
  if (p_singletransitions!=NULL) { delete p_singletransitions; p_singletransitions=NULL;  }
  if (p_doubletransitions!=NULL) { delete p_doubletransitions; p_doubletransitions=NULL;  }
  if (p_softclusters!=NULL)      { delete p_softclusters;      p_softclusters=NULL;       }
}

void Hadronisation_Parameters::Init(string dir,string file)
{
  msg_Tracking()<<"In Hadronisation_Parameters::Init("<<dir<<file<<")"<<endl;
  ReadParameters(dir,file);

  p_constituents      = new Constituents(false);

  if (msg_LevelIsTracking()) 
    p_constituents->PrintConstituents();

  p_multiplets        = new All_Hadron_Multiplets();
  if (msg_LevelIsTracking()) {
    p_multiplets->PrintWaveFunctions(); 
    p_multiplets->PrintMultiplets();
  }

  p_singletransitions = new Single_Transitions();
  if (msg_LevelIsTracking()) 
    p_singletransitions->PrintSingleTransitions(); 

  p_doubletransitions = new Double_Transitions();
  if (msg_LevelIsTracking()) 
    p_doubletransitions->PrintDoubleTransitions(); 

  p_softclusters      = new Soft_Cluster_Handler(m_ana);
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
  m_parametermap[string("pt02")]                 = 
    dataread.GetValue<double>("PT^2_0",1.00);
  m_parametermap[string("G2QQ_Exponent")]        = 
    dataread.GetValue<double>("G2QQ_EXPONENT",0.50);
  m_parametermap[string("G2QQ_LeadExponent")]    = 
    dataread.GetValue<double>("G2QQ_LEADEXPONENT",1.00);
  m_parametermap[string("ptmax")]                = 
    dataread.GetValue<double>("PT_MAX",1.00);
  m_parametermap[string("ptmax_factor")]         = 
    dataread.GetValue<double>("PT_MAX_FACTOR",1.0);
}

void Hadronisation_Parameters::
ReadClusterDecayParameters(Data_Reader & dataread) {
  m_parametermap[string("MaxNumberOfPairs")]  =
    dataread.GetValue<int>("MAX_PAIRS",1);
  m_parametermap[string("SplitExponent")]     = 
    dataread.GetValue<double>("SPLIT_EXPONENT",0.90);
  m_parametermap[string("SplitLeadExponent")] = 
    dataread.GetValue<double>("SPLIT_LEADEXPONENT",2.00);
  m_parametermap[string("SpectExponent")]     = 
    dataread.GetValue<double>("SPECT_EXPONENT",0.0);
  m_parametermap[string("SpectLeadExponent")] = 
    dataread.GetValue<double>("SPECT_LEADEXPONENT",0.5);
}

void Hadronisation_Parameters::
ReadDeprecatedParameters(Data_Reader & dataread) {
  m_parametermap[string("colour_reconnection_strength")] = 
    dataread.GetValue<double>("COLOUR_RECONNECTION_STRENGTH",0.23);
}

void Hadronisation_Parameters::
ReadClusterToMesonParameters(Data_Reader & dataread) {
  m_parametermap[string("Offset_C->H")]          =
    dataread.GetValue<double>("TRANSITION_OFFSET",0.8);
  m_parametermap[string("MassExponent_C->H")]    =
    dataread.GetValue<double>("TRANSITION_EXPONENT",0.15);
  m_parametermap[string("WidthExponent_C->H")]   =
    dataread.GetValue<double>("TRANSITION_EXPONENT2",-0.32);
  m_parametermap[string("Offset_C->HH")]         =
    dataread.GetValue<double>("DECAY_OFFSET",1.00);
  m_parametermap[string("MassExponent_C->HH")]   =
    dataread.GetValue<double>("DECAY_EXPONENT",6.00);
}

void Hadronisation_Parameters::ReadMesonWeights(Data_Reader & dataread) 
{
  m_parametermap[string("Singlet_Suppression")]   = 
    dataread.GetValue<double>("SINGLET_SUPPRESSION",0.57);
  m_parametermap[string("Multiplet_Meson_L0R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_PSEUDOSCALARS",1.00);
  m_parametermap[string("Multiplet_Meson_L0R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_VECTORS",0.75);
  m_parametermap[string("Multiplet_Meson_L0R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS2",0.30);
  m_parametermap[string("Multiplet_Meson_L0R0S3")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS3",0.00);
  m_parametermap[string("Multiplet_Meson_L0R0S4")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS4",0.00);
  m_parametermap[string("Multiplet_Meson_L1R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_SCALARS",1.00);
  m_parametermap[string("Multiplet_Meson_L1R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_AXIALVECTORS",0.75);
  m_parametermap[string("Multiplet_Meson_L1R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_TENSORS2",0.00);
  m_parametermap[string("Multiplet_Meson_L2R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L2R0_VECTORS",0.80);
  m_parametermap[string("Multiplet_Meson_L3R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L3R0_VECTORS",0.00);
  m_parametermap[string("Multiplet_Meson_L0R1S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_SCALARS",0.50);
  m_parametermap[string("Multiplet_Meson_L0R1S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_AXIALVECTORS",0.00);
  m_parametermap[string("Multiplet_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N_1/2",1.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N*_1/2",1.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_1/2",0.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_3/2",0.00);
  m_parametermap[string("Multiplet_Delta_L0R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_DELTA_3/2",0.45);
  m_parametermap[string("Multiplet_exc_Delta_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_DELTA*_3/2",0.00);
  m_parametermap[string("Mixing_Angle_0+")]    = 
    dataread.GetValue<double>("Mixing_0+",-0.31416)-M_PI/2.;
  m_parametermap[string("Mixing_Angle_1-")]    = 
    dataread.GetValue<double>("Mixing_1-",0.61075);
  m_parametermap[string("Mixing_Angle_2+")]    = 
    dataread.GetValue<double>("Mixing_2+",0.4887);
  m_parametermap[string("Mixing_Angle_3-")]    = 
    dataread.GetValue<double>("Mixing_3-",0.5411);
  m_parametermap[string("Mixing_Angle_4+")]    = 
    dataread.GetValue<double>("Mixing_4+",0.6283);
}

void Hadronisation_Parameters::ReadPoppingParameters(Data_Reader & dataread) 
{
  m_parametermap[string("Strange_fraction")]     =
    dataread.GetValue<double>("STRANGE_FRACTION",0.63);
  m_parametermap[string("Baryon_fraction")]      = 
    dataread.GetValue<double>("BARYON_FRACTION",1.00);
  m_parametermap[string("Heavy_Baryon_Enhancement")]    = 
    dataread.GetValue<double>("HEAVY_BARYON_ENHANCEMENT",4.);
  m_parametermap[string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",0.3);
  m_parametermap[string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",0.01);
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",1.0);
}

void Hadronisation_Parameters::ReadMassParameters(Data_Reader & dataread) 
{
  m_parametermap[string("minmass2")]             = 
    dataread.GetValue<double>("MIN_MASS2",0.0);
  m_parametermap[string("Mass_glue")]          = 
    dataread.GetValue<double>("M_GLUE",0.);
  Flavour(kf_gluon).SetHadMass(m_parametermap["Mass_glue"]);
  double mud = m_parametermap[string("Mass_updown")] = 
    dataread.GetValue<double>("M_UP_DOWN",0.3);
  double ms = m_parametermap[string("Mass_strange")] = 
    dataread.GetValue<double>("M_STRANGE",0.45);
  double mc = m_parametermap[string("Mass_charm")]         = 
    dataread.GetValue<double>("M_CHARM",1.8);
  double mb = m_parametermap[string("Mass_bottom")]        = 
    dataread.GetValue<double>("M_BOTTOM",5.1);
  double mdiq = m_parametermap[string("Mass_diquark")]           = 
    dataread.GetValue<double>("M_DIQUARK_OFFSET",0.12);
  double bind0 = m_parametermap[string("Mass_bind0")]           = 
    dataread.GetValue<double>("M_BIND_0",0.12);
  double bind1 = m_parametermap[string("Mass_bind1")]           = 
    dataread.GetValue<double>("M_BIND_1",0.25);
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
