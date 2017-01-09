#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Default_Reader.H"
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
  Default_Reader reader;
  reader.SetInputPath(dir);
  reader.SetInputFile(file);
  ReadGeneralSwitches(reader);
  ReadMassParameters(reader);
  ReadPoppingParameters(reader);
  ReadMesonWeights(reader);
  ReadGluonSplittingParameters(reader);
  ReadClusterDecayParameters(reader);
  ReadClusterToMesonParameters(reader);
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
ReadGluonSplittingParameters(Default_Reader & reader) {
  m_parametermap[string("pt02")]                 = 
    reader.Get<double>("PT^2_0",1.562);
  m_parametermap[string("G2QQ_Exponent")]        = 
    reader.Get<double>("G2QQ_EXPONENT",1.08);
  m_parametermap[string("G2QQ_LeadExponent")]    = 
    reader.Get<double>("G2QQ_LEADEXPONENT",0.00);
  m_parametermap[string("ptmax")]                = 
    reader.Get<double>("PT_MAX",1.00);
  m_parametermap[string("ptmax_factor")]         = 
    reader.Get<double>("PT_MAX_FACTOR",1.0);
}

void Hadronisation_Parameters::
ReadClusterDecayParameters(Default_Reader & reader) {
  m_parametermap[string("MaxNumberOfPairs")]  =
    reader.Get<int>("MAX_PAIRS",1);
  m_parametermap[string("SplitExponent")]     = 
    reader.Get<double>("SPLIT_EXPONENT",0.1608);
  m_parametermap[string("SplitLeadExponent")] = 
    reader.Get<double>("SPLIT_LEADEXPONENT",8);
  m_parametermap[string("SpectExponent")]     = 
    reader.Get<double>("SPECT_EXPONENT",1.739);
  m_parametermap[string("SpectLeadExponent")] = 
    reader.Get<double>("SPECT_LEADEXPONENT",8);
}

void Hadronisation_Parameters::
ReadDeprecatedParameters(Default_Reader & reader) {
  m_parametermap[string("colour_reconnection_strength")] = 
    reader.Get<double>("COLOUR_RECONNECTION_STRENGTH",0.23);
}

void Hadronisation_Parameters::
ReadClusterToMesonParameters(Default_Reader & reader) {
  m_parametermap[string("Offset_C->H")]          =
    reader.Get<double>("TRANSITION_OFFSET",0.8);
  m_parametermap[string("MassExponent_C->H")]    =
    reader.Get<double>("TRANSITION_EXPONENT",0.15);
  m_parametermap[string("WidthExponent_C->H")]   =
    reader.Get<double>("TRANSITION_EXPONENT2",-0.32);
  m_parametermap[string("Offset_C->HH")]         =
    reader.Get<double>("DECAY_OFFSET",1.202);
  m_parametermap[string("MassExponent_C->HH")]   =
    reader.Get<double>("DECAY_EXPONENT",2.132);
}

void Hadronisation_Parameters::ReadMesonWeights(Default_Reader & reader) 
{
  m_parametermap[string("Singlet_Suppression")]   = 
    reader.Get<double>("SINGLET_SUPPRESSION",0.57);
  m_parametermap[string("Multiplet_Meson_L0R0S0")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_PSEUDOSCALARS",1.00);
  m_parametermap[string("Multiplet_Meson_L0R0S1")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_VECTORS",0.75);
  m_parametermap[string("Multiplet_Meson_L0R0S2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_TENSORS2",0.30);
  m_parametermap[string("Multiplet_Meson_L0R0S3")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_TENSORS3",0.00);
  m_parametermap[string("Multiplet_Meson_L0R0S4")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_TENSORS4",0.00);
  m_parametermap[string("Multiplet_Meson_L1R0S0")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_SCALARS",1.00);
  m_parametermap[string("Multiplet_Meson_L1R0S1")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_AXIALVECTORS",0.75);
  m_parametermap[string("Multiplet_Meson_L1R0S2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_TENSORS2",0.00);
  m_parametermap[string("Multiplet_Meson_L2R0S1")]   = 
    reader.Get<double>("MULTI_WEIGHT_L2R0_VECTORS",0.80);
  m_parametermap[string("Multiplet_Meson_L3R0S1")]   = 
    reader.Get<double>("MULTI_WEIGHT_L3R0_VECTORS",0.00);
  m_parametermap[string("Multiplet_Meson_L0R1S0")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R1_SCALARS",0.50);
  m_parametermap[string("Multiplet_Meson_L0R1S1")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R1_AXIALVECTORS",0.00);
  m_parametermap[string("Multiplet_Nucleon_L0R0S1/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_N_1/2",1.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L0R0S1/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_N*_1/2",1.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S1/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_N*_1/2",0.00);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S3/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_N*_3/2",0.00);
  m_parametermap[string("Multiplet_Delta_L0R0S3/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L0R0_DELTA_3/2",0.45);
  m_parametermap[string("Multiplet_exc_Delta_L1R0S3/2")]   = 
    reader.Get<double>("MULTI_WEIGHT_L1R0_DELTA*_3/2",0.00);
  m_parametermap[string("Mixing_Angle_0+")]    = 
    reader.Get<double>("Mixing_0+",-0.31416)-M_PI/2.;
  m_parametermap[string("Mixing_Angle_1-")]    = 
    reader.Get<double>("Mixing_1-",0.61075);
  m_parametermap[string("Mixing_Angle_2+")]    = 
    reader.Get<double>("Mixing_2+",0.4887);
  m_parametermap[string("Mixing_Angle_3-")]    = 
    reader.Get<double>("Mixing_3-",0.5411);
  m_parametermap[string("Mixing_Angle_4+")]    = 
    reader.Get<double>("Mixing_4+",0.6283);
}

void Hadronisation_Parameters::ReadPoppingParameters(Default_Reader & reader) 
{
  m_parametermap[string("Strange_fraction")]     =
    reader.Get<double>("STRANGE_FRACTION",0.6049);
  m_parametermap[string("Baryon_fraction")]      = 
    reader.Get<double>("BARYON_FRACTION",1.00);
  m_parametermap[string("Heavy_Baryon_Enhancement")]    = 
    reader.Get<double>("HEAVY_BARYON_ENHANCEMENT",4.);
  m_parametermap[string("P_qs_by_P_qq")]       = 
    reader.Get<double>("P_{QS}/P_{QQ}",0.3);
  m_parametermap[string("P_ss_by_P_qq")]       = 
    reader.Get<double>("P_{SS}/P_{QQ}",0.01);
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    reader.Get<double>("P_{QQ_1}/P_{QQ_0}",1.0);
}

void Hadronisation_Parameters::ReadMassParameters(Default_Reader & reader) 
{
  m_parametermap[string("minmass2")]             = 
    reader.Get<double>("MIN_MASS2",0.0);
  m_parametermap[string("Mass_glue")]          = 
    reader.Get<double>("M_GLUE",0.);
  Flavour(kf_gluon).SetHadMass(m_parametermap["Mass_glue"]);
  double mud = m_parametermap[string("Mass_updown")] = 
    reader.Get<double>("M_UP_DOWN",0.3);
  double ms = m_parametermap[string("Mass_strange")] = 
    reader.Get<double>("M_STRANGE",0.45);
  double mc = m_parametermap[string("Mass_charm")]         = 
    reader.Get<double>("M_CHARM",1.8);
  double mb = m_parametermap[string("Mass_bottom")]        = 
    reader.Get<double>("M_BOTTOM",5.1);
  double mdiq = m_parametermap[string("Mass_diquark")]           = 
    reader.Get<double>("M_DIQUARK_OFFSET",0.12);
  double bind0 = m_parametermap[string("Mass_bind0")]           = 
    reader.Get<double>("M_BIND_0",0.12);
  double bind1 = m_parametermap[string("Mass_bind1")]           = 
    reader.Get<double>("M_BIND_1",0.25);
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

void Hadronisation_Parameters::ReadGeneralSwitches(Default_Reader & reader) 
{
  // General switches for operational modes
  m_ana = reader.Get<int>("FRAGMENTATION_ANALYSIS",0)==1;
  m_parametermap[string("colour_reconnections")] = 
    reader.Get<int>("COLOUR_RECONNECTIONS",0);
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
