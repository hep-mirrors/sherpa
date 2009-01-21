#include "Hadronisation_Parameters.H"
#include "Soft_Cluster_Handler.H"
#include "Flavour.H"
#include "Momenta_Stretcher.H"
#include "MathTools.H"
#include "Data_Reader.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Hadronisation_Parameters AHADIC::hadpars;


Hadronisation_Parameters::Hadronisation_Parameters() :
  p_constituents(NULL),p_multiplets(NULL),
  p_singletransitions(NULL),p_doubletransitions(NULL),
  p_softclusters(NULL), 
  m_asform(asform::constant), p_coupling(NULL), p_splitter(NULL)
{ }

Hadronisation_Parameters::~Hadronisation_Parameters() {
  if (p_constituents!=NULL)      { delete p_constituents;      p_constituents=NULL;       }
  if (p_multiplets!=NULL)        { delete p_multiplets;        p_multiplets=NULL;         }
  if (p_singletransitions!=NULL) { delete p_singletransitions; p_singletransitions=NULL;  }
  if (p_doubletransitions!=NULL) { delete p_doubletransitions; p_doubletransitions=NULL;  }
  if (p_splitter!=NULL)          { delete p_splitter;          p_splitter=NULL;           }
  if (p_coupling!=NULL)          { delete p_coupling;          p_coupling=NULL;           }
  if (p_softclusters!=NULL)      { delete p_softclusters;      p_softclusters=NULL;       }
}

void Hadronisation_Parameters::Init(string dir,string file)
{
  msg_Tracking()<<"In Hadronisation_Parameters::Init("<<dir<<file<<")"<<endl;
  ReadParameters(dir,file);
  p_constituents      = new Constituents(false);
  if (msg_LevelIsTracking()) p_constituents->PrintConstituents();
  p_multiplets        = new All_Hadron_Multiplets();
  if (msg_LevelIsTracking()) p_multiplets->PrintWaveFunctions(); 

  p_singletransitions = new Single_Transitions();
  if (msg_LevelIsTracking()) p_singletransitions->PrintSingleTransitions(); 

  p_doubletransitions = new Double_Transitions();
  if (msg_LevelIsTracking()) p_doubletransitions->PrintDoubleTransitions(); 

  if (m_parametermap[string("pt02")]<=0.) m_asform = asform::GDH_inspired;
  p_coupling          = new Strong_Coupling(m_asform,dabs(m_parametermap[string("pt02")]));
  p_splitter          = new Dipole_Splitter(p_coupling,m_parametermap[string("ptmax")]);

  p_softclusters      = new Soft_Cluster_Handler(p_coupling,p_singletransitions,p_doubletransitions,
						 m_parametermap[string("Offset_C->H")],
						 m_parametermap[string("Offset_C->HH")],
						 m_parametermap[string("C->H_Transition_Factor")],
						 m_parametermap[string("C->HH_Decay_Exponent")],
						 m_parametermap[string("C->HH_Decay_Angle")],
						 m_parametermap[string("Photon_Energy")],
						 m_parametermap[string("leading_particles")],
						 true);
}
  
void Hadronisation_Parameters::ReadParameters(string dir,string file)
{
  Data_Reader dataread(" ",";","!","=");
  dataread.AddWordSeparator("\t");
  dataread.SetInputPath(dir);
  dataread.SetInputFile(file);
  m_parametermap[string("leading_particles")]  = 
    dataread.GetValue<int>("LEADING",0);
  m_parametermap[string("pt02")]               = 
    dataread.GetValue<double>("PT^2_0",-0.36);
  m_parametermap[string("ptmax")]              = 
    dataread.GetValue<double>("PT_MAX",1.0);
  m_parametermap[string("asfix")]              = 
    dataread.GetValue<double>("AS_FIX",1.0);
  m_parametermap[string("Offset_C->H")] =
    dataread.GetValue<double>("TRANSITION_OFFSET",-1.0);      
  m_parametermap[string("Offset_C->HH")] =
    dataread.GetValue<double>("DECAY_OFFSET",0.2);      
  m_parametermap[string("C->H_Transition_Factor")]       = 
    dataread.GetValue<double>("C->H_TRANSITION_FACTOR",1.);
  m_parametermap[string("C->HH_Decay_Exponent")]       = 
    dataread.GetValue<double>("C->HH_DECAY_EXPONENT",8.);
  m_parametermap[string("C->HH_Decay_Angle")]              = 
    dataread.GetValue<double>("C->HH_DECAY_THETA_EXPONENT",2.);
  m_parametermap[string("Photon_Energy")] =
    dataread.GetValue<double>("PHOTON_ENERGY",0.005);      
  m_parametermap[string("Strange_fraction")] =
    dataread.GetValue<double>("STRANGE_FRACTION",0.37);      
  m_parametermap[string("Baryon_fraction")]  = 
    dataread.GetValue<double>("BARYON_FRACTION",0.2);
  m_parametermap[string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",0.40);
  m_parametermap[string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",0.05);    
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",0.25);
  m_parametermap[string("Singlet_Suppression")]   = 
    dataread.GetValue<double>("SINGLET_SUPPRESSION",0.4);
  m_parametermap[string("Multiplet_Meson_L0R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_PSEUDOSCALARS",1.);
  m_parametermap[string("Multiplet_Meson_L0R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_VECTORS",0.55);
  m_parametermap[string("Multiplet_Meson_L0R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS2",0.12);
  m_parametermap[string("Multiplet_Meson_L0R0S3")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS3",0.05);
  m_parametermap[string("Multiplet_Meson_L0R0S4")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS4",0.);
  m_parametermap[string("Multiplet_Meson_L1R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_SCALARS",0.1);
  m_parametermap[string("Multiplet_Meson_L1R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_AXIALVECTORS",0.05);
  m_parametermap[string("Multiplet_Meson_L1R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_TENSORS2",0.01);
  m_parametermap[string("Multiplet_Meson_L2R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L2R0_VECTORS",0.);
  m_parametermap[string("Multiplet_Meson_L3R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L3R0_VECTORS",0.);
  m_parametermap[string("Multiplet_Meson_L0R1S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_SCALARS",0.);
  m_parametermap[string("Multiplet_Meson_L0R1S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_AXIALVECTORS",0.);
  m_parametermap[string("Multiplet_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N_1/2",1.);
  m_parametermap[string("Multiplet_exc_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N*_1/2",0.05);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_1/2",0.);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_3/2",0.);
  m_parametermap[string("Multiplet_Delta_L0R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_DELTA_3/2",0.15);
  m_parametermap[string("Multiplet_exc_Delta_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_DELTA*_3/2",0.0);
  m_parametermap[string("Mass_glue")]          = 
    dataread.GetValue<double>("M_GLUE",0.);
  m_parametermap[string("Mass_down")]          = 
    dataread.GetValue<double>("M_DOWN",0.3);
  m_parametermap[string("Mass_up")]            = 
    dataread.GetValue<double>("M_UP",0.3);
  m_parametermap[string("Mass_strange")]       = 
    dataread.GetValue<double>("M_STRANGE",0.5);
  m_parametermap[string("Mass_charm")]         = 
    dataread.GetValue<double>("M_CHARM",1.8);
  m_parametermap[string("Mass_bottom")]        = 
    dataread.GetValue<double>("M_BOTTOM",5.2);
  m_parametermap[string("Mass_dd1")]           = 
    dataread.GetValue<double>("M_DD_1",2.*Get("Mass_down"));
  m_parametermap[string("Mass_ud0")]           = 
    dataread.GetValue<double>("M_UD_0",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[string("Mass_ud1")]           = 
    dataread.GetValue<double>("M_UD_1",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[string("Mass_uu1")]           = 
    dataread.GetValue<double>("M_UU_1",2.*Get("Mass_up"));
  m_parametermap[string("Mass_sd0")]           = 
    dataread.GetValue<double>("M_SD_0",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[string("Mass_sd1")]           = 
    dataread.GetValue<double>("M_SD_1",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[string("Mass_su0")]           = 
    dataread.GetValue<double>("M_SU_0",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[string("Mass_su1")]           = 
    dataread.GetValue<double>("M_SU_1",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[string("Mass_ss1")]           = 
    dataread.GetValue<double>("M_SS_1",2.*Get("Mass_strange"));
  m_parametermap[string("Mass_cd0")]           = 
    dataread.GetValue<double>("M_CD_0",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[string("Mass_cd1")]           = 
    dataread.GetValue<double>("M_CD_1",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[string("Mass_cu0")]           = 
    dataread.GetValue<double>("M_CU_0",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[string("Mass_cu1")]           = 
    dataread.GetValue<double>("M_CU_1",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[string("Mass_cs0")]           = 
    dataread.GetValue<double>("M_CS_0",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[string("Mass_cs1")]           = 
    dataread.GetValue<double>("M_CS_1",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[string("Mass_cc1")]           = 
    dataread.GetValue<double>("M_CC_1",2.*Get("Mass_charm"));
  m_parametermap[string("Mass_bd0")]           = 
    dataread.GetValue<double>("M_BD_0",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[string("Mass_bd1")]           = 
    dataread.GetValue<double>("M_BD_1",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[string("Mass_bu0")]           = 
    dataread.GetValue<double>("M_BU_0",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[string("Mass_bu1")]           = 
    dataread.GetValue<double>("M_BU_1",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[string("Mass_bs0")]           = 
    dataread.GetValue<double>("M_BS_0",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[string("Mass_bs1")]           = 
    dataread.GetValue<double>("M_BS_1",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[string("Mass_bc0")]           = 
    dataread.GetValue<double>("M_BC_0",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[string("Mass_bc1")]           = 
    dataread.GetValue<double>("M_BC_1",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[string("Mass_bb1")]           = 
    dataread.GetValue<double>("M_BB_1",2.*Get("Mass_bottom"));
  m_parametermap[string("Mixing_Angle_0+")]    = 
    dataread.GetValue<double>("Mixing_0+",-0.301885);
  m_parametermap[string("Mixing_Angle_1-")]    = 
    dataread.GetValue<double>("Mixing_1-",0.95531);
  m_parametermap[string("Mixing_Angle_2+")]    = 
    dataread.GetValue<double>("Mixing_2+",0.4887);
  m_parametermap[string("Mixing_Angle_3-")]    = 
    dataread.GetValue<double>("Mixing_3-",0.5411);
  m_parametermap[string("Mixing_Angle_4+")]    = 
    dataread.GetValue<double>("Mixing_4+",0.6283);
}

double Hadronisation_Parameters::Get(string keyword) 
{
  m_piter = m_parametermap.find(keyword);
  if (m_piter!=m_parametermap.end()) return m_piter->second;
  msg_Tracking()<<"Error in Hadronisation_Parameters::Get("<<keyword<<") in "<<m_parametermap.size()<<endl
		<<"   Keyword not found. Return 0 and hope for the best."<<endl;
  return 0.;
}

bool Hadronisation_Parameters::AdjustMomenta(const int n,ATOOLS::Vec4D * moms,const double * masses)
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
      msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		    <<"   Total mass = "<<mass<<", total E = "<<sqrt(cms.Abs2())<<std::endl;
      for (int i=0;i<n;i++) {
	msg_Tracking()<<"   "<<i<<"th mass = "<<masses[i]<<std::endl;
      }
      msg_Tracking()<<"   Will possibly lead to retrying the event."<<std::endl;
      return false;
    }
    if (prepare) success = success && stretcher.ZeroThem(0,n,moms,1.e-10);
    success = success && stretcher.MassThem(0,n,moms,masses);
    if (boost) {
      for (int i=0;i<n;i++) rest.BoostBack(moms[i]);
    }
  } 
  else {
    success = stretcher.MassThem(0,n,moms,masses,1.e-10);
  }
  if (!success) {
    msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		  <<"   Could not shift particles on new shells."<<std::endl;
    for (int i=0;i<n;i++) {
      msg_Tracking()<<"   "<<i<<"th mass = "<<masses[i]<<std::endl;
    }
    msg_Tracking()<<"   Will lead to retrying the event."<<std::endl;
  }
  return success;
}
