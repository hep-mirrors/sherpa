#include "Hadronisation_Parameters.H"
#include "Flavour.H"
#include "Momenta_Stretcher.H"
#include "MathTools.H"
#include "Data_Read.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Hadronisation_Parameters AHADIC::hadpars;


Hadronisation_Parameters::Hadronisation_Parameters() :
  p_constituents(NULL),p_multiplets(NULL),
  p_singletransitions(NULL),p_doubletransitions(NULL)
{ }

Hadronisation_Parameters::~Hadronisation_Parameters() {
  if (p_constituents!=NULL)      { delete p_constituents;      p_constituents=NULL;       }
  if (p_multiplets!=NULL)        { delete p_multiplets;        p_multiplets=NULL;         }
  if (p_singletransitions!=NULL) { delete p_singletransitions; p_singletransitions=NULL;  }
  if (p_doubletransitions!=NULL) { delete p_doubletransitions; p_doubletransitions=NULL;  }
  if (p_popper!=NULL)            { delete p_popper;            p_popper=NULL;             }
}

void Hadronisation_Parameters::Init(string dir,string file)
{
  msg.Tracking()<<"In Hadronisation_Parameters::Init("<<dir<<file<<")"<<endl;
  ReadParameters(dir,file);
  p_constituents = new Constituents(false);
  // if (msg.LevelIsTracking()) p_constituents->PrintConstituents();

  p_multiplets   = new All_Hadron_Multiplets();
  //if (msg.LevelIsTracking()) p_multiplets->PrintWaveFunctions(); 

  p_singletransitions  = new Single_Transitions();
  //if (msg.LevelIsTracking()) p_singletransitions->PrintSingleTransitions(); 

  p_doubletransitions  = new Double_Transitions();
  //if (msg.LevelIsTracking()) p_doubletransitions->PrintDoubleTransitions(); 

  p_popper       = new Pair_Popper();
}
  
void Hadronisation_Parameters::ReadParameters(string dir,string file)
{
  Data_Read dataread(dir+file);
  m_parametermap[string("Tension")]            = 
    dataread.GetValue<double>("COLOUR_TENSION",0.33);
  m_parametermap[string("<pt shift>")]         = 
    dataread.GetValue<double>("<PT_SHIFT>",0.5);
  m_parametermap[string("<Y>")]                = 
    dataread.GetValue<double>("<Y*_SHIFT>",0.5);
  m_parametermap[string("Y*_WIDTH")]           = 
    dataread.GetValue<double>("Y*_WIDTH",dabs(Get("<Y>"))>1.e-3?Get("<Y>"):0.5);
  m_parametermap[string("C->HH_Decay_Exponent")]       = 
    dataread.GetValue<double>("C->HH_DECAY_EXPONENT",0.5);
  m_parametermap[string("MassFraction")]       = 
    dataread.GetValue<double>("CLUSTER_MASS_FRACTION",0.5);
  m_parametermap[string("Offset_C->H")] =
    dataread.GetValue<double>("TRANSITION_OFFSET",0.75);      
  m_parametermap[string("Offset_C->HH")] =
    dataread.GetValue<double>("DECAY_OFFSET",0.5);      
  m_parametermap[string("Photon_Energy")] =
    dataread.GetValue<double>("PHOTON_ENERGY",hadpars.Get((string("Offset_C->HH"))));      
  m_parametermap[string("Strange_fraction")] =
    dataread.GetValue<double>("STRANGE_FRACTION",0.2);      
  m_parametermap[string("Baryon_fraction")]  = 
    dataread.GetValue<double>("BARYON_FRACTION",0.25);
  m_parametermap[string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",0.5);
  m_parametermap[string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",0.1);    
  m_parametermap[string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",1.);
  m_parametermap[string("Multiplet_Meson_L0R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_PSEUDOSCALARS",1.);
  m_parametermap[string("Multiplet_Meson_L0R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_VECTORS",0.5);
  m_parametermap[string("Multiplet_Meson_L0R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS2",0.1);
  m_parametermap[string("Multiplet_Meson_L0R0S3")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS3",0.);
  m_parametermap[string("Multiplet_Meson_L0R0S4")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_TENSORS4",0.);
  m_parametermap[string("Multiplet_Meson_L1R0S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_SCALARS",0.05);
  m_parametermap[string("Multiplet_Meson_L1R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_AXIALVECTORS",0.025);
  m_parametermap[string("Multiplet_Meson_L1R0S2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_TENSORS2",0.);
  m_parametermap[string("Multiplet_Meson_L2R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L2R0_VECTORS",0.01);
  m_parametermap[string("Multiplet_Meson_L3R0S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L3R0_VECTORS",0.);
  m_parametermap[string("Multiplet_Meson_L0R1S0")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_SCALARS",0.);
  m_parametermap[string("Multiplet_Meson_L0R1S1")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R1_AXIALVECTORS",0.);
  m_parametermap[string("Multiplet_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N_1/2",1.);
  m_parametermap[string("Multiplet_exc_Nucleon_L0R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_N*_1/2",0.);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S1/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_1/2",0.);
  m_parametermap[string("Multiplet_exc_Nucleon_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_N*_3/2",0.);
  m_parametermap[string("Multiplet_Delta_L0R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L0R0_DELTA_3/2",0.6);
  m_parametermap[string("Multiplet_exc_Delta_L1R0S3/2")]   = 
    dataread.GetValue<double>("MULTI_WEIGHT_L1R0_DELTA*_3/2",0.6);
  m_parametermap[string("Mass_glue")]          = 
    dataread.GetValue<double>("M_GLUE",0.75);
  m_parametermap[string("Mass_down")]          = 
    dataread.GetValue<double>("M_DOWN",0.32);
  m_parametermap[string("Mass_up")]            = 
    dataread.GetValue<double>("M_UP",0.32);
  m_parametermap[string("Mass_strange")]       = 
    dataread.GetValue<double>("M_STRANGE",0.45);
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
    dataread.GetValue<double>("Mixing_0+",-0.4294);
  m_parametermap[string("Mixing_Angle_1-")]    = 
    dataread.GetValue<double>("Mixing_1-",0.6283);
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
  msg.Error()<<"Error in Hadronisation_Parameters::Get("<<keyword<<") in "<<m_parametermap.size()<<endl
	     <<"   Keyword not found. Return 0 and hope for the best."<<endl;
  return 0.;
}

bool Hadronisation_Parameters::AdjustMomenta(const int n,ATOOLS::Vec4D * moms,const double * masses)
{
  Momenta_Stretcher stretcher;
  if (n==1) return false;
  if (n!=2) {
    bool  prepare=false,boost=false,success=true;
    Poincare rest;
    Vec4D cms = Vec4D(0.,0.,0.,0.);
    for (int i=0;i<n;i++) {
      cms += moms[i];
      if (dabs(moms[i].Abs2())>1.e-6) prepare = true;
    } 
    if (Vec3D(cms).Abs()>1.e-3) { 
      boost = true;
      rest  = Poincare(cms);
      for (int i=0;i<n;i++) rest.Boost(moms[i]);
    }
    if (prepare) success = success && stretcher.ZeroThem(0,n,moms);
    success = success && stretcher.MassThem(0,n,moms,masses);
    if (boost) {
      for (int i=0;i<n;i++) rest.BoostBack(moms[i]);
    }
    return success;
  } 
  else return stretcher.MassThem(0,n,moms,masses);
}
