#include "YFS/Main/YFS_Base.H"
#include "ATOOLS/Math/Random.H" 
// #include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"



using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;

YFS_Base::YFS_Base()
{
  // p_yfsFormFact = new YFS::YFS_Form_Factor();
  RegisterDefaults();
  RegisterSettings();  
}

YFS_Base::~YFS_Base() 
{
}


void YFS_Base::RegisterDefaults(){
  m_s = sqr(rpa->gen.Ecms());
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  s["ISR_MODE"].SetDefault(0);
  s["BETA"].SetDefault(0);
  s["SEMI"].SetDefault(0);
  s["VMIN1"].SetDefault(1e-9);
  s["VMAX"].SetDefault(1);
  s["VMIN"].SetDefault(1e-6);
  s["ISR_CUT"].SetDefault(1e-6);
  s["DELTA"].SetDefault(1e-2);
  s["PHOTON_MAX"].SetDefault(100);
  s["LOOP_TOOL"].SetDefault(false);
  s["RS"].SetDefault(false);
  s["FSR"].SetDefault(0);
  s["SET_MASSES"].SetDefault(false);
  s["FILL_BLOB"].SetDefault(true);
  s["FSR_DEBUG"].SetDefault(false);
  s["ISR_DEBUG"].SetDefault(false);
  s["DEBUG_DIR_ISR"].SetDefault("ISR_DEBUG");
  s["DEBUG_DIR_FSR"].SetDefault("FSR_DEBUG");
  s["DEBUG_DIR_NLO"].SetDefault("YFS_NLO_Hist");
  s["TChannel-Cut"].SetDefault(0);
  s["COULOMB"].SetDefault(false);
  s["FSR_CONST_WEIGHT"].SetDefault(false);
  s["HIDE_PHOTONS"].SetDefault(1);
  s["FULL_FORM"].SetDefault(1);
  s["WW_FORM"].SetDefault(0);
  s["WW_BETAT"].SetDefault(0.382);
  s["INT"].SetDefault(0);
  s["CHECK_MASS_REG"].SetDefault(0);
  s["CHECK_POLES"].SetDefault(0);
  s["CHECK_REAL"].SetDefault(0);
  s["CHECK_REAL_REAL"].SetDefault(0);
  s["CHECK_VIRT_BORN"].SetDefault(1);
  s["VIRTUAL_ONLY"].SetDefault(0);
  s["REAL_ONLY"].SetDefault(0);
  s["USE_MODEL_ALPHA"].SetDefault(0);
  s["FSR_BETA"].SetDefault(1);
  s["KKMC_ANG"].SetDefault(1);
  s["FIXED_WEIGHT"].SetDefault(0);
  s["GRIFFIN_MODE"].SetDefault(0);
  s["GRIFFIN_ORDER"].SetDefault(2);
  s["HARD_MIN"].SetDefault(0.);
  s["PHOTON_MASS"].SetDefault(0.1);
  s["CEEX"].SetDefault(0);
  s["Resonance_Max"].SetDefault(10);
  s["N_Photons"].SetDefault(-1);
  s["TChannel"].SetDefault(0);
  //fcc defaults
  s["BES"].SetDefault(0);
  s["BES_SIG_X"].SetDefault(0.132e-2);
  s["BES_SIG_Y"].SetDefault(0.132e-2);
  s["BES_rho"].SetDefault(-0.745e0);
  s["No_Born"].SetDefault(0);
  s["No_Sub"].SetDefault(0);
  s["Sub_Mode"].SetDefault(submode::local);
  s["Resonance_Mode"].SetDefault(0);
  s["No_Flux"].SetDefault(0);
  s["Flux_Mode"].SetDefault(1);
  s["Pole_Flux"].SetDefault(10);
  s["IFI_Sub"].SetDefault(0);
  s["Massless_Sub"].SetDefault(0);
  s["Check_Real_Sub"].SetDefault(0);
  s["Check_RR_Sub"].SetDefault(0);
  s["Min_Photon_Energy"].SetDefault(0);
}

void YFS_Base::RegisterSettings(){
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  m_betaorder = s["BETA"].Get<int>();
  m_mode = s["ISR_MODE"].Get<int>();
  m_isrcut   = s["ISR_CUT"].Get<double>();
  m_vmin = s["VMIN"].Get<double>();
  m_vmin1 = s["VMIN1"].Get<double>();
  m_vmax = s["VMAX"].Get<double>();
  m_nmax  = s["PHOTON_MAX"].Get<int>();
  m_fillblob  = s["FILL_BLOB"].Get<bool>();
  m_looptool  = s["LOOP_TOOL"].Get<bool>();
  m_YFS_RS  = s["RS"].Get<bool>();
  m_setmass = s["SET_MASSES"].Get<bool>();
  m_RealPhotons = s["NO_PHOTONS"].SetDefault(false).Get<bool>();
  m_NReal       = s["REAL_PHOTONS"].SetDefault(-1).Get<int>();
  m_fsrmode  = s["FSR"].Get<int>();
  m_debugDIR_ISR = s["DEBUG_DIR_ISR"].Get<std::string>();
  m_debugDIR_FSR = s["DEBUG_DIR_FSR"].Get<std::string>();
  m_debugDIR_NLO = s["DEBUG_DIR_NLO"].Get<std::string>();
  m_fsr_debug = s["FSR_DEBUG"].Get<bool>();
  m_isr_debug = s["ISR_DEBUG"].Get<bool>();
  m_deltacut = s["DELTA"].Get<double>()*m_isrcut;
  m_coulomb = s["COULOMB"].Get<bool>();
  m_constfsrW = s["FSR_CONST_WEIGHT"].Get<bool>();
  m_hidephotons=s["HIDE_PHOTONS"].Get<int>();
  m_fullform = s["FULL_FORM"].Get<int>();
  m_formWW = s["WW_FORM"].Get<int>();
  m_betatWW = s["WW_BETAT"].Get<double>();
  m_useint = s["INT"].Get<int>();
  m_check_mass_reg = s["CHECK_MASS_REG"].Get<int>();
  m_check_poles = s["CHECK_POLES"].Get<int>();
  m_check_real = s["CHECK_REAL"].Get<int>();
  m_check_rr = s["CHECK_REAL_REAL"].Get<int>();
  m_check_virt_born = s["CHECK_VIRT_BORN"].Get<int>();
  m_virtual_only = s["VIRTUAL_ONLY"].Get<bool>();
  m_real_only = s["REAL_ONLY"].Get<bool>();
  m_use_model_alpha = s["USE_MODEL_ALPHA"].Get<bool>();
  m_use_fsr_beta = s["FSR_BETA"].Get<int>();
  m_semiyfs = s["SEMI"].Get<int>();
  m_kkmcAngles =  s["KKMC_ANG"].Get<int>();
  m_fixed_weight = s["FIXED_WEIGHT"].Get<double>();
  m_griff = s["GRIFFIN_MODE"].Get<int>();
  m_hardmin = s["HARD_MIN"].Get<double>();
  m_photonMass = s["PHOTON_MASS"].Get<double>();
  m_useceex = s["CEEX"].Get<int>();
  m_resonace_max = s["Resonance_Max"].Get<double>();
  m_fixed_photons = s["N_Photons"].Get<int>();
  m_beam_rho = s["BES_rho"].Get<double>();
  m_beam_sigx = s["BES_SIG_X"].Get<double>();
  m_beam_sigy = s["BES_SIG_Y"].Get<double>();
  m_no_born = s["No_Born"].Get<int>();
  m_no_subtraction = s["No_Sub"].Get<int>();
  m_submode = s["Sub_Mode"].Get<submode::code>();
  m_resonance_mode = s["Resonance_Mode"].Get<int>();
  m_tchannel = s["TChannel"].Get<int>();
  m_noflux = s["No_Flux"].Get<int>();
  m_flux_mode=s["Flux_Mode"].Get<int>();
  m_pole_fac = s["Pole_Flux"].Get<double>();
  m_min_photon_E = s["Min_Photon_Energy"].Get<double>();
  m_ifisub = s["IFI_Sub"].Get<int>();
  m_massless_sub = s["Massless_Sub"].Get<int>();
  m_check_real_sub = s["Check_Real_Sub"].Get<bool>();
  m_check_rr_sub = s["Check_RR_Sub"].Get<bool>();
  m_CalForm = false;
  m_realtool = false;
  //update when beamstrahlung is added
  m_beamspread=s["BES"].Get<int>();;
  m_isrinital=true;
  m_g = 0;
  m_gp = 0;
  if(m_use_model_alpha) m_alpha = s_model->ScalarConstant("alpha_QED");
  else m_alpha  = 1./s["1/ALPHAQED(0)"].SetDefault(137.03599976).Get<double>(); 
  if (m_use_model_alpha) m_rescale_alpha = 1;
  else m_rescale_alpha = m_alpha / s_model->ScalarConstant("alpha_QED");
  m_alpi = m_alpha/M_PI;
  m_beam_sigx*=sqrt(m_s)/2.;
  m_beam_sigy*=sqrt(m_s)/2.;
  // p_yfsFormFact = new YFS::YFS_Form_Factor();

}

std::istream &YFS::operator>>(std::istream &str,submode::code &sub)
{
  std::string tag;
  str>>tag;
  sub=submode::local;
  if      (tag.find("Off")!=std::string::npos)  sub=submode::off;
  else if (tag.find("0")!=std::string::npos)    sub=submode::off;
  else if (tag.find("Local")!=std::string::npos)sub=submode::local;
  else if (tag.find("1")!=std::string::npos)   sub=submode::local;
  else if (tag.find("Global")!=std::string::npos) sub=submode::global;
  else if (tag.find("2")!=std::string::npos)    sub=submode::global;
  return str;
}


double YFS_Base::Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double YFS_Base::EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  // return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
  return m_alpha/(4*M_PI*M_PI)*(2*p1*p2/((p1*k)*(p2*k)));
}
