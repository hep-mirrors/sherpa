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
  s["MODE"].ResetDefault().SetDefault(yfsmode::off);
  s["BETA"].SetDefault(3);
  s["VMAX"].SetDefault(0);
  s["IR_CUTOFF"].ResetDefault().SetDefault(1e-4);
  s["DELTA"].SetDefault(1e-2);
  s["PHOTON_MAX"].SetDefault(100);
  s["LOOP_TOOL"].SetDefault(false);
  s["FILL_BLOB"].SetDefault(true);
  s["FSR_DEBUG"].SetDefault(false);
  s["ISR_DEBUG"].SetDefault(false);
  s["DEBUG_DIR_ISR"].SetDefault("ISR_DEBUG");
  s["DEBUG_DIR_FSR"].SetDefault("FSR_DEBUG");
  s["DEBUG_DIR_NLO"].SetDefault("YFS_NLO_Hist");
  s["TChannel-Cut"].SetDefault(0);
  s["COULOMB"].SetDefault(false);
  s["HIDE_PHOTONS"].SetDefault(1);
  s["FULL_FORM"].SetDefault(0);
  s["WW_FORM"].SetDefault(0);
  s["WW_BETAT"].SetDefault(0.382);
  s["CHECK_MASS_REG"].SetDefault(0);
  s["CHECK_POLES"].SetDefault(0);
  s["CHECK_REAL"].SetDefault(0);
  s["CHECK_REAL_REAL"].SetDefault(0);
  s["CHECK_VIRT_BORN"].SetDefault(0);
  s["VIRTUAL_ONLY"].SetDefault(0);
  s["REAL_ONLY"].SetDefault(0);
  s["USE_MODEL_ALPHA"].SetDefault(0);
  s["KKMC_ANG"].SetDefault(0);
  s["FIXED_WEIGHT"].SetDefault(0);
  s["HARD_MIN"].SetDefault(0.);
  s["PHOTON_MASS"].SetDefault(0.1);
  s["CEEX"].SetDefault(0);
  s["Collinear_Real"].SetDefault(0);
  s["CLUSTERING_THRESHOLD"].SetDefault(10);
  s["TChannel"].SetDefault(0);
  s["No_Born"].SetDefault(0);
  s["No_Sub"].SetDefault(0);
  s["Sub_Mode"].SetDefault(submode::global);
  s["No_Flux"].SetDefault(0);
  s["Flux_Mode"].SetDefault(1);
  s["IFI_Sub"].SetDefault(0);
  s["Massless_Sub"].SetDefault(0);
  s["Check_Real_Sub"].SetDefault(0);
  s["Integrate_NLO"].SetDefault(0);
}

void YFS_Base::RegisterSettings(){
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  m_betaorder = s["BETA"].Get<int>();
  m_mode = s["MODE"].Get<yfsmode::code>();
  m_isrcut   = s["IR_CUTOFF"].Get<double>();
  m_isrcut = m_isrcut/sqrt(m_s); // dimensionless units
  m_vmax = s["VMAX"].Get<double>();
  m_vmax = 1.-sqr(m_vmax)/m_s;
  m_fillblob  = s["FILL_BLOB"].Get<bool>();
  m_looptool  = s["LOOP_TOOL"].Get<bool>();
  m_debugDIR_ISR = s["DEBUG_DIR_ISR"].Get<std::string>();
  m_debugDIR_FSR = s["DEBUG_DIR_FSR"].Get<std::string>();
  m_debugDIR_NLO = s["DEBUG_DIR_NLO"].Get<std::string>();
  m_fsr_debug = s["FSR_DEBUG"].Get<bool>();
  m_isr_debug = s["ISR_DEBUG"].Get<bool>();
  m_deltacut = s["DELTA"].Get<double>()*m_isrcut;
  m_coulomb = s["COULOMB"].Get<bool>();
  m_hidephotons=s["HIDE_PHOTONS"].Get<int>();
  m_fullform = s["FULL_FORM"].Get<int>();
  m_formWW = s["WW_FORM"].Get<int>();
  m_betatWW = s["WW_BETAT"].Get<double>();
  m_check_mass_reg = s["CHECK_MASS_REG"].Get<int>();
  m_check_poles = s["CHECK_POLES"].Get<int>();
  m_check_real = s["CHECK_REAL"].Get<int>();
  m_check_virt_born = s["CHECK_VIRT_BORN"].Get<int>();
  m_virtual_only = s["VIRTUAL_ONLY"].Get<bool>();
  m_real_only = s["REAL_ONLY"].Get<bool>();
  m_use_model_alpha = s["USE_MODEL_ALPHA"].Get<bool>();
  m_kkmcAngles =  s["KKMC_ANG"].Get<int>();
  m_fixed_weight = s["FIXED_WEIGHT"].Get<double>();
  m_hardmin = s["HARD_MIN"].Get<double>();
  m_photonMass = s["PHOTON_MASS"].Get<double>();
  m_useceex = s["CEEX"].Get<int>();
  m_coll_real = s["Collinear_Real"].Get<bool>();
  m_resonace_max = s["CLUSTERING_THRESHOLD"].Get<double>();
  m_no_born = s["No_Born"].Get<int>();
  m_no_subtraction = s["No_Sub"].Get<int>();
  m_submode = s["Sub_Mode"].Get<submode::code>();
  m_tchannel = s["TChannel"].Get<int>();
  m_noflux = s["No_Flux"].Get<int>();
  m_flux_mode=s["Flux_Mode"].Get<int>();
  m_ifisub = s["IFI_Sub"].Get<int>();
  m_massless_sub = s["Massless_Sub"].Get<int>();
  m_check_real_sub = s["Check_Real_Sub"].Get<bool>();
  m_photon_split = s["PHOTON_SPLITTER_MODE"].ResetDefault().SetDefault(0).Get<bool>();
  m_int_nlo = s["Integrate_NLO"].SetDefault(0).Get<bool>();
  m_CalForm = false;
  m_realtool = false;
  //update when beamstrahlung is added
  m_isrinital=true;
  m_g = 0;
  m_gp = 0;

  if(m_use_model_alpha) m_alpha = s_model->ScalarConstant("alpha_QED");
  else m_alpha  = 1./s["1/ALPHAQED(0)"].SetDefault(137.03599976).Get<double>(); 
  if (m_use_model_alpha) m_rescale_alpha = 1;
  else m_rescale_alpha = m_alpha / s_model->ScalarConstant("alpha_QED");
  m_alpi = m_alpha/M_PI;

}

std::istream &YFS::operator>>(std::istream &str,submode::code &sub)
{
  std::string tag;
  str>>tag;
  sub=submode::local;
  if      (tag.find("Off")!=std::string::npos)    sub=submode::off;
  else if (tag.find("0")!=std::string::npos)      sub=submode::off;
  else if (tag.find("Local")!=std::string::npos)  sub=submode::local;
  else if (tag.find("1")!=std::string::npos)      sub=submode::local;
  else if (tag.find("Global")!=std::string::npos) sub=submode::global;
  else if (tag.find("2")!=std::string::npos)      sub=submode::global;
  return str;
}

std::istream &YFS::operator>>(std::istream &str, yfsmode::code &mode)
{
  std::string tag;
  str>>tag;
  mode=yfsmode::off;
  if      (tag.find("Off")!=std::string::npos)    mode=yfsmode::off;
  else if (tag.find("ISRFSR")!=std::string::npos) mode=yfsmode::isrfsr;
  else if (tag.find("ISR")!=std::string::npos)    mode=yfsmode::isr;
  else if (tag.find("FSR")!=std::string::npos)    mode=yfsmode::fsr;
  // else THROW(fatal_error, "Unknown YFS: MODE for Lepton Colliders")
  return str;
}

std::ostream &YFS::operator<<(std::ostream &str,const yfsmode::code &ym)
{
  if      (ym==yfsmode::off)    return str<<"Off";
  else if (ym==yfsmode::isr)    return str<<"ISR";
  else if (ym==yfsmode::isrfsr) return str<<"ISR+FSR";
  else if (ym==yfsmode::fsr)    return str<<"FSR";
  return str<<"unknown";
}

double YFS_Base::Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double YFS_Base::EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  // return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
  return m_alpha/(4*M_PI*M_PI)*(2*p1*p2/((p1*k)*(p2*k)));
}
