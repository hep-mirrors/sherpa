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


double YFS_Base::ExternalFormFactor(const Vec4D_Vector &p,
                                    const Flavour_Vector &fl,
                                    const size_t &nin)
{
  // Stub. The pion form factor lands in its own merge request; until then this
  // is the identity, which is what the real implementation returns whenever the
  // form factor is switched off.
  return 1.;
}


void YFS_Base::RegisterDefaults(){
  m_s = sqr(rpa->gen.Ecms());
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  s["MODE"].ResetDefault().SetDefault(yfsmode::off);
  s["BETA"].SetDefault(2);
  s["VMAX"].SetDefault(0);
  s["IR_CUTOFF"].ResetDefault().SetDefault(1e-6);
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
  s["FULL_FORM"].SetDefault(1);
  s["WW_FORM"].SetDefault(0);
  // Which YFS scheme to use for a W+W- -> 4f final state.
  //   flat : exponentiate the four external fermions. Correct for photons
  //          softer than Gamma_W, where the W lives too briefly to be resolved.
  //   pole : exponentiate the pole expansion -- a production stage with the
  //          W's as charged emitters and a decay stage per W. Correct for
  //          photons harder than Gamma_W. This is the YFSWW3 picture.
  s["WW_Scheme"].SetDefault(wwscheme::flat);
  // Put the W's exactly on M_W before building the pole-expansion dipoles,
  // rather than letting them carry the invariant mass of their own decay
  // products. The leading-pole approximation strictly wants on shell; the
  // difference is an LPA ambiguity worth measuring rather than choosing.
  s["WW_OnShell"].SetDefault(0);
  // Whether the pole expansion also drives photon emission, or only the form
  // factor. With it off the production-stage W pair does not radiate, so the
  // pole scheme changes the exponent alone and the two halves of the scheme
  // can be measured separately. Only read when WW_Scheme is pole.
  s["WW_Pole_Emission"].SetDefault(1);
  s["WW_BETAT"].SetDefault(0.382);
  s["CHECK_MASS_REG"].SetDefault(0);
  s["CHECK_POLES"].SetDefault(0);
  s["CHECK_REAL"].SetDefault(0);
  s["CHECK_RV"].SetDefault(0);
  s["RV_Hard_Photon"].SetDefault(0);
  s["CEEX_Compare"].SetDefault(0);
  s["VV_Approx_Uncertainty"].SetDefault(1.0);
  s["RV_TEST_PHOTON_X"].SetDefault(0.3);
  s["RV_TEST_PHOTON_THETA"].SetDefault(M_PI / 2.);
  s["RV_TEST_PHOTON_PHI"].SetDefault(0.);
  s["RV_SOFT_CUT"].SetDefault(0.);
  s["RR_SOFT_CUT"].SetDefault(0.);
  s["RV_CANCEL_EPS"].SetDefault(0.);
  s["RV_CANCEL_HIST"].SetDefault(0);
  s["RV_ME_MAX_RATIO"].SetDefault(0.);
  s["CHECK_REAL_REAL"].SetDefault(0);
  s["CHECK_VIRT_BORN"].SetDefault(1);
  s["VIRTUAL_ONLY"].SetDefault(0);
  s["REAL_ONLY"].SetDefault(0);
  s["USE_MODEL_ALPHA"].SetDefault(1);
  s["KKMC_ANG"].SetDefault(0);
  s["WEIGHT_MODE"].SetDefault(wgt::full);
  s["HARD_MIN"].SetDefault(0.);
  s["PHOTON_MASS"].SetDefault(0.1);
  s["CEEX"].SetDefault(0);
  s["Collinear_Real"].SetDefault(0);
  s["CLUSTERING_THRESHOLD"].SetDefault(10);
  s["TChannel"].SetDefault(0);
  s["NLO_Weight_Breakdown"].SetDefault(0);
  s["Ladder_Weights"].SetDefault(0);
  s["Dump_Dipoles"].SetDefault(0);
  s["CHECK_INVARIANTS"].SetDefault(0);
  s["No_Born"].SetDefault(0);
  s["No_Sub"].SetDefault(0);
  s["Sub_Mode"].SetDefault(submode::global);
  s["No_Flux"].SetDefault(0);
  s["Flux_Mode"].SetDefault(1);
  s["IFI_Sub"].SetDefault(1);
  // Emission-side initial-final interference. Off by default: with it off the
  // IF form factor holds the whole soft integral (omega = sqrt(s)/2) and is
  // cutoff-independent, which is the safe inclusive approximation. With it on,
  // Define_Dipoles::RealIFWeight() reweights every generated photon by the IF
  // radiation function and the form factor's cutoff drops to IFI_Omega, so the
  // two must be switched together - see Define_Dipoles::IFIOmega().
  s["IFI_Real"].SetDefault(0);
  // Soft cutoff shared by the IF form factor and the photon reweighting, in
  // GeV. Only read when IFI_Real is on.
  //
  // ISR and FSR generate against different cutoffs (v > IR_CUTOFF/sqrt(s) in
  // the CMS vs FSR_CUT in the dipole rest frame), but that is not a mismatch
  // to work around: FSR::YFS_FORM() carries the Piatek translation
  // m_DelYFS = Btilda(dipole, m_Emin) - Btilda(Q frame, m_EminQ) (FSR.C:488),
  // exactly as KKMC does in KKarFin/KKceex.cxx:288-296, which moves the FSR
  // bookkeeping onto the single scale m_Emin. KKMC then uses that one Emin for
  // Yisr, Yfsr AND Yint alike, so the IF cutoff belongs on the same scale.
  //
  // 0 therefore means "use FSR::Initialize()'s m_Emin", which is IR_CUTOFF/2.
  s["IFI_Omega"].SetDefault(0.);
  // Diagnostic clamp on the per-photon IF reweight, OFF by default (<= 0).
  // RealIFWeight cancels against beta_1 exactly, so clamping is not a safety
  // net - it injects a residue exactly where it fires. Only for bisecting.
  s["IFI_RClip"].SetDefault(0.);
  s["Massless_Sub"].SetDefault(0);
  s["Check_Real_Sub"].SetDefault(0);
  s["Check_RR_Sub"].SetDefault(0);
  s["Integrate_NLO"].SetDefault(1);
  s["Collinear_Virtual"].SetDefault(0);
  s["Virtual_Sub"].SetDefault(1);
  s["Dim_Reg"].SetDefault(1);
  s["IR_SCALE"].SetDefault(100);
  s["NLO_CUTS"].SetDefault(false);
  s["Fixed_Order"].SetDefault(fixed_order::full);
  s["SKIP_NEG_WEIGHTS"].SetDefault(false);
  s["NLO_FSR_PHOTONS"].SetDefault(true);
  s["MIN_PHOTON"].SetDefault<int>(-1);
  s["FB_Analysis"].SetDefault(false);
  s["FB_Analysis_KF"].SetDefault<int>(0);
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
  m_wwscheme = s["WW_Scheme"].Get<wwscheme::code>();
  m_wwonshell = s["WW_OnShell"].Get<int>();
  m_wwpoleemission = s["WW_Pole_Emission"].Get<int>();
  m_betatWW = s["WW_BETAT"].Get<double>();
  m_check_mass_reg = s["CHECK_MASS_REG"].Get<int>();
  m_check_poles = s["CHECK_POLES"].Get<int>();
  m_check_real = s["CHECK_REAL"].Get<int>();
  m_check_rv = s["CHECK_RV"].Get<int>();
  m_rv_hard_photon = s["RV_Hard_Photon"].Get<int>();
  m_ceex_compare = s["CEEX_Compare"].Get<int>();
  // Parse once. 16 tags is comfortably above any YFS multiplicity (2 -> n + a
  // few photons); unused ones cost nothing.
  m_vv_approx_unc = s["VV_Approx_Uncertainty"].Get<double>();
  m_rv_test_x = s["RV_TEST_PHOTON_X"].Get<double>();
  m_rv_test_theta = s["RV_TEST_PHOTON_THETA"].Get<double>();
  m_rv_test_phi = s["RV_TEST_PHOTON_PHI"].Get<double>();
  // Dimensionless: photon energy fraction E/sqrt(s). <=0 disables the guard.
  m_rv_soft_cut = s["RV_SOFT_CUT"].Get<double>();
  // Same, applied to each photon of the double-real pair. <=0 disables.
  m_rr_soft_cut = s["RR_SOFT_CUT"].Get<double>();
  // Dimensionless relative-cancellation threshold for the RV beta_1^1. <=0
  // disables. Recommended O(1e-10..1e-12); validate with RV_CANCEL_HIST.
  m_rv_cancel_eps = s["RV_CANCEL_EPS"].Get<double>();
  m_rv_cancel_hist = s["RV_CANCEL_HIST"].Get<int>();
  // Skip RV photons whose loop ME exceeds max(|rvsub|,|aB|) by this factor
  // (unstable one-loop provider). Recommended O(1e6); <=0 disables.
  m_rv_me_max_ratio = s["RV_ME_MAX_RATIO"].Get<double>();
  m_check_virt_born = s["CHECK_VIRT_BORN"].Get<int>();
  m_virtual_only = s["VIRTUAL_ONLY"].Get<bool>();
  m_real_only = s["REAL_ONLY"].Get<bool>();
  m_use_model_alpha = s["USE_MODEL_ALPHA"].Get<bool>();
  m_kkmcAngles =  s["KKMC_ANG"].Get<int>();
  m_fixed_weight = s["WEIGHT_MODE"].Get<wgt::code>();
  m_hardmin = s["HARD_MIN"].Get<double>();
  m_photonMass = s["PHOTON_MASS"].Get<double>();
  m_useceex = s["CEEX"].Get<int>();
  m_coll_real = s["Collinear_Real"].Get<bool>();
  m_resonace_max = s["CLUSTERING_THRESHOLD"].Get<double>();
  m_nlo_weight_breakdown = s["NLO_Weight_Breakdown"].Get<int>();
  m_ladder_weights = s["Ladder_Weights"].Get<int>();
  m_dump_dipoles = s["Dump_Dipoles"].Get<int>();
  m_check_invariants = s["CHECK_INVARIANTS"].Get<bool>();
  m_no_born_setting = s["No_Born"].Get<int>();
  m_no_born = m_no_born_setting;
  m_no_subtraction = s["No_Sub"].Get<int>();
  m_submode = s["Sub_Mode"].Get<submode::code>();
  m_tchannel = s["TChannel"].Get<int>();
  m_noflux = s["No_Flux"].Get<int>();
  m_flux_mode=s["Flux_Mode"].Get<int>();
  m_ifisub = s["IFI_Sub"].Get<int>();
  m_ifireal = s["IFI_Real"].Get<int>();
  m_ifiomega = s["IFI_Omega"].Get<double>();
  // Default to FSR::Initialize()'s m_Emin, the scale the Piatek term m_DelYFS
  // translates the FSR bookkeeping onto. Kept as the same expression, not a
  // rederived one, so the two cannot drift apart.
  // Only defaulted when IFI_Real is on; an explicit IFI_Omega is honoured
  // either way so the exponent and the restoration can be varied separately.
  if (m_ifireal && m_ifiomega <= 0.) m_ifiomega = 0.5*sqrt(m_s)*m_isrcut;
  m_ifi_rclip = s["IFI_RClip"].Get<double>();
  m_ifi_clipped = 0;
  m_massless_sub = s["Massless_Sub"].Get<int>();
  // 0 = off, 1 = one-shot energy-scan sub check (CheckReal[Real]Sub, exits),
  // 2 = accumulating angle/energy scatter (RecordSubScatter, no exit)
  m_check_real_sub = s["Check_Real_Sub"].Get<int>();
  m_check_rr_sub = s["Check_RR_Sub"].Get<int>();
  m_photon_split = s["PHOTON_SPLITTER_MODE"].ResetDefault().SetDefault(0).Get<bool>();
  m_int_nlo = s["Integrate_NLO"].Get<bool>();
  m_eex_virt = s["Collinear_Virtual"].Get<int>();
  m_virt_sub = s["Virtual_Sub"].Get<int>();
  m_dim_reg = s["Dim_Reg"].Get<bool>();
  m_irscale = s["IR_SCALE"].Get<double>();
  m_nlocuts = s["NLO_CUTS"].Get<bool>();
  m_fixedOrder = s["Fixed_Order"].Get<fixed_order::code>();
  m_skipNegWeights = s["SKIP_NEG_WEIGHTS"].Get<bool>();
  m_nlo_fsr_photons = s["NLO_FSR_PHOTONS"].Get<bool>();
  m_mingammaN = s["MIN_PHOTON"].Get<int>();
  m_fb_analysis = s["FB_Analysis"].Get<bool>();
  m_fb_kf = s["FB_Analysis_KF"].Get<int>();
  m_CalForm = false;
  m_realtool = false;
  //update when beamstrahlung is added
  m_isrinital=true;
  m_g = 0;
  m_gp = 0;
  m_failcut=false;
  double alpha0 = 1./s["1/ALPHAQED(0)"].SetDefault(137.03599976).Get<double>();
  if(m_use_model_alpha) m_alpha = s_model->ScalarConstant("alpha_QED");
  else m_alpha  = alpha0;
  if (m_use_model_alpha) m_rescale_alpha = 1.;//m_rescale_alpha = alpha0/m_alpha;
  else m_rescale_alpha = m_alpha / s_model->ScalarConstant("alpha_QED");
  m_alpi = m_alpha/M_PI;
}

std::istream &YFS::operator>>(std::istream &str, wwscheme::code &sc)
{
  std::string tag;
  str>>tag;
  sc=wwscheme::flat;
  if      (tag.find("Pole")!=std::string::npos) sc=wwscheme::pole;
  else if (tag.find("pole")!=std::string::npos) sc=wwscheme::pole;
  else if (tag.find("1")!=std::string::npos)    sc=wwscheme::pole;
  else if (tag.find("Flat")!=std::string::npos) sc=wwscheme::flat;
  else if (tag.find("flat")!=std::string::npos) sc=wwscheme::flat;
  else if (tag.find("0")!=std::string::npos)    sc=wwscheme::flat;
  else THROW(fatal_error, "Unknown YFS WW_Scheme '"+tag+"'; use flat or pole.");
  return str;
}

std::ostream &YFS::operator<<(std::ostream &str, const wwscheme::code &sc)
{
  return str<<(sc==wwscheme::pole?"pole":"flat");
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

std::ostream &YFS::operator<<(std::ostream &str,const wgt::code &wm)
{
  if      (wm==wgt::off)    return str<<"Off";
  else if (wm==wgt::full)   return str<<"Full";
  else if (wm==wgt::mass)   return str<<"Mass";
  else if (wm==wgt::hide)   return str<<"Hidden";
  else if (wm==wgt::jacob)  return str<<"Jacobian";
  return str<<"unknown";
}

std::istream &YFS::operator>>(std::istream &str, wgt::code &mode)
{
  std::string tag;
  str>>tag;
  // mode=wgt::off;
  if      (tag.find("Off")!=std::string::npos)    mode=wgt::off;
  else if (tag.find("Full")!=std::string::npos)   mode=wgt::full;
  else if (tag.find("Mass")!=std::string::npos)   mode=wgt::mass;
  else if (tag.find("Hidden")!=std::string::npos) mode=wgt::hide;
  else if (tag.find("Jacobian")!=std::string::npos) mode=wgt::jacob;
  else THROW(fatal_error, "Unknown YFS: WEIGHT_MODE")
  return str;
}

std::istream &YFS::operator>>(std::istream &str, fixed_order::code &mode)
{
  std::string tag;
  str>>tag;
  // mode=wgt::off;
  if      (tag.find("Off")!=std::string::npos)    mode=fixed_order::off;
  else if (tag.find("Full")!=std::string::npos)   mode=fixed_order::full;
  else if (tag.find("NLO")!=std::string::npos) mode=fixed_order::nlo;
  else if (tag.find("NNLO")!=std::string::npos) mode=fixed_order::nnlo;
  else THROW(fatal_error, "Unknown YFS: Fixed_Order")
  return str;
}

std::ostream &YFS::operator<<(std::ostream &str,const fixed_order::code &wm)
{
  if      (wm==fixed_order::off)    return str<<"Off";
  else if (wm==fixed_order::full)   return str<<"Full";
  else if (wm==fixed_order::lo)   return str<<"LO";
  else if (wm==fixed_order::nlo)   return str<<"NLO";
  else if (wm==fixed_order::nnlo)  return str<<"NNLO";
  return str<<"unknown";
}

double YFS_Base::Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double YFS_Base::EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  // return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
  return m_alpha/(4*M_PI*M_PI)*(2*p1*p2/((p1*k)*(p2*k)));
}
