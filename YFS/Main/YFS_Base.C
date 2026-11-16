#include <strings.h>
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
  // Stub: identity, as returned when the form factor is off.
  return 1.;
}


void YFS_Base::RegisterDefaults(){
  m_s = sqr(rpa->gen.Ecms());
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  s["MODE"].ResetDefault().SetDefault(yfsmode::off);
  // The paper settings are the defaults, so a card needs only ACRAIC: 1.
  s["BETA"].SetDefault(1);
  s["VMAX"].SetDefault(0);
  // lower bound on |t| of the exchange lines, GeV^2 (0: off);
  // see YFS_Handler::ExchangeLinesPassTMin
  s["EXCHANGE_T_MIN"].SetDefault(0.);
  // lower bound on |t| of the exchange lines of the physical event (full beams
  // against outgoing legs after all radiation), GeV^2 (0: off). Pure phase-space veto.
  s["EXCHANGE_T_MIN_FINAL"].SetDefault(0.);
  // degrees between Born photons and the reduced beams (0: off)
  s["BORN_PHOTON_REDUCED_ANGLE"].SetDefault(0.);
  // lab-frame photon resolution for n-photon samples (0: off)
  s["RESOLVED_PHOTON_E"].SetDefault(0.);
  s["RESOLVED_PHOTON_ANGLE"].SetDefault(0.);
  // paper value; results are cutoff independent
  s["IR_CUTOFF"].ResetDefault().SetDefault(1e-5);
  s["DELTA"].SetDefault(1e-2);
  s["PHOTON_MAX"].SetDefault(100);
  s["LOOP_TOOL"].SetDefault(false);
  s["FSR_DEBUG"].SetDefault(false);
  s["ISR_DEBUG"].SetDefault(false);
  s["DEBUG_DIR_ISR"].SetDefault("ISR_DEBUG");
  s["DEBUG_DIR_FSR"].SetDefault("FSR_DEBUG");
  s["DEBUG_DIR_NLO"].SetDefault("YFS_NLO_Hist");
  s["COULOMB"].SetDefault(false);
  s["COULOMB_HIGHER_ORDER"].SetDefault(coulombho::none);
  // LOOP_SOFT_FLOOR (x = 2E/sqrt s): below it a photon set gets no multi-photon
  // one-loop evaluation (the massive (n+2)-point loop loses precision there)
  s["LOOP_SOFT_FLOOR"].SetDefault(1e-4);
  s["HIDE_PHOTONS"].SetDefault(true);
  s["FULL_FORM"].SetDefault(1);
  s["WW_FORM"].SetDefault(false);
  // W+W- -> 4f scheme: flat exponentiates the four external fermions (photons
  // softer than Gamma_W), pole exponentiates the pole expansion (production
  // stage plus a decay stage per W, photons harder than Gamma_W).
  s["WW_Scheme"].SetDefault(wwscheme::flat);
  // Put the W's exactly on M_W before building the pole-expansion dipoles.
  s["WW_OnShell"].SetDefault(false);
  // Pole expansion also drives photon emission (off: exponent only).
  // Only read when WW_Scheme is pole.
  s["WW_Pole_Emission"].SetDefault(true);
  s["WW_BETAT"].SetDefault(0.382);
  s["CHECK_MASS_REG"].SetDefault(0);
  s["CHECK_POLES"].SetDefault(0);
  s["CHECK_REAL"].SetDefault(0);
  s["CHECK_RV"].SetDefault(0);
  s["RV_Hard_Photon"].SetDefault(0);
  s["ACRAIC_Compare"].SetDefault(0);
  // diagnostic: per-event ingredients of ACRAIC and YFS.NLO_EW as named
  // columns YFS.DB_* (YFS_Handler::BuildNamedWeights); off = bit-identical
  s["DIFF_BUDGET_STREAMS"].SetDefault(false);
  s["VV_Approx_Uncertainty"].SetDefault(1.0);
  s["RV_TEST_PHOTON_X"].SetDefault(0.3);
  s["RV_TEST_PHOTON_THETA"].SetDefault(M_PI / 2.);
  s["RV_TEST_PHOTON_PHI"].SetDefault(0.);
  s["RV_SOFT_CUT"].SetDefault(1e-4);
  // skip the real-virtual of a photon collinear to a charged leg i,
  // min_i (p_i.k)/s below this (0: off); see NLO_Base::CalculateRealVirtual
  s["RV_COLLINEAR_CUT"].SetDefault(0.);
  s["RR_SOFT_CUT"].SetDefault(1e-5);
  s["RV_CANCEL_EPS"].SetDefault(0.);
  s["RV_CANCEL_HIST"].SetDefault(0);
  s["RV_ME_MAX_RATIO"].SetDefault(0.);
  s["CHECK_REAL_REAL"].SetDefault(0);
  s["CHECK_VIRT_BORN"].SetDefault(1);
  s["VIRTUAL_ONLY"].SetDefault(false);
  s["REAL_ONLY"].SetDefault(false);
  // photon emission with alpha(0), the hard process in the model scheme
  s["USE_MODEL_ALPHA"].SetDefault(false);
  s["KKMC_ANG"].SetDefault(0);
  s["WEIGHT_MODE"].SetDefault(wgt::full);
  s["PHOTON_MASS"].SetDefault(0.1);
  // ACRAIC is the nominal weight wherever YFS is on. The default
  // does NOT switch YFS on: only an explicit ACRAIC: true implies MODE (below).
  s["ACRAIC"].SetDefault(true);
  /*
    Which O(alpha) correction the event carries as its NOMINAL weight
    (tristate, or the old -1/0/1):
      auto (-1, default): follow ACRAIC - with ACRAIC on, ACRAIC (the
        amplitude-level ACRAIC weight) is the nominal and YFS.NLO (squared
        level) becomes the named column YFS.NLO_EW; with ACRAIC off, YFS.NLO as
        before.
      on (1) / off (0): force ACRAIC / YFS.NLO.
    Every other named column keeps its absolute weight either way.
  */
  s["ACRAIC_WEIGHT"].SetDefault(tristate::automatic);
  /*
    Where ACRAIC's virtual comes from. auto (default), per process with V in
    NLO_Part: the external loop provider when the process names a
    Loop_Generator, else ACRAIC's own (e+e- -> f fbar only; refused for other
    processes, YFS_Handler::MakeACRAIC). Without V in NLO_Part there is no
    virtual whatever this says. acraic / external / helicity / none force it.
  */
  s["ACRAIC_Virtual"].SetDefault(acraicvirt::automatic);
  s["Collinear_Real"].SetDefault(false);
  s["CLUSTERING_THRESHOLD"].SetDefault(10);
  // End-of-run counts of why CalculateFSR returned false or set a zero
  // weight, per dipole and failure kind (YFS_Handler::CountFSRFailure).
  s["FSR_FAILURE_STATS"].SetDefault(0);
  // Pole scheme: subtract the decay dipoles' virtual B from the loop
  // (Define_Dipoles::DecayVirtualSubtraction). false (default), unvalidated.
  s["WW_DECAY_VIRTUAL_SUB"].SetDefault(false);
  /*
    beta_0^(2) (Beta0_2.H, NLO_Beta02.C): YFS: BETA0_2 for YFS.NLO, ACRAIC:
    BETA0_2 for ACRAIC. Shared keys: BETA0_2_PAIRS (all|same_side),
    BETA0_2_CROSS (pair products in delta_f^2/2), BETA0_2_WEAK_PRODUCT
    (delta_f w), BETA0_2_ALPHA (alpha0|loop), BETA0_2_BOXES (two_loop boxes,
    Beta0_2_Boxes.H), BETA0_2_DIPHOTON (two_loop gamma gamma hard function,
    hep-ph/0201274) and BETA0_2_DIPHOTON_POWER_CUT (m_e^2/min(|t|,|u|) above
    which it is skipped). BETA0_2_CROSSED_VERTEX (high_energy|exact): the
    two-loop vertex of a line through a space-like 2 -> 2 (mu e -> mu e),
    Beta0_2_Vertex.H. BETA0_2_W_SQUARE (true): add (1/2) w^2, w = delta -
    delta_f, so that w is exponentiated with delta_f and delta^(2) = delta^2/2
    + the logarithmic two-loop terms (the exact vertex structure); needs
    BETA0_2_WEAK_PRODUCT. BETA0_2_MUE_BOXES (photonic|none, default photonic):
    the photonic two-loop remainder of l e -> l e (l = mu, tau) beyond the
    vertices (Beta0_2_Boxes.H), tabulated from the matrix element of
    arXiv:2212.06481; it checks the process and the table's range per event and
    does nothing elsewhere. BETA0_2 on retires the VV_EEX estimate.
  */
  s["BETA0_2"].SetDefault(beta02::none);
  s["BETA0_2_PAIRS"].SetDefault(beta02pairs::all);
  s["BETA0_2_CROSSED_VERTEX"].SetDefault(beta02vertex::exact);
  s["BETA0_2_CROSS"].SetDefault(true);
  s["BETA0_2_WEAK_PRODUCT"].SetDefault(true);
  s["BETA0_2_W_SQUARE"].SetDefault(true);
  s["BETA0_2_ALPHA"].SetDefault(beta02alpha::alpha0);
  s["BETA0_2_BOXES"].SetDefault(true);
  s["BETA0_2_MUE_BOXES"].SetDefault(beta02mueboxes::photonic);
  //! YFS: BETA0_2_HIGHER none|exponentiated (none): two_loop plus the higher orders
  //! of exp(delta_f)(1 + NLL), the LL delta_f^3/6 and beyond
  s["BETA0_2_HIGHER"].SetDefault(beta02higher::none);
  s["BETA0_2_DIPHOTON"].SetDefault(true);
  Settings::GetMainSettings()["ACRAIC"]["BETA0_2"].SetDefault(beta02::none);
  /*
    YFS: DYSON_VP none|resummed|complete (none): the Dyson resummation of
    the photon vacuum polarisation on the Born-level amplitude of both
    columns (YFS/Main/Dyson_VP.H): resummed beyond a full-EW provider's first
    order, complete for a provider without the photon self-energy.
  */
  s["DYSON_VP"].SetDefault(dysonvp::none);
  //! YFS: DYSON_VP_LEPTONS one_loop|two_loop (one_loop): the lepton loops in its Pi
  s["DYSON_VP_LEPTONS"].SetDefault(dysonvpleptons::one_loop);
  //! YFS: BETA0_2_VP_INSERTIONS none|vertex|vertex_box (none): lepton loops in the one-loop vertices, boxes
  s["BETA0_2_VP_INSERTIONS"].SetDefault(beta02vp::none);
  /*
    YFS: M1_SQUARE averaged|helicity - the |m_1|^2 of YFS.NLO's beta_0^(2);
    ACRAIC: EW_HIGHER_ORDER_HELICITY none|averaged|helicity - EW_HIGHER_ORDER
    (GRIFFIN NNLO - NLO) in ACRAIC's helicity virtual.
  */
  s["M1_SQUARE"].SetDefault(m1square::averaged);
  Settings::GetMainSettings()["ACRAIC"]["EW_HIGHER_ORDER_HELICITY"].SetDefault(ewhohel::none);
  s["NLO_Weight_Breakdown"].SetDefault(false);
  s["NLO_Columns"].SetDefault(false);
  s["Ladder_Weights"].SetDefault(false);
  s["Dump_Dipoles"].SetDefault(0);
  /*
    NLO diagnostics. Off by default, read through Settings rather than the
    environment for the same reason as the ACRAIC ones: a run has to be
    reproducible from its card alone.
  */
  /*
    REAL_FSR_MAP, the (n+1)-body point of a final-state photon at n >= 2:
      dipole (2, default): the photon's own dipole recoils, other resonant
        pairs keep their momenta (NLO_Base::MapMomentaFSRDipole); else
        rescale_all
      rescale_all (1, also 4): all final legs rescaled together, Q fixed
      pre_emission (3): rescale_all from the pre-emission legs
      rest_frame (0): the rest-frame map (puts sqrt(s_j) above sqrt(s))
  */
  s["REAL_FSR_MAP"].SetDefault(realfsrmap::dipole);
  /*
    REAL_MAP, the (n+1)-body point of an initial-state photon at n >= 2:
      scaled (2, default): scaled reduction (KKMC convention);
      invariant (3): scaled, with (final + k)^2 kept;
      beam_axis (1): beam-axis reduction;
      rest_frame (0): rest-frame map (legacy).
  */
  s["REAL_MAP"].SetDefault(realmap::scaled);
  s["SOFT_SCAN"].SetDefault(0.0);      // rel. deviation that triggers SoftScan, @@@ SCAN
  s["No_Born"].SetDefault(false);
  s["No_Sub"].SetDefault(false);
  s["Sub_Mode"].SetDefault(submode::global);
  s["No_Flux"].SetDefault(false);
  /*
    Flux_Mode, the initial-state flux the real (and RV, RR) carries:
      mapped (1, default): the NLO dipoles' flux at the mapped photon;
      event (0): the generation dipoles' flux at the event's photon;
      average (2): the mean of the two (single real and RV; the double real
        treats it as event).
  */
  s["Flux_Mode"].SetDefault(fluxmode::mapped);
  /*
    YFS: IFI (default true): initial-final interference at all
  */
  s["IFI"].SetDefault(true);
  s["IFI_Sub"].SetDefault(true);
  s["IFI_Real"].SetDefault(true);
  s["IFI_Omega"].SetDefault(0.);
  /*
    YFS: IF_VIRTUAL tbvirt_exact_ir|exact|tbvirt|s_channel (tbvirt_exact_ir): the
    virtual B of an initial-final dipole, the same in the exponent and the
    subtraction (YFS_Options.H). The old switch IF_VIRTUAL_TCHANNEL, if set
    explicitly, still selects tbvirt (true) or s_channel (false).
  */
  s["IF_VIRTUAL"].SetDefault(ifvirtual::tbvirt_exact_ir);
  s["IF_VIRTUAL_TCHANNEL"].SetDefault(true);

  s["IFI_RClip"].SetDefault(0.);
  s["Massless_Sub"].SetDefault(false);
  s["Check_Real_Sub"].SetDefault(0);
  s["Check_RR_Sub"].SetDefault(0);
  s["Integrate_NLO"].SetDefault(false);
  s["Collinear_Virtual"].SetDefault(false);
  s["Virtual_Sub"].SetDefault(true);
  s["Dim_Reg"].SetDefault(true);
  s["IR_SCALE"].SetDefault(100);
  s["NLO_CUTS"].SetDefault(false);
  s["Fixed_Order"].SetDefault("Full");
  s["SKIP_NEG_WEIGHTS"].SetDefault(false);
  s["NLO_FSR_PHOTONS"].SetDefault(true);

  s["NLO_FSR_FROM_EVENT"].SetDefault(true);
}


static void RefuseOldAcraicKeys()
{
  const std::string old("CE" "EX");
  Settings &main(Settings::GetMainSettings());
  std::vector<std::string> bad;
  if (!main[old].GetKeys().empty() || main[old].IsSetExplicitly())
    bad.push_back(old + ":");
  for (const std::string &k : main["YFS"].GetKeys())
    if (k.compare(0, old.size(), old) == 0) bad.push_back("YFS: " + k);
  if (bad.empty()) return;
  std::string list;
  for (const std::string &b : bad) list += (list.empty() ? "" : ", ") + b;
  THROW(fatal_error, "the " + old + " settings were renamed to ACRAIC on "
        "2026-10-06; this card still has: " + list + ". Write ACRAIC in "
        "place of " + old + " in these keys (the block, YFS: ACRAIC, "
        "YFS: ACRAIC_*).");
}

void YFS_Base::RegisterSettings(){
  RefuseOldAcraicKeys();
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  m_betaorder = s["BETA"].Get<int>();
  m_mode = s["MODE"].Get<yfsmode::code>();
  // ACRAIC: 1 alone switches YFS on: ISRFSR here, narrowed to ISR in
  // YFS_Handler::SetFlavours when no final-state particle is charged. An
  // explicit MODE always wins, including MODE: off.
  if (!s["MODE"].IsSetExplicitly() && m_mode == yfsmode::off &&
      s["ACRAIC"].IsSetExplicitly() && s["ACRAIC"].Get<bool>()) { m_mode = yfsmode::isrfsr; m_mode_from_acraic = true; }
  m_isrcut   = s["IR_CUTOFF"].Get<double>();
  m_isrcut = m_isrcut/sqrt(m_s); // dimensionless units
  m_vmax = s["VMAX"].Get<double>();
  m_vmax = 1.-sqr(m_vmax)/m_s;
  m_exchange_tmin = s["EXCHANGE_T_MIN"].Get<double>();
  m_exchange_tmin_final = s["EXCHANGE_T_MIN_FINAL"].Get<double>();
  m_bornphoton_reduced_angle = s["BORN_PHOTON_REDUCED_ANGLE"].Get<double>();
  m_resolved_photon_e = s["RESOLVED_PHOTON_E"].Get<double>();
  m_resolved_photon_angle = s["RESOLVED_PHOTON_ANGLE"].Get<double>();
  m_looptool  = s["LOOP_TOOL"].Get<bool>();
  m_debugDIR_ISR = s["DEBUG_DIR_ISR"].Get<std::string>();
  m_debugDIR_FSR = s["DEBUG_DIR_FSR"].Get<std::string>();
  m_debugDIR_NLO = s["DEBUG_DIR_NLO"].Get<std::string>();
  m_fsr_debug = s["FSR_DEBUG"].Get<bool>();
  m_isr_debug = s["ISR_DEBUG"].Get<bool>();
  m_deltacut = s["DELTA"].Get<double>()*m_isrcut;
  m_coulomb = s["COULOMB"].Get<bool>();
  m_coulombho = s["COULOMB_HIGHER_ORDER"].Get<coulombho::code>();
  if (m_coulomb && m_coulombho != coulombho::none)
    THROW(fatal_error, "YFS: COULOMB (pole scheme, the whole Coulomb factor) and "
          "COULOMB_HIGHER_ORDER (flat scheme, beyond the one-loop virtual) double count");
  m_hidephotons = s["HIDE_PHOTONS"].Get<bool>();
  m_fullform = s["FULL_FORM"].Get<int>();
  m_formWW = s["WW_FORM"].Get<bool>();
  m_wwscheme = s["WW_Scheme"].Get<wwscheme::code>();
  m_wwonshell = s["WW_OnShell"].Get<bool>();
  m_wwpoleemission = s["WW_Pole_Emission"].Get<bool>();
  m_betatWW = s["WW_BETAT"].Get<double>();
  m_check_mass_reg = s["CHECK_MASS_REG"].Get<int>();
  m_check_poles = s["CHECK_POLES"].Get<int>();
  m_check_real = s["CHECK_REAL"].Get<int>();
  m_check_rr = s["CHECK_REAL_REAL"].Get<int>();
  m_check_rv = s["CHECK_RV"].Get<int>();
  m_rv_hard_photon = s["RV_Hard_Photon"].Get<int>();
  m_acraic_compare = s["ACRAIC_Compare"].Get<int>();
  m_vv_approx_unc = s["VV_Approx_Uncertainty"].Get<double>();
  m_rv_test_x = s["RV_TEST_PHOTON_X"].Get<double>();
  m_rv_test_theta = s["RV_TEST_PHOTON_THETA"].Get<double>();
  m_rv_test_phi = s["RV_TEST_PHOTON_PHI"].Get<double>();
  // Dimensionless: photon energy fraction E/sqrt(s). <=0 disables the guard.
  m_rv_soft_cut = s["RV_SOFT_CUT"].Get<double>();
  m_rv_collinear_cut = s["RV_COLLINEAR_CUT"].Get<double>();
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
  m_photonMass = s["PHOTON_MASS"].Get<double>();
  m_useacraic = m_mode != yfsmode::off && s["ACRAIC"].Get<bool>();
  const tristate::code acraicweight(s["ACRAIC_WEIGHT"].Get<tristate::code>());
  m_acraic_weight = (acraicweight == tristate::automatic ? m_useacraic
                                                     : acraicweight == tristate::on);
  m_acraicvirtsrc = s["ACRAIC_Virtual"].Get<acraicvirt::code>();
  m_coll_real = s["Collinear_Real"].Get<bool>();
  m_resonace_max = s["CLUSTERING_THRESHOLD"].Get<double>();
  m_nlo_weight_breakdown = s["NLO_Weight_Breakdown"].Get<bool>();
  m_nlo_columns = s["NLO_Columns"].Get<bool>();
  m_ladder_weights = s["Ladder_Weights"].Get<bool>();
  m_dump_dipoles = s["Dump_Dipoles"].Get<int>();
  m_no_born_setting = s["No_Born"].Get<bool>();
  m_no_born = m_no_born_setting;
  m_no_subtraction = s["No_Sub"].Get<bool>();
  m_submode = s["Sub_Mode"].Get<submode::code>();
  m_noflux = s["No_Flux"].Get<bool>();
  m_flux_mode = s["Flux_Mode"].Get<fluxmode::code>();
  m_ifi = s["IFI"].Get<bool>();
  m_ifisub = m_ifi && s["IFI_Sub"].Get<bool>();
  m_ifireal = m_ifi && s["IFI_Real"].Get<bool>();
  m_ifiomega = s["IFI_Omega"].Get<double>();
  m_ifvirtual = s["IF_VIRTUAL"].Get<ifvirtual::code>();
  if (s["IF_VIRTUAL_TCHANNEL"].IsSetExplicitly())
    m_ifvirtual = s["IF_VIRTUAL_TCHANNEL"].Get<bool>() ? ifvirtual::tbvirt : ifvirtual::s_channel;
  // Defaults to FSR::Initialize()'s m_Emin (same expression) when IFI_Real is
  // on; an explicit IFI_Omega is honoured either way.
  m_ifiomega_is_emin = (m_ifireal && m_ifiomega <= 0.);
  if (m_ifiomega_is_emin) m_ifiomega = 0.5*sqrt(m_s)*m_isrcut;
  m_ifi_rclip = s["IFI_RClip"].Get<double>();
  m_ifi_clipped = 0;
  m_massless_sub = s["Massless_Sub"].Get<bool>();
  // 0 = off, 1 = one-shot energy-scan sub check (CheckReal[Real]Sub, exits),
  // 2 = accumulating angle/energy scatter (RecordSubScatter, no exit)
  m_check_real_sub = s["Check_Real_Sub"].Get<int>();
  m_check_rr_sub = s["Check_RR_Sub"].Get<int>();
  m_photon_split = s["PHOTON_SPLITTER_MODE"].ResetDefault().SetDefault(0).Get<bool>();
  m_int_nlo = s["Integrate_NLO"].Get<bool>();
  m_eex_virt = s["Collinear_Virtual"].Get<bool>();
  m_virt_sub = s["Virtual_Sub"].Get<bool>();
  m_dim_reg = s["Dim_Reg"].Get<bool>();
  m_irscale = s["IR_SCALE"].Get<double>();
  m_nlocuts = s["NLO_CUTS"].Get<bool>();
  // Fixed_Order is removed: use ACRAIC: FIXED_ORDER_STREAM
  if (s["Fixed_Order"].Get<std::string>() != "Full")
    THROW(not_implemented, "YFS: Fixed_Order is removed. For fixed-order results"
          " use the named weights YFS.ACRAIC_FixedOrder (O(alpha)) and"
          " YFS.ACRAIC_FixedOrder2 (O(alpha^2), ACRAIC: ORDER 2) with ACRAIC:"
          " {FIXED_ORDER_STREAM: true}.");
  m_skipNegWeights = s["SKIP_NEG_WEIGHTS"].Get<bool>();
  m_nlo_fsr_photons = s["NLO_FSR_PHOTONS"].Get<bool>();
  m_nlo_fsr_from_event = s["NLO_FSR_FROM_EVENT"].Get<bool>();
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

std::istream &YFS::operator>>(std::istream &str, acraicvirt::code &sc)
{
  std::string tag;
  str>>tag;
  if      (tag=="acraic"     || tag=="ACRAIC"     || tag=="0") sc=acraicvirt::acraic;
  else if (tag=="external" || tag=="External" || tag=="1") sc=acraicvirt::external;
  else if (tag=="helicity" || tag=="Helicity" || tag=="2") sc=acraicvirt::helicity;
  else if (tag=="none"     || tag=="None"     || tag=="3") sc=acraicvirt::none;
  else if (tag=="auto"     || tag=="Auto"     || tag=="4") sc=acraicvirt::automatic;
  else THROW(fatal_error, "Unknown YFS ACRAIC_Virtual '"+tag
                          +"'; use acraic, external, helicity or none.");
  return str;
}

std::ostream &YFS::operator<<(std::ostream &str, const acraicvirt::code &sc)
{
  switch (sc) {
  case acraicvirt::acraic:     return str<<"acraic";
  case acraicvirt::external: return str<<"external";
  case acraicvirt::helicity: return str<<"helicity";
  case acraicvirt::none:     return str<<"none";
  case acraicvirt::automatic: return str<<"auto";
  }
  return str<<"unknown";
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

// YFS: MODE, matched whole and case-insensitively, by name only. Photons++
// reads the same key (None, Soft, Full or 0, 1, 2), so any other value,
// numbers included, leaves YFS off rather than being an error.
static const std::vector<YFS::Option_Name> s_yfsmode{
  {"Off", yfsmode::off}, {"ISR", yfsmode::isr}, {"ISRFSR", yfsmode::isrfsr},
  {"ISR+FSR", yfsmode::isrfsr}, {"FSR", yfsmode::fsr}};

std::istream &YFS::operator>>(std::istream &str, yfsmode::code &mode)
{
  std::string tag;
  str>>tag;
  mode=yfsmode::off;
  for (const YFS::Option_Name &n : s_yfsmode)
    if (strcasecmp(tag.c_str(), n.name) == 0)
      mode=static_cast<yfsmode::code>(n.value);
  return str;
}

std::ostream &YFS::operator<<(std::ostream &str,const yfsmode::code &ym)
{
  return str<<OptionName(ym, s_yfsmode);
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


double YFS_Base::Eikonal(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return -m_alpha / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double YFS_Base::EikonalMassless(const Vec4D &k, const Vec4D &p1, const Vec4D &p2) {
  return m_alpha/(4*M_PI*M_PI)*(2*p1*p2/((p1*k)*(p2*k)));
}
