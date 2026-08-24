#include "METOOLS/HadronCurrents/FormFactors/FF_0_PP.H"
#include "METOOLS/HadronCurrents/FormFactors/RChL_Functions.H"
#include "METOOLS/HadronCurrents/FormFactors/Kstar_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/A0_Decays.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


void FF_0_PP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi_plus)     m_mode = FF_0_PP_mode::pipi_plus;
  // K K-bar (tau -> K0bar K- nu): as with Kpi below, real events carry
  // K_S/K_L (the observable strangeness-eigenstate mixtures), not the
  // bare (unobservable) K^0/K0bar flavour eigenstate FixMode() used to
  // check exclusively - this was silently mis-tagging every physical
  // K_S/K_L event as "unknown" (mode 999), falling all the way through
  // to a flat constant form factor with NO resonance structure at all.
  // Fixed here; m_isKSKL flags the extra 1/sqrt(2) projection factor
  // needed in FixParameters() (see the member comment in FF_0_PP.H).
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_K ||
	    m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
	    m_flavs[m_pi[0]].Kfcode()==kf_K_L) &&
	   m_flavs[m_pi[1]].Kfcode()==kf_K_plus) {
    m_mode   = FF_0_PP_mode::KK_plus;
    m_isKSKL = (m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
		m_flavs[m_pi[0]].Kfcode()==kf_K_L);
  }
  // K-pi (Cabibbo-suppressed, tau -> K^- pi^0 nu / K0bar pi^- nu):
  // Kpi_plus mode already has working KS branches (vector T_K*, Sec.II
  // of Finkemeier-Mirkes hep-ph/9503474 = FM95) further down in this
  // file, but FixMode() never actually assigned this mode - it was
  // dead code, unreachable for any real K-pi decay. Fixed here.
  // Both K^-pi^0 and K0bar pi^- reduce to the SAME 1/sqrt(2) x
  // T_K*(Q^2) structure (FM95 Eq.5), so one mode value covers both
  // charge/neutral assignments; order (K,pi) or (pi,K) is accepted.
  // Same K_S/K_L caveat as KK_plus above applies here too.
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_K_plus ||
	    m_flavs[m_pi[0]].Kfcode()==kf_K ||
	    m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
	    m_flavs[m_pi[0]].Kfcode()==kf_K_L) &&
	   (m_flavs[m_pi[1]].Kfcode()==kf_pi_plus ||
	    m_flavs[m_pi[1]].Kfcode()==kf_pi)) {
    m_mode   = FF_0_PP_mode::Kpi_plus;
    m_isKSKL = (m_flavs[m_pi[0]].Kfcode()==kf_K_S ||
		m_flavs[m_pi[0]].Kfcode()==kf_K_L);
  }
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_pi_plus ||
	    m_flavs[m_pi[0]].Kfcode()==kf_pi) &&
	   (m_flavs[m_pi[1]].Kfcode()==kf_K_plus ||
	    m_flavs[m_pi[1]].Kfcode()==kf_K ||
	    m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
	    m_flavs[m_pi[1]].Kfcode()==kf_K_L)) {
    m_mode   = FF_0_PP_mode::Kpi_plus;
    m_isKSKL = (m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
		m_flavs[m_pi[1]].Kfcode()==kf_K_L);
  }
  // K-eta (tau -> K^- eta nu): per the cross-check note
  // (tau_two_meson_currents_KS_RChiT.tex), there is NO genuine original
  // Kuhn-Santamaria Keta model - what follows is a "KS-like" extension
  // (Eq.(eq:KS-Keta) of that note), built the same way as Kpi_plus but
  // with its own beta_K* and an eta-eta' mixing normalisation. See
  // Fplus_0_PiZeroPiPlus's Keta_plus branch below for details/flags.
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta) ||
	   (m_flavs[m_pi[0]].Kfcode()==kf_eta &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus))  m_mode = FF_0_PP_mode::Keta_plus;
  // K-eta' (tau -> K^- eta' nu): per the cross-check note Sec.5.7,
  // no KS-like model is given explicitly (only Keta gets Eq.(eq:KS-
  // Keta)), but the note's own Eq.(eq:Ketap-vector-relation) states
  // the SAME strange vector current controls Kpi, Keta and Keta' -
  // reusing Construct_Keta()'s K* tower here on that basis. The
  // genuine, dedicated treatment is "Level 1: bare RChiPT"
  // (Eq.(eq:Ketap-vector-twores)), identical in form to Kpi/Keta - see
  // Fplus_0_PiZeroPiPlus's Construct_Keta()/FF_RChiPT() for details.
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta_prime_958) ||
	   (m_flavs[m_pi[0]].Kfcode()==kf_eta_prime_958 &&
	    m_flavs[m_pi[1]].Kfcode()==kf_K_plus))  m_mode = FF_0_PP_mode::Ketaprime_plus;
  // pi-eta / pi-eta' (tau -> pi^- eta(') nu): genuine second-class
  // currents - per the cross-check note (Escribano, Gonzalez-Solis,
  // Roig, arXiv:1601.03989 = EGSR16), the normalized vector form
  // factor is IDENTICAL to the pi pi one (just an overall
  // isospin-breaking mixing factor epsilon_pieta(') out front), and
  // the scalar form factor is a dressed a0(980)[+a0(1450)] pole. See
  // Fplus/Fzero_0_PiZeroPiPlus's etapi_plus/etaprimepi_plus branches.
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta) ||
	   (m_flavs[m_pi[0]].Kfcode()==kf_eta &&
	    m_flavs[m_pi[1]].Kfcode()==kf_pi_plus))  m_mode = FF_0_PP_mode::etapi_plus;
  else if ((m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[m_pi[1]].Kfcode()==kf_eta_prime_958) ||
	   (m_flavs[m_pi[0]].Kfcode()==kf_eta_prime_958 &&
	    m_flavs[m_pi[1]].Kfcode()==kf_pi_plus))  m_mode = FF_0_PP_mode::etaprimepi_plus;
}

Complex FF_0_PP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
  // Fallback policy: ff_model::none is an explicit, deliberate "this
  // form factor is identically zero" (e.g. F0 for pipi_plus/KK_plus,
  // vanishing by SU(3)/G-parity) and stays hard zero. Every other case
  // where a specific channel+model combination has no real dynamical
  // implementation falls back to a CONSTANT form factor equal to
  // m_norm (i.e. F=1 times whatever CKM/isospin/normalisation factors
  // are already folded into m_norm for this channel) rather than
  // silently vanishing - see ff_model::point below for the same
  // convention made explicit.
  switch (m_ffmodel) {
  case ff_model::none:    return Complex(0.,0.);
  case ff_model::point:   return Complex(m_norm,0.);
  // KS(100)/GS(101)/KS_CLEO(102): try the specifically-requested
  // model's own propagator first; if that channel/model combination
  // has nothing (p_props==NULL), fall back to the family-base
  // (KS=100) construction if IT exists for this channel
  // (p_props_fallback, request #2); only fall back further to the
  // flat constant if even that isn't available.
  case ff_model::KS:
  case ff_model::GS:
  case ff_model::KS_CLEO:
    if (p_props!=NULL)          return m_norm * (*p_props)(Q2);
    if (p_props_fallback!=NULL) return m_norm * (*p_props_fallback)(Q2);
    return Complex(m_norm,0.);
  case ff_model::RChiPT:  return m_norm * FF_RChiPT(Q2);
  case ff_model::combRChL:
    msg_Error()<<"Error in "<<METHOD<<": combRChL form factor requested for "
	       <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<", but it is not "
	       <<"implemented yet (missing dispersion integral / B22 loop "
	       <<"functions - see FormFactor_Base.H). Falling back to a "
	       <<"constant (isospin/CKM-only) form factor.\n";
    return Complex(m_norm,0.);
  case ff_model::RChL2012: return m_norm * FF_RChL2012(Q2);
  case ff_model::unknown:
  default:
    break;
  }
  return Complex(m_norm,0.);
}

Complex FF_0_PP_Base::FF_RChiPT(const double & Q2) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<" form factor.\n"
	     <<"   Falling back to a constant (isospin/CKM-only) form "
	     <<"factor (F=1 times whatever is already in m_norm).\n";
  return Complex(1.,0.);
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// pi^0 pi^+ form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS, 1: (Kuehn-Santamaria model):
//   * pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
//   * pi pi (TODO: KS with Gounaris-Sakurai form factors):
//     from Belle measurement & analysis:
//     (https://arxiv.org/pdf/0805.3773) 
//   * K pi: Z.Phys.C 69 (1996) 243 (vector form factor)
//     (https://doi.org/10.1007/s002880050024)
// - RChT, 2: Resonance Chiral Perturbation Theory
//   * pi pi/KK fit: Eur.Phys.J.C 79 (2019) 5, 436
//     (https://doi.org/10.1140/epjc/s10052-019-6943-9)
//           - We could also have some alternative dispersive approach, 
//             maybe as form factor model 3, based on RchT
//             (https://arxiv.org/pdf/1112.0962)
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class Fplus_0_PiZeroPiPlus : public FF_0_PP_Base {
  double  m_fpi, m_mpi2, m_mK2;
  double  m_mV, m_mV2, m_mVp, m_mVp2, m_mVpp, m_mVpp2, m_GVp, m_GVpp;
  Complex m_gamma, m_delta;
  // Family-base (ff_model::KS=100-equivalent) mixing weights, ALWAYS
  // computed regardless of which model was actually requested. Used
  // both as the direct KS(100) result and, per request #2, as the
  // second-tier fallback for GS(101)/KS_CLEO(102) whenever THEIR own
  // specific construction has nothing for this channel (e.g. GS has
  // no published Kpi variant, so Kpi_plus+GS falls back to the plain
  // FM95 K*(892)+K*(1410) tune rather than straight to a flat constant).
  Complex m_gamma_base, m_delta_base;
  // Kpi_plus/Keta_plus only: true iff we are actually building the
  // ff_model::KS_CLEO-specific K*(892)+K*(1680) propagator (as opposed
  // to the family-base FM95 K*(892)+K*(1410) one). See the cross-check
  // note tau_two_meson_currents_KS_RChiT.tex, which flags these as TWO
  // genuinely different, non-interchangeable historical tunes that
  // should never be silently conflated under one label.
  bool    m_useCLEO;

  // --- RChiPT Kpi_plus "Level 1: bare RChiPT" (Jamin, Pich, Portoles;
  //     tau_two_meson_currents_KS_RChiT.tex Eq.(eq:RchiT-Kpi-vector)) ---
  // Separate members from the RChL2012 ones below (m_MKst etc.), even
  // though physically the "same" K*(892)/K*(1410), because the two
  // models fit genuinely different mass/gamma values (see
  // Construct_Kpi()).
  double  m_MKst_RChiPT, m_MKst2_RChiPT, m_MKstp_RChiPT, m_MKstp2_RChiPT;
  Complex m_gammaKpi_RChiPT;

  // --- RChL2012 (arXiv:1203.3955, Sec.2.4) specific members ---
  // F is the RChL pion decay constant (distinct normalisation from
  // m_fpi above, which is used by the RChiPT branch only).
  double  m_F_rchl, m_mpi0_2, m_mK0_2, m_meta2, m_mu2_rchl;
  double  m_MKst, m_MKst2, m_MKstp, m_MKstp2;
  Complex m_gammaKpi;
  Complex m_gammaRC, m_deltaRC; // gamma*e^{i phi1}, delta*e^{i phi2}, Eq.(24)
  // NOTE: rather than re-deriving the analytic energy-dependent widths
  // of Eqs.(28)-(34) of 1203.3955, we reuse Sherpa's own dynamic
  // Total_Width_Base widths for rho(770/1450/1700) and K*(892/1410),
  // already available via Line_Shapes. Both are physically motivated
  // running-width prescriptions; they are NOT guaranteed to be
  // numerically identical to the paper's formulas, so this branch
  // should be cross-checked before being used for precision fits. The
  // RChiPT Kpi_plus branch above reuses the SAME p_GKst/p_GKstp
  // pointers (built once in Construct_Kpi(), shared between both
  // models).
  Total_Width_Base * p_GRho, * p_GRhop, * p_GRhopp, * p_GKst, * p_GKstp;

  // Master dispatcher (called once from the constructor) + one
  // self-contained method per channel, each of which sets m_norm/
  // model-specific parameters AND builds whatever propagators that
  // channel+model combination needs, all in one place - replacing the
  // old split FixParameters()/Construct() pair, which had grown long
  // enough (and duplicated the same per-channel/per-model dispatch
  // structure twice, once for parameters and once for propagators)
  // that keeping the two in sync was itself a source of bugs.
  void    Construct(const FF_Parameters & params);
  void    Construct_PiPi(const FF_Parameters & params);
  void    Construct_KK(const FF_Parameters & params);
  void    Construct_EtaPi(const FF_Parameters & params);
  void    Construct_Kpi(const FF_Parameters & params);
  void    Construct_Keta(const FF_Parameters & params);
  // Shared propagator-builder helpers, used by several of the
  // per-channel methods above.
  Summed_Propagator * BuildRhoTower(bool useGS,
				     const Complex & gamma,
				     const Complex & delta);
  Summed_Propagator * BuildKstarTower(bool useCLEO, const Complex & gamma);

  Complex FF_RChiPT(const double & s);
  Complex FF_RChL2012(const double & s);
  Complex A(const double & m2,const double & s,const double & mu2);
  double  Gamma_V(const double & s);
  double  Gamma_Vp(const double & s);
  double  Gamma_Vpp(const double & s);
public :
  Fplus_0_PiZeroPiPlus(const FF_Parameters & params);
};

Fplus_0_PiZeroPiPlus::Fplus_0_PiZeroPiPlus(const FF_Parameters & params)  :
  FF_0_PP_Base(params),
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
  m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
  m_mK2(sqr(Flavour(kf_K_plus).HadMass())),
  m_useCLEO(false),
  m_MKst_RChiPT(0.), m_MKst2_RChiPT(0.),
  m_MKstp_RChiPT(0.), m_MKstp2_RChiPT(0.), m_gammaKpi_RChiPT(0.,0.),
  p_GRho(NULL), p_GRhop(NULL), p_GRhopp(NULL), p_GKst(NULL), p_GKstp(NULL)
{
  Construct(params);
}

///////////////////////////////////////////////////////////////////////////
//
// Shared propagator-builder helpers.
//
///////////////////////////////////////////////////////////////////////////

Summed_Propagator * Fplus_0_PiZeroPiPlus::
BuildRhoTower(bool useGS, const Complex & gamma, const Complex & delta) {
  resonance_type type = (useGS ? resonance_type::GS : resonance_type::running);
  Propagator_Base * rho770  =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),  type);
  Propagator_Base * rho1450 =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)), type);
  Propagator_Base * rho1700 =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_1700_plus)), type);
  Summed_Propagator * props = new Summed_Propagator();
  props->Add(rho770,  Complex(1., 0.));
  props->Add(rho1450, gamma);
  props->Add(rho1700, delta);
  return props;
}

Summed_Propagator * Fplus_0_PiZeroPiPlus::
BuildKstarTower(bool useCLEO, const Complex & gamma) {
  Propagator_Base * Kstar892 =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Summed_Propagator * props = new Summed_Propagator();
  props->Add(Kstar892, Complex(1., 0.));
  if (!useCLEO) {
    Propagator_Base * Kstar1410 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    props->Add(Kstar1410, gamma);
  }
  else {
    // CLEO/TAUOLA-type tune: K*(892)+w*K*(1680) (kf 30313/30323,
    // confirmed against the real ATOOLS/Phys/Flavour_Tags.H).
    Total_Width_Base * wKst1680 = LineShapes->Get(Flavour(kf_K_star_1680_plus));
    if (wKst1680!=NULL) {
      props->Add(new BreitWigner(wKst1680), gamma);
    }
    else {
      msg_Error()<<"Error in "<<METHOD<<": missing K*(1680) lineshape for "
		 <<"the CLEO/TAUOLA-type vector form factor - falling back "
		 <<"to plain K*(892) (dropping the w term).\n";
    }
  }
  return props;
}

///////////////////////////////////////////////////////////////////////////
//
// Master dispatcher + one self-contained method per channel.
//
///////////////////////////////////////////////////////////////////////////

void Fplus_0_PiZeroPiPlus::Construct(const FF_Parameters & params) {
  switch (m_mode) {
  case FF_0_PP_mode::pipi_plus:       Construct_PiPi(params);  break;
  case FF_0_PP_mode::KK_plus:         Construct_KK(params);    break;
  case FF_0_PP_mode::etapi_plus:
  case FF_0_PP_mode::etaprimepi_plus: Construct_EtaPi(params); break;
  case FF_0_PP_mode::Kpi_plus:        Construct_Kpi(params);   break;
  case FF_0_PP_mode::Keta_plus:
  case FF_0_PP_mode::Ketaprime_plus:  Construct_Keta(params);  break;
  default: break;
  }
  // Diagnostic dump (request #1), common to every channel.
  std::string label = std::string("Fplus_0_PiZeroPiPlus, mode=")+
                       std::to_string(int(m_mode));
  DumpPropagatorStructure(label, int(m_ffmodel),
			   p_props!=NULL ? p_props : p_props_fallback);
}

void Fplus_0_PiZeroPiPlus::Construct_PiPi(const FF_Parameters & params) {
  // The isospin factor for <pi- pi0| dbar gamma^mu u |0> is sqrt(2):
  //   H^mu = sqrt(2) Vud (p_pi- - p_pi0)^mu Ftilde_+^pipi(s),
  // with Ftilde_+(0)=1. The note fixes this internally rather than by
  // quoting C_pipi directly: Eq.(eq:pieta-current-pm) gives c_piP^V =
  // sqrt(2) for the pi eta / pi eta' currents, and
  // Eq.(eq:pieta-vff-data) states F_+^{pi eta} = eps_pieta *
  // Ftilde_+^{pipi}, i.e. those channels carry the SAME normalised pipi
  // form factor and differ only by the isospin-breaking eps. So pi pi
  // must carry the identical sqrt(2), with eps -> 1 - see
  // Construct_EtaPi below, which already has it.
  // This sqrt(2) was previously missing here, making Gamma(tau -> pi-
  // pi0 nu) exactly a factor 2 too small (measured 0.499 x PDG before
  // this fix, for the highest-branching-fraction hadronic mode there is).
  m_norm = (*params.p_model)("Vud", Tools::Vud) * sqrt(2.);
  // "CLEO/KS-type" 3-resonance tune (rho(770)+rho(1450)+rho(1700),
  // weights -0.167/+0.050) - confirmed against
  // tau_two_meson_currents_KS_RChiT.tex Table "KSparams", row
  // "pipi (CLEO/KS-type)": M_rho'=1.408 GeV/Gamma=0.502 GeV,
  // M_rho''=1.700 GeV/Gamma=0.235 GeV. This is a DIFFERENT tune from
  // Kuhn-Santamaria's/Finkemeier-Mirkes' own original pipi fit
  // (beta_rho=-0.145, single rho(1370,0.510) term, no third resonance
  // - see the same note) - the two should not be conflated under one
  // "KS" label. This is the family-base (KS=100) weight set, ALWAYS
  // computed (request #2's fallback target for GS=101).
  m_gamma_base = ReadComplexParam(params.p_model,
				   "gammaMag_pipi2_100",-0.167,"gammaPhase_pipi2_100");
  m_delta_base = ReadComplexParam(params.p_model,
				   "deltaMag_pipi2_100", 0.050,"deltaPhase_pipi2_100");
  p_props_fallback = BuildRhoTower(false, m_gamma_base, m_delta_base);

  if (m_ffmodel==ff_model::KS) {
    p_props = p_props_fallback;
    p_props_fallback = NULL;
  }
  else if (m_ffmodel==ff_model::GS) {
    // TODO: defaults below are placeholder values, identical to the
    // KS fit. Replace with the actual Gounaris-Sakurai beta/gamma from
    // the Belle fit (Fujikawa et al., Phys.Rev.D78 (2008) 072006, the
    // fit arXiv:1509.09140 reports as best (chi^2=95.65)) once we
    // have the numbers on hand - do NOT trust these for physics
    // results yet. Given their own "_GS"-suffixed parameter names
    // (distinct from "_KS"), so the real fit can be entered later
    // without disturbing the KS values it currently mirrors.
    m_gamma = ReadComplexParam(params.p_model,
			       "gammaMag_pipi_GS",-0.167,"gammaPhase_pipi_GS");
    m_delta = ReadComplexParam(params.p_model,
			       "deltaMag_pipi_GS", 0.050,"deltaPhase_pipi_GS");
    p_props = BuildRhoTower(true, m_gamma, m_delta);
  }
  else if (m_ffmodel==ff_model::RChiPT) {
    // Resonance masses/widths now read from the particle database
    // (Flavour(kf_xxx).HadMass()/.Width()) rather than hardcoded fit-
    // specific literals, at the cost of some precision (these will
    // generally differ somewhat from the actual Guerrero-Pich-style
    // fitted values) - overridable via the parameters below if you
    // have the real fit numbers and want them back.
    m_mV    = Flavour(kf_rho_770_plus).HadMass();   m_mV2   = sqr(m_mV);
    m_mVp   = (*params.p_model)("MVp_pipi",
				 Flavour(kf_rho_1450_plus).HadMass());
    m_mVp2  = sqr(m_mVp);
    m_GVp   = (*params.p_model)("GVp_pipi",
				 Flavour(kf_rho_1450_plus).Width());
    m_mVpp  = (*params.p_model)("MVpp_pipi",
				 Flavour(kf_rho_1700_plus).HadMass());
    m_mVpp2 = sqr(m_mVpp);
    m_GVpp  = (*params.p_model)("GVpp_pipi",
				 Flavour(kf_rho_1700_plus).Width());
    m_gamma = ReadComplexParam(params.p_model,
			       "gammaMag_pipi2_200",-0.15,"gammaPhase_pipi2_200",-0.36);
    m_delta = ReadComplexParam(params.p_model,
			       "deltaMag_pipi2_200", 0.12,"deltaPhase_pipi2_200",-0.02);
  }
  else if (m_ffmodel==ff_model::RChL2012) {
    // Defaults from Table 4 of 1203.3955, channel 1 (pi- pi0); masses
    // now default to the particle database rather than hardcoded
    // literals.
    m_F_rchl   = (*params.p_model)("F_rchl",0.0924);
    m_mu2_rchl = sqr(Flavour(kf_rho_770_plus).HadMass()); // mu=Mrho, fn.42
    m_mpi0_2   = sqr(Flavour(kf_pi).HadMass());
    m_mK0_2    = sqr(Flavour(kf_K).HadMass());
    m_mVp      = (*params.p_model)("Mrhop_rchl",
				    Flavour(kf_rho_1450_plus).HadMass());
    m_mVp2     = sqr(m_mVp);
    m_mVpp     = (*params.p_model)("Mrhopp_rchl",
				    Flavour(kf_rho_1700_plus).HadMass());
    m_mVpp2    = sqr(m_mVpp);
    m_gammaRC = ReadComplexParam(params.p_model,
				 "gammaMag_pipi_RChL2012",0.14199,
				 "gammaPhase_pipi_RChL2012",-0.17377);
    m_deltaRC = ReadComplexParam(params.p_model,
				 "deltaMag_pipi_RChL2012",-0.12623,
				 "deltaPhase_pipi_RChL2012",0.27632);
    p_GRho   = LineShapes->Get(Flavour(kf_rho_770_plus));
    p_GRhop  = LineShapes->Get(Flavour(kf_rho_1450_plus));
    p_GRhopp = LineShapes->Get(Flavour(kf_rho_1700_plus));
  }
}

void Fplus_0_PiZeroPiPlus::Construct_KK(const FF_Parameters & params) {
  m_norm = (*params.p_model)("Vud", Tools::Vud)/sqrt(2.);
  // K_S/K_L projection (see the m_isKSKL comment in FF_0_PP.H): an
  // ADDITIONAL 1/sqrt(2) on top of the isospin factor above, since the
  // current only produces K0bar (never K^0), and <K_S,L|K0bar> =
  // -+1/sqrt(2). Applies equally to K_S and K_L (same |amplitude|^2).
  if (m_isKSKL) m_norm /= sqrt(2.);
  // Same "CLEO/KS-type" isovector rho tower as pipi_plus above -
  // tau_two_meson_currents_KS_RChiT.tex Table "KSparams", row
  // "KK (CLEO/KS-type)" quotes IDENTICAL numbers to the pipi row.
  m_gamma_base = ReadComplexParam(params.p_model,
				   "gammaMag_KK_100",-0.167,"gammaPhase_KK_100");
  m_delta_base = ReadComplexParam(params.p_model,
				   "deltaMag_KK_100", 0.050,"deltaPhase_KK_100");
  p_props_fallback = BuildRhoTower(false, m_gamma_base, m_delta_base);

  if (m_ffmodel==ff_model::KS) {
    p_props = p_props_fallback;
    p_props_fallback = NULL;
  }
  else if (m_ffmodel==ff_model::GS) {
    // TODO: placeholder, see note above for pipi_plus/GS. Own "_GS"
    // parameter names so the eventual real fit doesn't collide with KS.
    m_gamma = ReadComplexParam(params.p_model,
			       "gammaMag_KK_GS",-0.167,"gammaPhase_KK_GS");
    m_delta = ReadComplexParam(params.p_model,
			       "deltaMag_KK_GS", 0.050,"deltaPhase_KK_GS");
    p_props = BuildRhoTower(true, m_gamma, m_delta);
  }
  else if (m_ffmodel==ff_model::RChiPT) {
    m_mV    = Flavour(kf_rho_770_plus).HadMass();   m_mV2   = sqr(m_mV);
    m_mVp   = (*params.p_model)("MVp_KK",
				 Flavour(kf_rho_1450_plus).HadMass());
    m_mVp2  = sqr(m_mVp);
    m_GVp   = (*params.p_model)("GVp_KK",
				 Flavour(kf_rho_1450_plus).Width());
    m_mVpp  = (*params.p_model)("MVpp_KK",
				 Flavour(kf_rho_1700_plus).HadMass());
    m_mVpp2 = sqr(m_mVpp);
    m_GVpp  = (*params.p_model)("GVpp_KK",
				 Flavour(kf_rho_1700_plus).Width());
    // Best-guess fix for a magnitude/sign-inconsistent legacy value -
    // still NOT verified against the original fit source; please
    // confirm or override.
    m_gamma = ReadComplexParam(params.p_model,
			       "gammaMag_KK_200",-0.15,"gammaPhase_KK_200",-1.88);
    // delta (rho'' mixing) defaults to 0 here - the KK channel's
    // RChiPT fit apparently doesn't need a third resonance - but is
    // now overridable like everywhere else, rather than hardcoded.
    m_delta = ReadComplexParam(params.p_model,
			       "deltaMag_KK_200",0.,"deltaPhase_KK_200",0.);
  }
  else if (m_ffmodel==ff_model::RChL2012) {
    // Eq.(26): the "default" (FFKKVEC=0) two-kaon form factor has NO
    // excited-resonance terms, only the loop-function exponential on
    // top of a plain rho(770) pole - m_gammaRC/m_deltaRC default to 0
    // here, matching that, but are now overridable like every other
    // channel/model combination in case a nonzero KK RChL2012 tune is
    // ever wanted.
    m_F_rchl   = (*params.p_model)("F_rchl",0.0924);
    m_mu2_rchl = sqr(Flavour(kf_rho_770_plus).HadMass());
    m_mpi0_2   = sqr(Flavour(kf_pi).HadMass());
    m_mK0_2    = sqr(Flavour(kf_K).HadMass());
    m_gammaRC  = ReadComplexParam(params.p_model,
				  "gammaMag_KK_RChL2012",0.,"gammaPhase_KK_RChL2012",0.);
    m_deltaRC  = ReadComplexParam(params.p_model,
				  "deltaMag_KK_RChL2012",0.,"deltaPhase_KK_RChL2012",0.);
    p_GRho   = LineShapes->Get(Flavour(kf_rho_770_plus));
    p_GRhop  = LineShapes->Get(Flavour(kf_rho_1450_plus));
    p_GRhopp = LineShapes->Get(Flavour(kf_rho_1700_plus));
  }
}

void Fplus_0_PiZeroPiPlus::Construct_EtaPi(const FF_Parameters & params) {
  // pi-eta/pi-eta' vector form factor, EGSR16 (arXiv:1601.03989)
  // Eqs.(eq:pieta-current-pm),(eq:pieta-vff-data): second-class
  // current with c_{piP}^V=sqrt(2) and an overall isospin-breaking
  // mixing normalisation epsilon_pieta(') multiplying the SAME
  // (normalized) rho tower used for pi-pi. Only ff_model::KS is
  // implemented; other models fall back to the family-base/generic
  // constant policy.
  double eps = (m_mode==FF_0_PP_mode::etapi_plus ?
		(*params.p_model)("eps_pieta", 0.017) :
		(*params.p_model)("eps_pietaprime", 0.004));
  m_norm = (*params.p_model)("Vud", Tools::Vud) * sqrt(2.) * eps;
  // Deliberately reuses the pipi_plus 2-body KS names (gammaMag_pipi2_100 etc.,
  // not new etapi-specific ones) - EGSR16 states this is literally the
  // SAME normalized rho tower as pipi, not just a similar one, so an
  // override entered for one applies to both by design.
  m_gamma_base = ReadComplexParam(params.p_model,
				   "gammaMag_pipi2_100",-0.167,"gammaPhase_pipi2_100");
  m_delta_base = ReadComplexParam(params.p_model,
				   "deltaMag_pipi2_100", 0.050,"deltaPhase_pipi2_100");
  p_props_fallback = BuildRhoTower(false, m_gamma_base, m_delta_base);
  if (m_ffmodel==ff_model::KS) {
    p_props = p_props_fallback;
    p_props_fallback = NULL;
  }
}

void Fplus_0_PiZeroPiPlus::Construct_Kpi(const FF_Parameters & params) {
  // This one method serves two physically DIFFERENT charge states, whose
  // isospin Clebsch factors are not the same:
  //   K^- pi^0   : Vus/sqrt(2)
  //   K0bar pi^- : Vus            (a factor sqrt(2) LARGER)
  // Isospin gives A(K0bar pi^-) = sqrt(2) A(K^- pi^0), which the measured
  // rates confirm: BR(K0bar pi^-) = 2 BR(K_S pi^-) = 0.0084 against
  // BR(K^- pi^0) = 0.0045, i.e. a ratio of 1.87 ~ 2.
  // The K0bar pi^- final state is then observed as K_S or K_L, each
  // carrying an ADDITIONAL 1/sqrt(2) amplitude projection. Those two
  // factors cancel exactly, so the K_S/K_L modes end up with the SAME
  // m_norm as K^- pi^0 - they must NOT get a further 1/sqrt(2) on top of
  // the K^- pi^0 Clebsch, which is what this code used to do, making
  // those rates a factor 2 too small (measured 0.553 x PDG for pi- K_L
  // and 0.415 for K_S pi-, against 1.073 for K- pi0 on the same model).
  m_norm = (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
  // Family-base (KS=100): FM95 (hep-ph/9503474), Eq.(FM95-Kpi):
  // beta_K*=-0.135+-0.025. ALWAYS computed - this is both the KS(100)
  // result itself and, per request #2, the fallback used whenever
  // GS(101) (no published Kpi variant) is requested for this channel.
  m_gamma_base = ReadComplexParam(params.p_model,
				   "gammaMag_Kpi_100",-0.135,"gammaPhase_Kpi_100");
  p_props_fallback = BuildKstarTower(false, m_gamma_base);

  if (m_ffmodel==ff_model::KS) {
    p_props = p_props_fallback;
    p_props_fallback = NULL;
  }
  else if (m_ffmodel==ff_model::KS_CLEO) {
    // CLEO/TAUOLA-type tune: K*(892)+w*K*(1680), w~=-0.038. The note
    // quotes M=1.700 GeV, Gamma=0.235 GeV for K*(1680) here, which may
    // or may not match whatever Sherpa's particle table has for that
    // kf-code (30313/30323) - see the identical caveat for the FM95
    // 3-body T_K*^(2) usage of the same resonance in FF_0_PPP.C.
    m_gamma   = ReadComplexParam(params.p_model,
				 "gammaMag_Kpi_102",-0.038,"gammaPhase_Kpi_102");
    m_useCLEO = true;
    p_props   = BuildKstarTower(true, m_gamma);
  }
  else if (m_ffmodel==ff_model::RChiPT) {
    // "Level 1: bare RChiPT" two-resonance vector form factor (Jamin,
    // Pich, Portoles - JPP), cross-check note
    // tau_two_meson_currents_KS_RChiT.tex Eq.(eq:RchiT-Kpi-vector):
    //   F~_+^{Kpi}(s) = [M_K*^2+gamma s]/D_K*(s) - gamma s/D_K*'(s),
    //   D_X(s) = M_X^2 - s - i M_X Gamma_X(s).
    // Automatically normalized to 1 at s=0 - no explicit 1/(1+gamma)
    // prefactor needed, unlike KS. Deliberately does NOT include the
    // chiral-loop unitarization exponential the pipi_plus branch
    // uses - that belongs to the note's "Level 2: unitarized"
    // treatment (needs the unequal-mass Kpi loop function
    // Htilde_Kpi(s)) - NOT implemented here; flag if wanted.
    // Masses now default to the particle database rather than the
    // note's "Combined disp. 2014" fit values (0.89203/1.304 GeV) -
    // the K*(1410) default in particular will come out noticeably
    // higher than that fit's shifted pole (1.304 GeV), a real loss of
    // precision accepted here for simplicity; override via the two
    // parameters below if you want the fitted values back. gamma is
    // a real number in JPP's own Eq. (no phase), but is now Complex-
    // capable like every other mixing weight in this file - the
    // default phase of 0 exactly reproduces the original real value.
    m_MKst_RChiPT  = (*params.p_model)("MKst",
					Flavour(kf_K_star_892_plus).HadMass());
    m_MKst2_RChiPT = sqr(m_MKst_RChiPT);
    m_MKstp_RChiPT = (*params.p_model)("MKstp",
					Flavour(kf_K_star_1410_plus).HadMass());
    m_MKstp2_RChiPT= sqr(m_MKstp_RChiPT);
    m_gammaKpi_RChiPT = ReadComplexParam(params.p_model,
					 "gammaMag_Kpi_200",-0.034,"gammaPhase_Kpi_200");
    p_GKst  = LineShapes->Get(Flavour(kf_K_star_892_plus));
    p_GKstp = LineShapes->Get(Flavour(kf_K_star_1410_plus));
  }
  else if (m_ffmodel==ff_model::RChL2012) {
    // Eq.(27), default (FFKPIVEC=1) parameters from Table 4; masses
    // now default to the particle database rather than the paper's
    // own fitted literals (0.8953/1.307 GeV) - see the identical
    // precision caveat for the RChiPT branch above. gammaKpi is real
    // in the paper's own Eq.(27) but, as above, is now Complex-capable
    // with a default phase of 0.
    m_F_rchl    = (*params.p_model)("F_rchl",0.0924);
    m_MKst      = (*params.p_model)("MKst_rchl",
				     Flavour(kf_K_star_892_plus).HadMass());
    m_MKst2     = sqr(m_MKst);
    m_MKstp     = (*params.p_model)("MKstp_rchl",
				     Flavour(kf_K_star_1410_plus).HadMass());
    m_MKstp2    = sqr(m_MKstp);
    m_gammaKpi  = ReadComplexParam(params.p_model,
				   "gammaMag_Kpi_RChL2012",-0.043,"gammaPhase_Kpi_RChL2012");
    m_mu2_rchl  = m_MKst2; // mu=MK* for the Kpi channel, footnote 42
    m_meta2     = sqr(Flavour(kf_eta).HadMass());
    p_GKst  = LineShapes->Get(Flavour(kf_K_star_892_plus));
    p_GKstp = LineShapes->Get(Flavour(kf_K_star_1410_plus));
  }
}

void Fplus_0_PiZeroPiPlus::Construct_Keta(const FF_Parameters & params) {
  // "KS-like" Keta vector form factor - per the cross-check note,
  // there is NO genuine original Kuhn-Santamaria/FM Keta model; this
  // is a natural KS-STYLE extension (same K* resonance content as
  // Kpi_plus), not a literal historical result. N_Keta/N_Ketaprime
  // below are pure normalisation placeholders (default 1) standing in
  // for the eta-eta' mixing factor that should multiply each channel -
  // FIXME: please supply the mixing convention/angle you want used
  // (e.g. the standard octet-singlet theta_P=-13.3+-0.5deg or
  // quark-flavour phi_P=41.4+-0.5deg given in the note) so these can
  // be replaced by the correct cos/sin factors. The note itself flags
  // this explicitly for Keta' (Eq.(eq:Ketap-current)'s C_{Keta'}):
  // "should be stored explicitly ... rather than absorbed silently
  // into a Clebsch factor" - same reasoning applies to Keta.
  bool isPrime = (m_mode==FF_0_PP_mode::Ketaprime_plus);
  m_norm = (*params.p_model)(isPrime ? "N_Ketaprime" : "N_Keta", 1.) *
           (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
  // Family-base (KS=100): same placeholder -0.135 as Kpi_plus, own
  // parameter names per channel (Keta vs Ketaprime) so an eventual
  // distinct fit for either doesn't collide with the other.
  m_gamma_base = ReadComplexParam(params.p_model,
				   isPrime ? "gammaMag_Ketaprime_100" : "gammaMag_Keta_100",
				   -0.135,
				   isPrime ? "gammaPhase_Ketaprime_100" : "gammaPhase_Keta_100");
  p_props_fallback = BuildKstarTower(false, m_gamma_base);

  if (m_ffmodel==ff_model::KS) {
    p_props = p_props_fallback;
    p_props_fallback = NULL;
  }
  else if (m_ffmodel==ff_model::KS_CLEO) {
    // CLEO/TAUOLA-type Keta tune, Eq.(eq:CLEO-Keta) of the note:
    // K*(892) - 0.038*K*(1680) (the sign is now folded directly into
    // the default magnitude below, rather than negated in code). No
    // dedicated CLEO-tune source is given for Keta' in the note, but
    // the same K*(892)/K*(1680) structure is used here by the same
    // "same strange vector current" argument the note itself makes
    // for Keta - flag if you'd rather this fell back to constant
    // instead for Ketaprime_plus specifically.
    m_gamma   = ReadComplexParam(params.p_model,
				 isPrime ? "gammaMag_Ketaprime_102" : "gammaMag_Keta_102",
				 -0.038,
				 isPrime ? "gammaPhase_Ketaprime_102" : "gammaPhase_Keta_102");
    m_useCLEO = true;
    p_props   = BuildKstarTower(true, m_gamma);
  }
  else if (m_ffmodel==ff_model::RChiPT) {
    // "Level 1: bare RChiPT" (Jamin, Pich, Portoles), cross-check note
    // Eqs.(eq:RchiT-Keta-relation)/(eq:Ketap-vector-relation): the
    // SAME normalized K*(892)+K*(1410) two-resonance vector form
    // factor as Kpi (Eq.(eq:RchiT-Kpi-vector)) governs Keta and Keta'
    // too, "up to the channel-specific thresholds entering the
    // running widths" - which the note explicitly says should still
    // use the Kpi-dominated K*(1410) running width rather than an ad
    // hoc pure-Keta(')-width replacement (hence reusing p_GKst/
    // p_GKstp unchanged below, exactly as for Kpi_plus). The
    // "Combined 2014" fit quoted in the note's Keta parameter table is
    // in fact a JOINT Kpi-Keta fit, so the numeric defaults below are
    // identical to Kpi_plus's own RChiPT defaults (own parameter
    // names per channel, in case a dedicated Keta-only or Keta'-only
    // fit ever supersedes the joint one).
    m_MKst_RChiPT  = (*params.p_model)("MKst",
					Flavour(kf_K_star_892_plus).HadMass());
    m_MKst2_RChiPT = sqr(m_MKst_RChiPT);
    m_MKstp_RChiPT = (*params.p_model)("MKstp",
					Flavour(kf_K_star_1410_plus).HadMass());
    m_MKstp2_RChiPT= sqr(m_MKstp_RChiPT);
    m_gammaKpi_RChiPT = ReadComplexParam(params.p_model,
					 isPrime ? "gammaMag_Ketaprime_200" : "gammaMag_Keta_200",
					 -0.034,
					 isPrime ? "gammaPhase_Ketaprime_200" : "gammaPhase_Keta_200");
    p_GKst  = LineShapes->Get(Flavour(kf_K_star_892_plus));
    p_GKstp = LineShapes->Get(Flavour(kf_K_star_1410_plus));
  }
}


Complex Fplus_0_PiZeroPiPlus::FF_RChiPT(const double & s) {
  if (m_mode==FF_0_PP_mode::Kpi_plus ||
      m_mode==FF_0_PP_mode::Keta_plus ||
      m_mode==FF_0_PP_mode::Ketaprime_plus) {
    // Level-1 bare RChiPT (JPP), Eq.(eq:RchiT-Kpi-vector) for Kpi_plus;
    // Eqs.(eq:RchiT-Keta-relation)/(eq:Ketap-vector-relation) state the
    // SAME formula (with its own per-channel gamma, same M_K*/M_K*')
    // governs Keta_plus/Ketaprime_plus too - see Construct_Keta()'s
    // comment for the full derivation/caveats.
    if (p_GKst==NULL || p_GKstp==NULL) return Complex(1.,0.); // constant fallback
    double GKst_s  = (*p_GKst)(s);
    double GKstp_s = (*p_GKstp)(s);
    // D_X(s) = M_X^2 - s - i*M_X*Gamma_X(s), per Eq.(eq:Kpi-D-resummed)
    // (the simpler bare-Level-1 form uses the same sign convention).
    Complex D1(m_MKst2_RChiPT -s, -m_MKst_RChiPT *GKst_s);
    Complex D2(m_MKstp2_RChiPT-s, -m_MKstp_RChiPT*GKstp_s);
    Complex term1 = (Complex(m_MKst2_RChiPT,0.) + m_gammaKpi_RChiPT*s) / D1;
    Complex term2 = (m_gammaKpi_RChiPT*s) / D2;
    return term1 - term2;
  }
  // pipi_plus/KK_plus: original 3-resonance (rho+rho'+rho'') chiral-
  // loop-unitarized form, unchanged.
  Complex V   = ( (m_mV2 - s*(m_gamma+m_delta))                       *
		  Complex(m_mV2-s, m_mV*Gamma_V(s))/
		  (sqr(m_mV2-s)+sqr(m_mV*Gamma_V(s)))             *
		  exp(-s/(96*sqr(M_PI*m_fpi)) *
		      (A(m_mpi2,s,m_mV2)+1./2.*A(m_mK2,s,m_mV2)).real()) );
  Complex Vp  = ( s*m_gamma                                             *
		  Complex(m_mVp2-s, m_mVp*Gamma_Vp(s))/
		  (sqr(m_mVp2-s)+sqr(m_mVp*Gamma_Vp(s)))           *
		  exp(-s*Gamma_Vp(m_mVp2)/
		      (M_PI*pow(m_mVp,3.)*
		       pow(1.-4.*m_mpi2/m_mVp2,3./2.)) *
		      A(m_mpi2,s,m_mV2).real())                           );
  Complex Vpp = ( s*m_delta                                             *
		  Complex(m_mVpp2-s, m_mVpp*Gamma_Vpp(s))/
		  (sqr(m_mVpp2-s)+sqr(m_mVpp*Gamma_Vpp(s)))        *
		  exp(-s*Gamma_Vpp(m_mVpp2)/
		      (M_PI*pow(m_mVpp,3.)*
		       pow(1.-4.*m_mpi2/m_mVpp2,3./2.)) *
		      A(m_mpi2,s,m_mV2).real())                          );
  return V + Vp + Vpp;
}

Complex Fplus_0_PiZeroPiPlus::A(const double & m2,const double & s,const double & mu2) {
  Complex sigma = csqrt(1.-4.*m2/s);
  return log(m2/mu2)+8.*m2/s-5./3.+pow(sigma,3.)*log((sigma+1.)/(sigma-1.));
}

///////////////////////////////////////////////////////////////////////////
//
// RChL two-meson vector form factor, arXiv:1203.3955 Sec.2.4.
//  - pipi_plus: Eq.(24) (rho + rho' + rho'' with complex mixing).
//  - KK_plus:   Eq.(26) (rho only - m_gammaRC/m_deltaRC are 0 by
//               construction for this mode, see Construct_KK() above,
//               so the rho'/rho'' terms below vanish automatically).
//  - Kpi_plus:  Eq.(27) (K*(892) + K*'(1410)).
// The A_PQ(s)/A(s) loop-function exponentials use RChL::A_PQ (see that
// file's header for the one known approximation, Jbar'(0)->0) and the
// already-existing A(m2,s,mu2) helper (=Eq.(50)'s equal-mass Api(s)).
//
///////////////////////////////////////////////////////////////////////////

Complex Fplus_0_PiZeroPiPlus::FF_RChL2012(const double & s) {
  if (m_mode==FF_0_PP_mode::pipi_plus || m_mode==FF_0_PP_mode::KK_plus) {
    if (p_GRho==NULL) return Complex(0.,0.);
    double GammaRho_s = (*p_GRho)(s);
    // Eq.(23)'s exponent: -s/(96 pi^2 F^2) [ReA_pi-pi0(s) + 1/2 ReA_K-K0(s)]
    Complex Apipi = RChL::A_PQ(s,m_mpi0_2,m_mpi2,m_F_rchl,m_mu2_rchl);
    Complex AKK   = RChL::A_PQ(s,m_mK0_2, m_mK2, m_F_rchl,m_mu2_rchl);
    double loopExp = exp( -s/(96.*sqr(M_PI)*sqr(m_F_rchl)) *
			  (Apipi.real() + 0.5*AKK.real()) );
    Complex leadTerm = ( Complex(m_mV2,0.) + s*(m_gammaRC+m_deltaRC) ) /
                       Complex(m_mV2-s,-sqrt(m_mV2)*GammaRho_s) * loopExp;
    Complex rhopTerm(0.,0.), rhoppTerm(0.,0.);
    if (norm(m_gammaRC)>0. && p_GRhop!=NULL) {
      double GammaRhop_s = (*p_GRhop)(s);
      double sigpi3 = pow(RChL::SigmaPhaseSpace(m_mVp2,m_mpi2),3.);
      double GRhopPole = (*p_GRhop)(m_mVp2);
      double expFac = ( sigpi3>0. ?
			exp(-s*GRhopPole/(M_PI*pow(m_mVp,3.)*sigpi3) *
			    A(m_mpi2,s,m_mV2).real()) : 0. );
      rhopTerm = - s*m_gammaRC / Complex(m_mVp2-s,-m_mVp*GammaRhop_s) * expFac;
    }
    if (norm(m_deltaRC)>0. && p_GRhopp!=NULL) {
      double GammaRhopp_s = (*p_GRhopp)(s);
      double sigpi3 = pow(RChL::SigmaPhaseSpace(m_mVpp2,m_mpi2),3.);
      double GRhoppPole = (*p_GRhopp)(m_mVpp2);
      double expFac = ( sigpi3>0. ?
			exp(-s*GRhoppPole/(M_PI*pow(m_mVpp,3.)*sigpi3) *
			    A(m_mpi2,s,m_mV2).real()) : 0. );
      rhoppTerm = - s*m_deltaRC / Complex(m_mVpp2-s,-m_mVpp*GammaRhopp_s) * expFac;
    }
    return leadTerm + rhopTerm + rhoppTerm;
  }
  else if (m_mode==FF_0_PP_mode::Kpi_plus) {
    if (p_GKst==NULL) return Complex(0.,0.);
    double GammaKst_s  = (*p_GKst)(s);
    double GammaKstp_s = (p_GKstp!=NULL ? (*p_GKstp)(s) : 0.);
    Complex term1 = (m_MKst2+s*m_gammaKpi) /
                    Complex(m_MKst2-s,-m_MKst*GammaKst_s);
    Complex term2 = ( p_GKstp!=NULL ?
		     s*m_gammaKpi / Complex(m_MKstp2-s,-m_MKstp*GammaKstp_s) :
		     Complex(0.,0.) );
    Complex AKpi = RChL::A_PQ(s,m_mK2, m_mpi2,m_F_rchl,m_mu2_rchl);
    Complex AKeta= RChL::A_PQ(s,m_mK2, m_meta2,m_F_rchl,m_mu2_rchl);
    double loopExp = exp( -s/(128.*sqr(M_PI)*sqr(m_F_rchl)) *
			  (AKpi.real()+AKeta.real()) );
    return (term1 - term2) * loopExp;
  }
  return Complex(0.,0.);
}

double Fplus_0_PiZeroPiPlus::Gamma_V(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  double pref = m_mV*s/(96.*M_PI*sqr(m_fpi));
  double arg  = (pow(1.-4.*m_mpi2/s,3./2.) + 
		 ((s>4.*m_mK2) ? 1./2. * pow(1.-4.*m_mK2/s,3./2.) : 0.));
  return pref*arg;
}

double Fplus_0_PiZeroPiPlus::Gamma_Vp(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  // Eq.(eq:RchiT-rhop-width-complete): Gamma(s) = Gamma * (s/M^2) *
  // sigma_pi^3(s)/sigma_pi^3(M^2), with sigma_pi(x)=sqrt(1-4mpi^2/x) -
  // i.e. the VELOCITY ratio, same convention as Gamma_V above. An
  // earlier version used the momentum ratio ((s-4mpi^2)/(M^2-4mpi^2))
  // ^(3/2), which differs by a factor (s/M^2)^(3/2): exact at the pole
  // (so it survives any on-shell width check) but a factor ~3 too
  // small at s=1 GeV^2, i.e. wrong right where the pi pi spectrum is.
  return (m_GVp * s/m_mVp2 *
	  pow(1.-4.*m_mpi2/s,3./2.)/pow(1.-4.*m_mpi2/m_mVp2,3./2.));
}

double Fplus_0_PiZeroPiPlus::Gamma_Vpp(const double & s) {
  if (s<4.*m_mpi2) return 0.;
  // Velocity ratio, exactly as in Gamma_Vp above - see the comment there.
  return (m_GVpp * s/m_mVpp2 *
	  pow(1.-4.*m_mpi2/s,3./2.)/pow(1.-4.*m_mpi2/m_mVpp2,3./2.));
}

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi & K K vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Scalar Kpi form factor F_S, Finkemeier & Mirkes hep-ph/9601275 (FM96),
// Eqs.(13),(16),(17). Only implemented for FORM_FACTOR=KS and the
// Kpi_plus mode (F_S vanishes identically for pipi_plus/KK_plus in the
// SU(3) flavour-symmetry limit m_K=m_pi, Sec.II - not literally 0 for
// real masses, but FM96 only derives/needs it for Kpi).
//
// F_S^(K0bar pi-)(Q^2) =
//   (m_K^2-m_pi^2)/Q^2 * cV/(1+beta_K*) *
//     [ (m_K*^2-Q^2)/m_K*^2   * BW_K*(Q^2)
//       + beta_K* (m_K*'^2-Q^2)/m_K*'^2 * BW_K*'(Q^2) ]
//   + (m_K^2-m_pi^2)/m_K*0^2 * cS * BW_K*0(Q^2)                    (17)
//
// with cV=1 (matched to ChPT, Eq.21), cS=1.7 (Eq.26, the paper itself
// stresses this number "should not be taken too seriously" given the
// large extrapolation from Q^2=0 to Q^2=m_K*0^2), beta_K*=-0.135
// (shared with the vector form factor, Sec.II of FM95/hep-ph/9503474).
// F^(K-pi0)_S = (1/sqrt2) F^(K0bar pi-)_S (Eq.8) - same Vus/sqrt2
// norm already used for the vector form factor in Fplus above.
//
///////////////////////////////////////////////////////////////////////////

class Fzero_0_PiZeroPiPlus : public FF_0_PP_Base {
  double  m_fpi, m_mpi2, m_mK2;
  double  m_mV, m_mV2, m_mVp, m_mVp2, m_mVpp, m_mVpp2, m_GVp, m_GVpp;
  Complex m_gamma, m_delta;

  // --- FM96 scalar Kpi (hep-ph/9601275) ---
  double  m_betaKst, m_cV, m_cS;
  double  m_mKst2, m_mKstp2, m_mKst02;
  Propagator_Base * p_BWKst, * p_BWKstp, * p_BWKst0;

  // --- pi-eta/pi-eta' scalar form factor, EGSR16 arXiv:1601.03989
  //     Eq.(eq:pieta-sff-one-dressed): dressed a0(980) pole. ---
  double  m_c0, m_metaP2, m_Ma0_2; // m_metaP2 = m_P^2 (eta or eta', mode-dependent)
  Propagator_Base * p_BWa0;

  // Master dispatcher + one self-contained method per channel - same
  // rationale as Fplus_0_PiZeroPiPlus above (sets m_norm/parameters
  // AND builds propagators together, replacing the old split
  // FixParameters()/Construct() pair).
  void    Construct(const FF_Parameters & params);
  void    Construct_Kpi(const FF_Parameters & params);
  void    Construct_Keta(const FF_Parameters & params);
  void    Construct_EtaPi(const FF_Parameters & params);
  Complex FF_RChiPT(const double & Q2) { return Complex(0.,0.); }
public :
  Fzero_0_PiZeroPiPlus(const FF_Parameters & params)  :
    FF_0_PP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
    m_mpi2(sqr(Flavour(kf_pi_plus).HadMass())),
    m_mK2(sqr(Flavour(kf_K_plus).HadMass())),
    p_BWKst(NULL), p_BWKstp(NULL), p_BWKst0(NULL),
    p_BWa0(NULL)
 {
    Construct(params);
  }
  ~Fzero_0_PiZeroPiPlus() {
    if (p_BWKst)  { delete p_BWKst;  p_BWKst  = NULL; }
    if (p_BWKstp) { delete p_BWKstp; p_BWKstp = NULL; }
    if (p_BWKst0) { delete p_BWKst0; p_BWKst0 = NULL; }
    if (p_BWa0)   { delete p_BWa0;   p_BWa0   = NULL; }
  }
  Complex operator()(const ATOOLS::Vec4D_Vector& moms);
};

void Fzero_0_PiZeroPiPlus::Construct(const FF_Parameters & params) {
  switch (m_mode) {
  case FF_0_PP_mode::Kpi_plus:        Construct_Kpi(params);   break;
  case FF_0_PP_mode::Keta_plus:
  case FF_0_PP_mode::Ketaprime_plus:  Construct_Keta(params);  break;
  case FF_0_PP_mode::etapi_plus:
  case FF_0_PP_mode::etaprimepi_plus: Construct_EtaPi(params); break;
  default: break;
  }
  // Diagnostic dump (request #1). Fzero's poles are single independent
  // Propagator_Base objects (K*_0(1430) for Kpi_plus, a0(980) for
  // etapi_plus/etaprimepi_plus), not a Summed_Propagator, so dump
  // whichever one is relevant for this mode.
  std::string label = std::string("Fzero_0_PiZeroPiPlus, mode=")+
                       std::to_string(int(m_mode));
  Propagator_Base * dumpme = (p_BWKst0!=NULL ? p_BWKst0 :
			      (p_BWa0!=NULL ? p_BWa0 : NULL));
  DumpPropagatorStructure(label, int(m_ffmodel), dumpme);
}

void Fzero_0_PiZeroPiPlus::Construct_Kpi(const FF_Parameters & params) {
  // Set m_norm regardless of model, so the constant-fallback path in
  // operator() (any model other than KS) has a sensible non-zero
  // value to fall back to.
  m_norm = (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
  // K_S/K_L projection - see Fplus_0_PiZeroPiPlus::Construct_Kpi and
  // the m_isKSKL member comment in FF_0_PP.H for the reasoning.
  if (m_isKSKL) m_norm /= sqrt(2.);
  if (m_ffmodel!=ff_model::KS) return;

  m_betaKst  = (*params.p_model)("betaKst_scalar",-0.135);    // Sec.II of FM95
  m_cV       = (*params.p_model)("cV_scalar",1.);             // Eq.(21)
  // cS=1.7 is FM96's own rough estimate from matching the K*_0(1430)
  // pole residue to the ChPT scalar radius, Eq.(26) - the paper
  // itself flags this as not to be trusted quantitatively. Treat as
  // an adjustable input, not a precision constant.
  m_cS       = (*params.p_model)("cS_scalar",1.7);

  Total_Width_Base * wKst  = LineShapes->Get(Flavour(kf_K_star_892_plus));
  Total_Width_Base * wKstp = LineShapes->Get(Flavour(kf_K_star_1410_plus));
  if (wKst!=NULL)  { p_BWKst  = new BreitWigner(wKst);
		      m_mKst2  = sqr(wKst->Flav().Mass(true)); }
  if (wKstp!=NULL) { p_BWKstp = new BreitWigner(wKstp);
		      m_mKstp2 = sqr(wKstp->Flav().Mass(true)); }
  // K*_0(1430) lineshape confirmed and registered (Kstar_Decays.H/.C,
  // kf_K_0_star_1430(_plus)=10311/10321 per the real
  // ATOOLS/Phys/Flavour_Tags.H - note the underscore placement,
  // "K_0_star" not "K_star_0").
  Total_Width_Base * wKst0 = LineShapes->Get(Flavour(kf_K_0_star_1430_plus));
  if (wKst0!=NULL) { p_BWKst0 = new BreitWigner(wKst0);
		      m_mKst02 = sqr(wKst0->Flav().Mass(true)); }
  else msg_Error()<<"Error in "<<METHOD<<": missing K*_0(1430) lineshape "
		  <<"for the scalar Kpi form factor (FM96 hep-ph/9601275) "
		  <<"- the K*_0(1430) pole term will be dropped (returns "
		  <<"only the off-shell-vector-meson scalar projection, "
		  <<"which FM96 finds to be a small effect compared to "
		  <<"the K*_0 pole - see their Fig.1a/2a). Add a "
		  <<"K*_0(1430) Lineshape to Line_Shapes::Init() to fix.\n";
}

void Fzero_0_PiZeroPiPlus::Construct_Keta(const FF_Parameters & params) {
  // No scalar Keta/Keta' form factor exists in the literature (see the
  // cross-check note) - just set the constant fallback normalisation,
  // matching whichever N_Keta/N_Ketaprime override was used for the
  // vector form factor (Fplus_0_PiZeroPiPlus::Construct_Keta).
  bool isPrime = (m_mode==FF_0_PP_mode::Ketaprime_plus);
  m_norm = (*params.p_model)(isPrime ? "N_Ketaprime" : "N_Keta", 1.) *
           (*params.p_model)("Vus", Tools::Vus)/sqrt(2.);
}

void Fzero_0_PiZeroPiPlus::Construct_EtaPi(const FF_Parameters & params) {
  // Regularity (EGSR16, general PP current, Sec.2.1) requires
  // F_0(0)=F_+(0); F_+(0)=Vud*sqrt(2)*epsilon_pieta(') for the KS
  // vector form factor (normalized rho tower ->1 at s=0), so use the
  // SAME m_norm here for consistency.
  double eps = (m_mode==FF_0_PP_mode::etapi_plus ?
		(*params.p_model)("eps_pieta", 0.017) :
		(*params.p_model)("eps_pietaprime", 0.004));
  m_norm = (*params.p_model)("Vud", Tools::Vud) * sqrt(2.) * eps;
  m_metaP2 = (m_mode==FF_0_PP_mode::etapi_plus ?
	      sqr(Flavour(kf_eta).HadMass()) :
	      sqr(Flavour(kf_eta_prime_958).HadMass()));
  if (m_ffmodel!=ff_model::KS) return;

  // c_0^{piP}, EGSR16 Eq.(eq:pieta-sff-one-dressed): the paper fits
  // this directly and we don't have its numeric value from the
  // cross-check note - FIXME/please supply if you have the EGSR16
  // fit table. Default here instead auto-enforces the regularity
  // condition F_0(0)=F_+(0)=m_norm at s=0, which is a reasonable
  // placeholder but is NOT necessarily the paper's own fitted c_0.
  m_Ma0_2 = sqr((*params.p_model)("Ma0_980", 0.980));
  double bracket0 = 1. + (m_mpi2 - m_metaP2) / m_Ma0_2;
  m_c0 = (*params.p_model)("c0_pieta", m_norm/bracket0);

  // a0(980) lineshape confirmed (kf 9000111/9000211, decaying to
  // eta pi 100%/eta' pi with the same coupling once open) - see
  // A0_Decays.H/.C, now registered in Line_Shapes::Init(). The note
  // (EGSR16, arXiv:1601.03989) itself still flags a0(980) as better
  // treated via its coupled-channel (pi eta-KKbar threshold) line
  // shape than a simple single-resonance Breit-Wigner (Sec. on
  // Omnes/coupled channels) - the S_PP-based BW here remains a
  // simplification, not the preferred treatment, but is no longer
  // blocked on missing input. The optional two-resonance extension
  // with a0(1450) (Eq.(eq:pieta-sff-two)) is available as a
  // lineshape (A0_Decays.H) but not yet wired into this form factor -
  // flag if you want that added; it would need the c_m/c_m' coupling
  // split (only their sum-of-squares constraint is given in the
  // note, not fitted individual values).
  Total_Width_Base * wA0 = LineShapes->Get(Flavour(kf_a_0_980_plus));
  if (wA0!=NULL) { p_BWa0 = new BreitWigner(wA0); }
  else msg_Error()<<"Error in "<<METHOD<<": missing a0(980) lineshape "
		  <<"for the pi-eta(') scalar form factor (EGSR16 "
		  <<"arXiv:1601.03989) - F0 falls back to the constant "
		  <<"m_norm.\n";
}


Complex Fzero_0_PiZeroPiPlus::operator()(const ATOOLS::Vec4D_Vector& moms) {
  // pipi_plus/KK_plus: F0 is genuinely (SU(3)-limit) forced to vanish
  // here - see the class-level comment - this is physics, not a
  // missing implementation, so it stays hard zero regardless of model.
  if (m_mode==FF_0_PP_mode::pipi_plus || m_mode==FF_0_PP_mode::KK_plus) {
    return Complex(0.,0.);
  }
  // pi-eta/pi-eta': dressed a0(980) pole, EGSR16 Eq.(eq:pieta-sff-
  // one-dressed). Only ff_model::KS implemented; falls back to the
  // constant m_norm (already set to satisfy F_0(0)=F_+(0)) otherwise
  // or if the a0(980) lineshape is unavailable.
  if (m_mode==FF_0_PP_mode::etapi_plus || m_mode==FF_0_PP_mode::etaprimepi_plus) {
    if (m_ffmodel!=ff_model::KS || p_BWa0==NULL) return Complex(m_norm,0.);
    double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
    // (*p_BWa0)(Q2) = M_S^2/(M_S^2-Q2-i*sqrt(Q2)*Gamma(Q2)) in this
    // codebase's standard (FM/TAUOLA, Xi=sqrt(s)) convention; dividing
    // by m_Ma0_2 below cancels the M_S^2 numerator, reproducing the
    // note's 1/(M_S^2-s-i*Xi*Gamma(s)) denominator up to that same
    // Xi-convention choice (kept consistent with every other lineshape
    // in this file, rather than mixing in the M_R-convention some
    // R-chi-T papers use for this same equation).
    return m_c0 * (Complex(1.,0.) +
		   (Q2 + m_mpi2 - m_metaP2) * (*p_BWa0)(Q2) / m_Ma0_2);
  }
  // Kpi_plus/Keta_plus: the dynamical (FM96) scalar form factor is
  // only implemented for ff_model::KS. Any other model requested here
  // has no dedicated implementation, so falls back to the constant
  // (isospin/CKM-only) form factor m_norm, per the general fallback
  // policy - see FF_0_PP_Base::operator().
  if (!(m_ffmodel==ff_model::KS &&
	(m_mode==FF_0_PP_mode::Kpi_plus || m_mode==FF_0_PP_mode::Keta_plus))) {
    return Complex(m_norm,0.);
  }
  if (m_mode!=FF_0_PP_mode::Kpi_plus) {
    // Keta_plus/KS: no genuine scalar Keta form factor exists in the
    // literature either (see the cross-check note) - fall back to the
    // constant too, rather than reusing the Kpi-specific K*_0(1430)
    // dynamics for a channel it was never derived for.
    return Complex(m_norm,0.);
  }
  double Q2 = (moms[m_pi[0]]+moms[m_pi[1]]).Abs2();
  double massSplit = m_mK2 - m_mpi2;
  Complex vectorProjection(0.,0.);
  if (p_BWKst!=NULL) {
    vectorProjection += ((m_mKst2-Q2)/m_mKst2) * (*p_BWKst)(Q2);
  }
  if (p_BWKstp!=NULL) {
    vectorProjection += m_betaKst * ((m_mKstp2-Q2)/m_mKstp2) * (*p_BWKstp)(Q2);
  }
  vectorProjection *= m_cV/(1.+m_betaKst);
  Complex scalarPole(0.,0.);
  if (p_BWKst0!=NULL) {
    scalarPole = (massSplit/m_mKst02) * m_cS * (*p_BWKst0)(Q2);
  }
  Complex FS = (Q2!=0. ? (massSplit/Q2)*vectorProjection : Complex(0.,0.))
             + scalarPole;
  return m_norm * FS;
}


DECLARE_FF_GETTER(FF_0_PP_Base,"FF_0_PP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PP_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<":\n";
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    if (params.m_flavs[params.m_pi[i]].IsMeson()) Nmesons++;
  }
  if (Nmesons!=2) return NULL;
  // Below a first round of decays/currents for which we have both
  // Kuehn-Santamaria and RChiPT parametrizations
  if (//   pi^+ pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      //   K^+ K^0 (note - this could be K_S or K_L - just to play it safe
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
       //   pi^+ K^0 (note - this could be K_S or K_L - just to play it safe
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K)          ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_S)        ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_L)        ||
       //   K^+ pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi)         ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
       //   K^+ eta (tau -> K^- eta nu_tau) - "KS-like"/RChiPT, see
       //   Fplus_0_PiZeroPiPlus's Keta_plus branch.
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_eta)        ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
       //   K^+ eta' (tau -> K^- eta' nu_tau) - RChiPT (and KS-like/
       //   KS_CLEO by the same-current analogy), see
       //   Fplus_0_PiZeroPiPlus's Ketaprime_plus branch.
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_eta_prime_958) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta_prime_958 &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus)     ||
       //   pi^+ eta / pi^+ eta' (tau -> pi^- eta(') nu_tau), EGSR16
       //   second-class currents - see the etapi_plus/etaprimepi_plus
       //   branches in Fplus/Fzero_0_PiZeroPiPlus.
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_eta)        ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)    ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_eta_prime_958) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta_prime_958 &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus)
      ) {
    msg_Out()<<METHOD<<"("<<params.m_flavs[params.m_pi[0]]<<" + "
	     <<params.m_flavs[params.m_pi[1]]<<"): "
	     <<"model = "<<int(params.m_ffmodel)<<".\n";
    if (params.m_name=="F+_0_PP") return new Fplus_0_PiZeroPiPlus(params);
    if (params.m_name=="F0_0_PP") return new Fzero_0_PiZeroPiPlus(params);
  }
  return NULL;
}

