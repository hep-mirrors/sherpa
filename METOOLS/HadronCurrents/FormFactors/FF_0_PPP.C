#include "METOOLS/HadronCurrents/FormFactors/FF_0_PPP.H"
#include "METOOLS/HadronCurrents/FormFactors/RChL_Functions.H"
#include "METOOLS/HadronCurrents/FormFactors/K1_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Kstar_Decays.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <algorithm>

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


FF_0_PPP_Base::FF_0_PPP_Base(const FF_Parameters & params) :
  FormFactor_Base(params),
  m_mode(FF_0_PPP_mode::unknown),
  m_norm(Complex(0.,0.)),
  m_isKSKL(false) {
  FixMode();
}

void FF_0_PPP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi &&
      m_flavs[m_pi[2]].Kfcode()==kf_pi_plus)      m_mode = FF_0_PPP_mode::piP_pi0_pi0;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi_plus) m_mode = FF_0_PPP_mode::piM_piP_piP;
  // K^+ K^- pi^-: convention fixed in VA_0_KKPi.H (0=K+, 1=K-, 2=pi-).
  // Only this exact ordering is recognised for now - other index
  // permutations of the same physical final state are not yet mapped,
  // mirroring the existing minimal coverage of the pi pi pi modes above.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi_plus) m_mode = FF_0_PPP_mode::KP_KM_piM;
  // pi^0 pi^0 K^-: Bose-symmetric pair (the two pi^0's) + one odd
  // meson (the charged kaon) - same structure as piP_pi0_pi0 above,
  // just with the odd meson being a kaon instead of a pion. Handled
  // by the SAME VA_0_PiPiPi Current class (see F1_0_PiPlusPiZeroPiZero
  // below, which now also covers this mode).
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K_plus)     m_mode = FF_0_PPP_mode::pi0_pi0_KM;
  //
  // Below: Finkemeier-Mirkes hep-ph/9503474 (FM95) "KS"-mode channels.
  // None of these have two identical bosons, so p1,p2,p3 are simply
  // assigned in FM95's own q1,q2,q3 order per channel (see the
  // per-channel comments at FF_0_PPP.H's enum and reproduced at each
  // FF class below). These are routed to the new VA_0_KPiK/VA_0_KPiPi
  // Current classes, not VA_0_PiPiPi/VA_0_KKPi.
  //
  // K^-(q1) pi^-(q2) K^+(q3) - FM95 Tab.I/II row 1. NOTE this is the
  // SAME physical final state as KP_KM_piM above but with a DIFFERENT
  // momentum-index convention (q1=K-,q2=pi-,q3=K+ here, vs VA_0_KKPi.H's
  // p1=K+,p2=K-,p3=pi-) - the two are only compatible if routed through
  // their own dedicated Current class (VA_0_KPiK), never mixed.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K_plus)     m_mode = FF_0_PPP_mode::KM_piM_KP;
  // K^0(q1) pi^-(q2) K0bar(q3) - row 2.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K)          m_mode = FF_0_PPP_mode::K0_piM_K0bar;
  // K_S(q1) pi^-(q2) K_S(q3) / K_S pi^- K_L / K_L pi^- K_L - Sec.VI.
  // Recognised directly via the K_S/K_L kf-codes (distinct from K^0).
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_S &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K_S)        m_mode = FF_0_PPP_mode::KS_piM_KS;
  else if (((m_flavs[m_pi[0]].Kfcode()==kf_K_S &&
	     m_flavs[m_pi[2]].Kfcode()==kf_K_L) ||
	    (m_flavs[m_pi[0]].Kfcode()==kf_K_L &&
	     m_flavs[m_pi[2]].Kfcode()==kf_K_S)) &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus)    m_mode = FF_0_PPP_mode::KS_piM_KL;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_L &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K_L)        m_mode = FF_0_PPP_mode::KL_piM_KL;
  // K^-(q1) pi^0(q2) K^0(q3) - row 5. K_S/K_L recognised (with the
  // additional 1/sqrt(2) projection factor - same reasoning as
  // FF_0_PP.C's Kpi_plus/KK_plus and this file's piM_K0bar_pi0): the
  // K^-K_S/K^-K_L pi^0 channels flagged earlier as "not yet
  // recognised, needs new isospin coefficients" were WRONG - they are
  // just this EXISTING KM_pi0_K0 mode (FM95 Tab.I/II row 5) with its
  // K0bar observed as K_S/K_L, not a genuinely different current.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi &&
	   (m_flavs[m_pi[2]].Kfcode()==kf_K ||
	    m_flavs[m_pi[2]].Kfcode()==kf_K_S ||
	    m_flavs[m_pi[2]].Kfcode()==kf_K_L)) {
    m_mode   = FF_0_PPP_mode::KM_pi0_K0;
    m_isKSKL = (m_flavs[m_pi[2]].Kfcode()==kf_K_S ||
		m_flavs[m_pi[2]].Kfcode()==kf_K_L);
  }
  // K^-(q1) pi^-(q2) pi^+(q3) - row 7.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi_plus)    m_mode = FF_0_PPP_mode::KM_piM_piP;
  // pi^-(q1) K0bar(q2) pi^0(q3) - row 8. K_S/K_L recognised (with the
  // additional 1/sqrt(2) projection factor, see F1_0_KPiPi/FS_0_KPiPi
  // Construct()), same fix as FF_0_PP.C's Kpi_plus/KK_plus.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	   (m_flavs[m_pi[1]].Kfcode()==kf_K ||
	    m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
	    m_flavs[m_pi[1]].Kfcode()==kf_K_L) &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi) {
    m_mode   = FF_0_PPP_mode::piM_K0bar_pi0;
    m_isKSKL = (m_flavs[m_pi[1]].Kfcode()==kf_K_S ||
		m_flavs[m_pi[1]].Kfcode()==kf_K_L);
  }
  // eta(q1) pi^-(q2) pi^0(q3) / eta'(q1) pi^-(q2) pi^0(q3) - see the
  // enum comment in FF_0_PPP.H. Exact order required (q1=eta(prime)),
  // not interchangeable.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_eta &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi)         m_mode = FF_0_PPP_mode::EtaPiPi_pi0;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_eta_prime_958 &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi)         m_mode = FF_0_PPP_mode::EtaprimePiPi_pi0;
  msg_Out()<<METHOD<<"("<<m_pi[0]<<"/"<<m_pi[1]<<"/"<<m_pi[2]<<") --> "<<int(m_mode)<<"\n";
}

Complex FF_0_PPP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  Vec4D p1    = moms[m_pi[0]],  p2 = moms[m_pi[1]],   p3  = moms[m_pi[2]];
  Vec4D q     = p1+p2+p3;
  double s123 = q.Abs2(),      s12 = (p1+p2).Abs2(), s13 = (p1+p3).Abs2();
  // Fallback policy: any model requested for a channel/form-factor
  // slot with no real dynamical implementation falls back to a
  // CONSTANT form factor (F=1, i.e. m_norm alone - whatever CKM/
  // isospin factors are already folded into it), rather than silently
  // vanishing. ff_model::none already implements exactly this (see
  // below); the bottom fallback for unmatched/unknown models now does
  // the same, as do the FF_KS/FF_RChiPT base-class defaults below.
  switch (m_ffmodel) {
  case ff_model::none:    return m_norm;
  case ff_model::KS:
  case ff_model::KS_flatte: return m_norm * FF_KS(s123,s12,s13);
  // KS_CLEO (102): CLEO-fitted alternative parametrizations - 3pi's
  // "CLEO/default-TAUOLA" current and Kpipi's "CLEO K1 data-driven
  // alternative" (see tau_two_meson_currents_KS_RChiT.tex). Routed
  // through the SAME FF_KS() as plain KS - each derived class checks
  // m_ffmodel internally to pick the KS vs KS_CLEO formula, exactly
  // as FF_0_PP.C's Construct_XXX methods already do.
  case ff_model::KS_CLEO: return m_norm * FF_KS(s123,s12,s13);
  // KS_f0/CLEO_f0 (103/104): routed through the same FF_KS() entry
  // point - only F1_0_KPiPi/FS_0_KPiPi currently check m_ffmodel for
  // these two values internally (Kppipi + added f0(500) admixture).
  case ff_model::KS_f0:   return m_norm * FF_KS(s123,s12,s13);
  case ff_model::CLEO_f0: return m_norm * FF_KS(s123,s12,s13);
  case ff_model::RChiPT:  return m_norm * FF_RChiPT(s123,s12,s13);
  case ff_model::RChL2012: return m_norm * FF_RChL2012(s123,s12,s13);
  case ff_model::unknown:
  default:
    break;
  }
  return m_norm;
}

Complex FF_0_PPP_Base::
FF_KS(const double & s123,const double & s1,const double & s2) {
  msg_Error()<<"Error in "<<METHOD<<": KS form factor not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Falling back to a constant (isospin/CKM-only) form "
	     <<"factor (F=1 times whatever is already in m_norm).\n";
  return Complex(1.,0.);
}

Complex FF_0_PPP_Base::
FF_RChiPT(const double & s123,const double & s1,const double & s2) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Falling back to a constant (isospin/CKM-only) form "
	     <<"factor (F=1 times whatever is already in m_norm).\n";
  return Complex(1.,0.);
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// {pi pi pi}^+, {K pi pi}^+, {K K pi}^+  form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS
//     * pi pi pi from
// - none, 0: no form factor
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
protected:
  bool    m_isF2;
  double  m_fpi;
  Complex m_alpha, m_gamma, m_delta;
  
  Summed_Propagator * p_a1s, * p_rhos;

  // --- pi0_pi0_KM (FM95 hep-ph/9503474, Tab.I row 6 / Tab.II row 6) ---
  // Uses T_K1^(a) (K1(1400)+xi*K1(1270) mix) in place of the a1(1260),
  // and T_K*^(1) (K*(892)+beta_K* K*(1410) mix) in place of T_rho in
  // BOTH G1 and G2 (unlike the pion-only modes, this channel is not
  // symmetric under s1<->s2 in the propagator choice - see FF_KS).
  // FS_0_PiPlusPiZeroPiZero additionally needs T_K*^(2) for the vector
  // current (nonzero here, unlike the pure-pion modes where G-parity
  // forbids it) - see that class below.
  Summed_Propagator * p_TK1a, * p_TKstar1_pi0K;
  Complex m_xiK1;

  // --- RChL2012 (1203.3955 Sec.2.1 + 1310.1053 sigma extension) ---
  // Defaults below are the BaBar-fitted values of Table I of
  // 1310.1053 (pi- pi- pi+ channel), used for BOTH 3pi charge modes
  // per that paper's explicit recommendation (Sec. "Case of pi0pi0pi-
  // mode"): only alpha_sigma, gamma_sigma differ (scaled by 0.63 for
  // the neutral mode), everything else stays as in Table I. These are
  // genuine global-fit-specific numbers (not "the mass of a standard
  // resonance"), so - unlike the KS-family weights below - they are
  // NOT switched to Flavour(kf_xxx)-based defaults; only the KS-style
  // rho'/rho''/omega mixing weights get that (and the unified
  // gamma/delta/alpha nomenclature) treatment here.
  double  m_F_rchl, m_FV_rchl, m_FA_rchl, m_GV_rchl;
  double  m_Mrho_rchl, m_Mrhop_rchl, m_Grhop_rchl, m_betarhop, m_Ma1_rchl;
  double  m_Msigma, m_Gsigma0, m_Rsigma2;
  double  m_alphasigma, m_betasigma, m_gammasigma, m_deltasigma;
  double  m_mpi2_iso; // isospin-averaged pion mass^2 used throughout RChL
  double  m_lambda_p, m_lambda_pp;
  // p_rhos (declared above, shared with the KS branch) is reused here:
  // built via RChL_BW propagators + Summed_Propagator's weighted-sum-
  // then-normalise behaviour, which reproduces Eq.(11)'s combined
  // rho+rho' propagator exactly (Add(rho,1)+Add(rhop,betarhop), then
  // Summed_Propagator divides by (1+betarhop) automatically). p_a1s
  // is NOT reused for RChL2012: the a1 term needs a q^2 numerator (not
  // a plain BreitWigner's M^2) and uses the Eq.(17)-of-1310.1053
  // polynomial width fit rather than a Total_Width_Base/Line_Shapes
  // lineshape, neither of which fits the existing Propagator_Base
  // hierarchy without a larger refactor - so it stays inline in
  // FF_RChL2012 below.

  // --- CLEO / default-TAUOLA 3pi current (KS_CLEO=102) ---
  // tau_two_meson_currents_KS_RChiT.tex Sec."CLEO / default-TAUOLA 3pi
  // current". A genuinely different, data-driven parametrization from
  // KS90/RChL2012 above - own resonance masses/widths (ALL different
  // from the KS90 ones above, even for "the same" rho/a1), its own
  // running-width convention (BW_L, Eq.(eq:CLEO-BW)), and 7 complex
  // fitted amplitudes (beta1..beta7). Implemented as self-contained
  // functions using the literal CLEO masses/widths as plain numbers -
  // NOT tied to Sherpa's LineShapes/Total_Width_Base machinery, since
  // reproducing this specific parametrization exactly requires its
  // own closed-form running-width prescription (Eq.(eq:CLEO-running-
  // width)), not Sherpa's generic decay-channel-based running widths.
  //
  // BASIS REDUCTION (previously flagged as a caveat, now resolved by
  // the note): the note gives F1^CLEO, F2^CLEO, F3^CLEO as three
  // analytic functions, but the general PPPcurrent basis states "only
  // two transverse axial structures are independent" - F1,F2,F3
  // multiplying (p2-p3),(p3-p1),(p1-p2) with coefficients c1,c2,c3 are
  // NOT all independent (those three vectors sum to zero identically).
  // The note's "Active transverse basis in the default TAUOLA-CLEO
  // current" now states the active coefficients explicitly:
  //   c1=+2sqrt2/3,  c2=-2sqrt2/3,  c3=0,  c4=c5=0,
  // so the reduction is just F1_eff=F1^CLEO, F2_eff=F2^CLEO with NO F3
  // term (Eq.(eq:CLEO-two-slot-map)); the generic fold
  // F_i-(c3/c_i)F3 (Eq.(eq:CLEO-generic-fold)) vanishes. An earlier
  // version of this class assumed c1=c2=c3=1 and subtracted F3 - the
  // note's "Implementation warning" says explicitly that this is not
  // the default TAUOLA-CLEO convention. See FF_KS below for the slot
  // assignment that goes with it. CLEO_F3term is kept (correct as
  // written, coefficient zero in this current, useful for comparing
  // against a redundant three-vector basis).
  double  m_MCLEOrho, m_GCLEOrho, m_MCLEOrhop, m_GCLEOrhop;
  double  m_MCLEOsigma, m_GCLEOsigma, m_MCLEOf0, m_GCLEOf0, m_MCLEOf2, m_GCLEOf2;
  double  m_MCLEOa1, m_GCLEOa1, m_MCLEOa1p, m_GCLEOa1p;
  Complex m_betaCLEOa1p; // dormant a1' term, 0 in the nominal current
  double  m_MCLEOKstar, m_mCLEOK, m_mCLEOpi;
  Complex m_beta1,m_beta2,m_beta3,m_beta4,m_beta5,m_beta6,m_beta7;

  static double  CLEO_P(const double & S,const double & m1,const double & m2);
  static double  CLEO_GammaL(const double & S,const double & M,const double & Gamma,
			     const double & m1,const double & m2,int L);
  static Complex CLEO_BWL(const double & S,const double & M,const double & Gamma,
			  const double & m1,const double & m2,int L);
  double  CLEO_w1(const double & x) const;
  double  CLEO_w2(const double & x) const;
  double  CLEO_wK(const double & x) const;
  double  CLEO_WGA(const double & Q2) const;
  Complex CLEO_Aa1(const double & Q2) const;
  // F1^CLEO(Q2;s1,s2,s3,m1_2,m2_2,m3_2) - F2^CLEO is obtained by
  // calling this with (s2,s1,s3,m2_2,m1_2,m3_2), per the note's own
  // "F2=F1(1<->2)" statement.
  Complex CLEO_F1term(const double & Q2,const double & s1,const double & s2,
		      const double & s3,const double & m1_2,const double & m2_2,
		      const double & m3_2) const;
  Complex CLEO_F3term(const double & Q2,const double & s1,const double & s2,
		      const double & s3,const double & m1_2,const double & m2_2,
		      const double & m3_2) const;
  void    Construct_3Pi_CLEO(const FF_Parameters & params);

  // Master dispatcher (called once from the constructor) + one
  // self-contained method per channel - same rationale as
  // Fplus_0_PiZeroPiPlus in FF_0_PP.C: sets m_norm/parameters AND
  // builds propagators together, replacing the old split
  // FixParameters()/Construct() pair.
  void    Construct(const FF_Parameters & params);
  void    Construct_3Pi(const FF_Parameters & params);   // piP_pi0_pi0, piM_piP_piP
  void    Construct_Pi0Pi0K(const FF_Parameters & params); // pi0_pi0_KM

  Complex FF_KS(const double & s123,const double & s1,const double & s2);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2);
  Complex FF_RChL2012(const double & s123,const double & s1,const double & s2);
public :
  F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params);
  ~F1_0_PiPlusPiZeroPiZero();
};

F1_0_PiPlusPiZeroPiZero::F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
  FF_0_PPP_Base(params),
  p_a1s(NULL), p_rhos(NULL),
  p_TK1a(NULL), p_TKstar1_pi0K(NULL), m_xiK1(0.33,0.),
  m_isF2(false), 
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  Construct(params);
}

F1_0_PiPlusPiZeroPiZero::~F1_0_PiPlusPiZeroPiZero() {
  if (p_a1s)  { delete p_a1s;  p_a1s  = NULL; }
  if (p_rhos) { delete p_rhos; p_rhos = NULL; }
  if (p_TK1a)         { delete p_TK1a;         p_TK1a         = NULL; }
  if (p_TKstar1_pi0K) { delete p_TKstar1_pi0K; p_TKstar1_pi0K = NULL; }
}

void F1_0_PiPlusPiZeroPiZero::Construct(const FF_Parameters & params) {
  switch (m_mode) {
  case FF_0_PPP_mode::piP_pi0_pi0:
  case FF_0_PPP_mode::piM_piP_piP: Construct_3Pi(params);      break;
  case FF_0_PPP_mode::pi0_pi0_KM:  Construct_Pi0Pi0K(params);  break;
  default: break;
  }
  // Diagnostic dump (mirrors FF_0_PP.C's request #1): print whichever
  // propagator structures this instance actually built.
  std::string label = std::string("F1_0_PiPlusPiZeroPiZero, mode=")+
                       std::to_string(int(m_mode))+
                       (m_isF2 ? " (F2)" : " (F1)");
  DumpPropagatorStructure(label+" [a1/K1a]", int(m_ffmodel),
			   p_a1s!=NULL ? p_a1s : p_TK1a);
  DumpPropagatorStructure(label+" [rho/K*]", int(m_ffmodel),
			   p_rhos!=NULL ? p_rhos : p_TKstar1_pi0K);
}

void F1_0_PiPlusPiZeroPiZero::Construct_3Pi(const FF_Parameters & params) {
  m_norm = Complex(0., -((2.*sqrt(2)*(*params.p_model)("Vud", Tools::Vud)) /
			  (3.*m_fpi) ));
  if (m_ffmodel==ff_model::KS) {
    // Unified nomenclature (see FF_0_PP.C): "gamma" = rho' weight,
    // "delta" = rho'' weight (0 by default - the original KS 3pi model
    // has no third resonance), "alpha" = rho-omega mixing weight
    // (piM_piP_piP only, Construct() below skips it for piP_pi0_pi0).
    m_gamma = ReadComplexParam(params.p_model,"gammaMag_3pi",-0.145,"gammaPhase_3pi");
    m_delta = ReadComplexParam(params.p_model,"deltaMag_3pi", 0.,   "deltaPhase_3pi");
    m_alpha = ReadComplexParam(params.p_model,"alphaMag_3pi", 0.00185,"alphaPhase_3pi");
  }
  else if (m_ffmodel==ff_model::KS_CLEO) {
    // CLEO/default-TAUOLA 3pi current - see the class-level comment
    // above for the F1/F2/F3 architecture caveat. m_norm above (the
    // KS90 value, i*2sqrt2*Vud/(3 fpi)) is REUSED here rather than
    // overridden: the note writes the CLEO F1/F2/F3 as instances of
    // the SAME master PPPcurrent equation (which carries the overall
    // N), and gives no separate explicit normalization for this
    // variant - if the CLEO beta-fit was actually performed against a
    // different absolute normalization convention, the lineshape
    // SHAPE here is still correct but the absolute rate may not be;
    // flagged rather than silently assumed exactly right.
    Construct_3Pi_CLEO(params);
  }
  else if (m_ffmodel==ff_model::RChL2012) {
    // m_norm above is the KS normalisation (i*2sqrt2*Vud/(3 fpi));
    // RChL2012 uses its own normalisation N=cos(theta_C)/F, Eq.(3.1)
    // of 1509.09140 / just below Eq.(1) of 1203.3955 - overwrite it.
    m_norm       = (*params.p_model)("Vud",Tools::Vud) /
                   (*params.p_model)("F_rchl3pi",0.091337);
    m_F_rchl     = (*params.p_model)("F_rchl3pi",0.091337);
    m_FV_rchl    = (*params.p_model)("FV_rchl3pi",0.168652);
    m_FA_rchl    = (*params.p_model)("FA_rchl3pi",0.131425);
    m_GV_rchl    = sqr(m_F_rchl)/m_FV_rchl; // GV=F^2/FV constraint
    m_Mrho_rchl  = (*params.p_model)("Mrho_rchl3pi",0.771849);
    m_Mrhop_rchl = (*params.p_model)("Mrhop_rchl3pi",1.350000);
    m_Grhop_rchl = (*params.p_model)("Grhop_rchl3pi",0.448379);
    m_betarhop   = (*params.p_model)("betarhop_rchl3pi",-0.318551);
    m_Ma1_rchl   = (*params.p_model)("Ma1_rchl3pi",1.091865);
    m_Msigma     = (*params.p_model)("Msigma_rchl3pi",0.487512);
    m_Gsigma0    = (*params.p_model)("Gsigma_rchl3pi",0.700000);
    m_Rsigma2    = sqr((*params.p_model)("Rsigma_rchl3pi",1.866913));
    m_mpi2_iso   = sqr((Flavour(kf_pi).HadMass()*1.+
			2.*Flavour(kf_pi_plus).HadMass())/3.);
    m_lambda_p   = RChL::Lambda_p(m_F_rchl,m_FA_rchl,m_GV_rchl);
    m_lambda_pp  = RChL::Lambda_pp(m_F_rchl,m_GV_rchl,m_lambda_p);
    if (m_mode==FF_0_PPP_mode::piM_piP_piP) {
      // Table I of 1310.1053: default (unscaled) sigma couplings.
      m_alphasigma = (*params.p_model)("alphasigma_rchl3pi",-8.795938);
      m_betasigma  = (*params.p_model)("betasigma_rchl3pi", 9.763701);
      m_gammasigma = (*params.p_model)("gammasigma_rchl3pi",1.264263);
      m_deltasigma = (*params.p_model)("deltasigma_rchl3pi",0.656762);
    }
    else {
      // pi0 pi0 pi-: Eq.(16) of 1310.1053 - only alpha0_sigma,
      // gamma0_sigma are used (single sigma term in s3, see FF_RChL2012),
      // with the 0.63 scaling factor already folded into these defaults.
      m_alphasigma = (*params.p_model)("alpha0sigma_rchl3pi",0.63*1.139486);
      m_gammasigma = (*params.p_model)("gamma0sigma_rchl3pi",0.63*0.889769);
      m_betasigma  = 0.; m_deltasigma = 0.;
      m_Msigma     = 0.550; m_Gsigma0 = 0.700;
      m_Rsigma2    = sqr(0.000013);
    }
  }

  if (m_ffmodel==ff_model::KS) {
    if (m_mode==FF_0_PPP_mode::piP_pi0_pi0) {
      Propagator_Base * a11260  =
	new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
      Propagator_Base * rho770  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
      Propagator_Base * rho1450 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      p_a1s  = new Summed_Propagator();
      p_a1s->Add(a11260,  Complex(1., 0.));
      p_rhos = new Summed_Propagator();
      p_rhos->Add(rho770,  Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
    }
    else if (m_mode==FF_0_PPP_mode::piM_piP_piP) {
      Propagator_Base * a11260  =
	new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
      Propagator_Base * rho770  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * rho770_1  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * omega782  =
	new BreitWigner(LineShapes->Get(Flavour(kf_omega_782)));
      Propagator_Base * rho1450 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      p_a1s   = new Summed_Propagator();
      p_a1s->Add(a11260,  Complex(1., 0.));
      Summed_Propagator     * rho   = new Summed_Propagator();
      rho->Add(rho770,  Complex(1.,0.));
      Multiplied_Propagator * inter = new Multiplied_Propagator(); 
      inter->Add(rho770_1, Complex(1.,0.));
      inter->Add(omega782, Complex(1.,0.));
      rho->Add(inter,  m_alpha);
      p_rhos  = new Summed_Propagator();
      p_rhos->Add(rho,     Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
    }
  }
  else if (m_ffmodel==ff_model::RChL2012) {
    // Reuses the SAME p_rhos member the KS branch above populates,
    // built the same way (RChL_BW propagators inside a
    // Summed_Propagator), so FF_RChL2012 below can just call
    // (*p_rhos)(s) exactly like FF_KS does. Summed_Propagator's
    // operator() computes [Add(rho,w1)*rho(s) + Add(rhop,w2)*rhop(s)]
    // / (w1+w2) - with w1=1, w2=betarhop this IS Eq.(11)'s combined
    // rho+rho' propagator, no extra plumbing needed. RChL_BW (see
    // Propagator.[C,H]) matches the literal denominator sign/scaling
    // printed in Eq.(6) of 1203.3955, which differs from BreitWigner's
    // own convention (used by the KS branch above) - that's why this
    // branch can't just reuse plain BreitWigner objects.
    //
    // Charged/neutral rho choice matches the KS branch above: for
    // pi0 pi0 pi-, the rho in each pi0-pi- sub-system is CHARGED
    // (rho- -> pi0 pi-); for pi- pi- pi+, the rho in the pi-pi+
    // sub-system is NEUTRAL (rho0 -> pi+ pi-, hence also mixing with
    // omega there in the KS branch).
    Flavour rhoFlav = (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ?
		      Flavour(kf_rho_770_plus) : Flavour(kf_rho_770));
    Propagator_Base * rho  = new RChL_BW(LineShapes->Get(rhoFlav));
    Propagator_Base * rhop = new RChL_BW(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    p_rhos = new Summed_Propagator();
    p_rhos->Add(rho,  Complex(1.,0.));
    p_rhos->Add(rhop, Complex(m_betarhop,0.));
  }
}

void F1_0_PiPlusPiZeroPiZero::Construct_Pi0Pi0K(const FF_Parameters & params) {
  // Finkemeier-Mirkes hep-ph/9503474 (FM95), Tab.I row 6:
  // A^(abc) = sin(theta_c)/4, G1=T_K1a(Q^2)T_K*^(1)(s2),
  // G2=T_K1a(Q^2)T_K*^(1)(s1). Eq.(23)/(24): F_i = (2sqrt2 A/3fpi) G_i,
  // Cabibbo-SUPPRESSED (sin, not cos) since this is a |Delta S|=1
  // (one-kaon) channel - use Vus, not Vud.
  if (m_ffmodel!=ff_model::KS) return;
  m_norm = Complex((2.*sqrt(2.)*(*params.p_model)("Vus",Tools::Vus)*
		    (*params.p_model)("sinThetaC_over_4_num",1.)/4.) /
		    (3.*m_fpi), 0.);
  // xi: relative K1(1270) admixture in T_K1^(a), Eq.(32)-(33) of
  // FM95: |xi|=0.33, sign preferred by data is xi=+0.33 (Sec.VII).
  // xi_K1 weights K1(1270) against K1(1400) in T_K1^(a), Eq.(9):
  //   T_K1^(a) = BW_K1(1400) + xi_K1 BW_K1(1270).
  // FM95 quote a real 0.33, but that was fitted alongside THEIR K1
  // line shapes.  Exposed as a complex here for two reasons: the two
  // physical K1 states are K1A/K1B (3P1/1P1) mixtures, so a relative
  // PHASE between them is physically expected and a real coefficient
  // cannot produce it; and the magnitude needs re-deriving whenever the
  // K1 line shapes change - as they did when K1(1270) was given its
  // PDG channels (see K1_Decays.C).
  //
  // The default is MODEL-DEPENDENT.  FM95's 0.33 is not a fit: Eq.(33)
  // DERIVES |xi| = 0.33 from Gamma(K1(1270)->K*pi)/Gamma(K1(1400)->K*pi) after
  // phase-space correction, against a fixed-width K1(1270).  Model 106 replaces
  // that K1(1270) with a Flatte, so the ratio the derivation rests on no longer
  // holds and xi has to be re-determined; 0.70 is the measured optimum against
  // BaBar 2007 / Belle 2010 / CLEO 2000.  Carrying FM95's 0.33 into 106 would
  // silently give a much worse model (chi2/ndf 40.9 vs 10.6 at 150k).
  m_xiK1 = ReadComplexParam(params.p_model,"xiK1Mag",
                            m_ffmodel==ff_model::KS_flatte ? 0.70 : 0.33,
                            "xiK1Phase");

  // T_K1^(a) = [BW_K1(1400) + xi*BW_K1(1270)] / (1+xi), Eq.(32).
  // K1(1270)/K1(1400) lineshapes built and registered (K1_Decays.H/.C,
  // kf 10313/10323 and 20313/20323, confirmed against the real
  // ATOOLS/Phys/Flavour_Tags.H), using the proper running widths from
  // their rho-K/K*-pi decay channels rather than FM95's own simple
  // fixed-width prescription - a reasonable, arguably better-motivated
  // substitution.
  Total_Width_Base * wK11270 = LineShapes->Get(Flavour(kf_K_1_1270_plus));
  Total_Width_Base * wK11400 = LineShapes->Get(Flavour(kf_K_1_1400_plus));
  if (wK11270==NULL || wK11400==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": missing K1(1270)/K1(1400) "
	       <<"lineshape(s) for pi0_pi0_KM (FM95) - T_K1^(a) will be "
	       <<"treated as identically zero.\n";
  }
  else {
    Propagator_Base * K11400 = new BreitWigner(wK11400);
    Propagator_Base * K11270 = new BreitWigner(wK11270);
    p_TK1a = new Summed_Propagator();
    p_TK1a->Add(K11400, Complex(1.,0.));
    p_TK1a->Add(K11270, m_xiK1);
  }
  // T_K*^(1) = [BW_K*(892) + beta_K* BW_K*(1410)]/(1+beta_K*), Eq.(10).
  double betaKst = (*params.p_model)("betaKstar_pi0pi0K",-0.135);
  Propagator_Base * Kstar892  =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1_pi0K = new Summed_Propagator();
  p_TKstar1_pi0K->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1_pi0K->Add(Kstar1410, Complex(betaKst,0.));
}

///////////////////////////////////////////////////////////////////////////
//
// CLEO/default-TAUOLA 3pi current, tau_two_meson_currents_KS_RChiT.tex
// Sec."CLEO/default-TAUOLA 3pi current". See the class-level comment
// for the F1/F2/F3 architecture caveat this implementation relies on.
//
///////////////////////////////////////////////////////////////////////////

void F1_0_PiPlusPiZeroPiZero::Construct_3Pi_CLEO(const FF_Parameters & params) {
  // Nominal TAUOLA-CLEO masses/widths (Table "Complete coding inputs
  // for the default TAUOLA-CLEO 3pi current"). All overridable, all
  // DISTINCT from the KS90 rho/a1 numbers above (same physical
  // resonances, genuinely different fitted values in this current).
  m_MCLEOrho    = (*params.p_model)("MCLEOrho",   0.7743);
  m_GCLEOrho    = (*params.p_model)("GCLEOrho",   0.1491);
  m_MCLEOrhop   = (*params.p_model)("MCLEOrhop",  1.370);
  m_GCLEOrhop   = (*params.p_model)("GCLEOrhop",  0.386);
  m_MCLEOsigma  = (*params.p_model)("MCLEOsigma", 0.860);
  m_GCLEOsigma  = (*params.p_model)("GCLEOsigma", 0.880);
  m_MCLEOf0     = (*params.p_model)("MCLEOf0",    1.186);
  m_GCLEOf0     = (*params.p_model)("GCLEOf0",    0.350);
  m_MCLEOf2     = (*params.p_model)("MCLEOf2",    1.275);
  m_GCLEOf2     = (*params.p_model)("GCLEOf2",    0.185);
  m_MCLEOa1     = (*params.p_model)("MCLEOa1",    1.275);
  m_GCLEOa1     = (*params.p_model)("GCLEOa1",    0.700);
  // Dormant a1' term - beta=0 in the nominal current, so its own
  // mass/width are irrelevant unless overridden together with beta.
  m_betaCLEOa1p = ReadComplexParam(params.p_model,
				   "betaMag_3pi_CLEOa1prime",0.,"betaPhase_3pi_CLEOa1prime");
  m_MCLEOa1p    = (*params.p_model)("MCLEOa1prime", 1.275);
  m_GCLEOa1p    = (*params.p_model)("GCLEOa1prime", 0.700);
  m_MCLEOKstar  = (*params.p_model)("MCLEOKstar",   0.894);
  m_mCLEOK      = (*params.p_model)("mCLEOK",       0.496);
  // pi0/pi+- mass average as used in the note's own GEV inputs - using
  // the physical per-index mass in the actual formula evaluation
  // (m_masses2[...]) instead; this m_mCLEOpi is only used inside the
  // K*K threshold/momentum piece of WGA, which needs a single generic
  // pion mass scale, not per-index masses.
  m_mCLEOpi     = 0.137;
  // Complex fitted amplitudes, Eq.(eq:CLEO-betas) (nominal TAUOLA-CLEO
  // values - the note also quotes a slightly different "CLEO
  // published" mass/width set, not implemented here as a separate
  // option; flag if you want that added too).
  m_beta1 = Complex(1.,0.); // real, fixed
  m_beta2 = ReadComplexParam(params.p_model,"beta2Mag_3pi_CLEO",0.12,"beta2Phase_3pi_CLEO",0.99*M_PI);
  m_beta3 = ReadComplexParam(params.p_model,"beta3Mag_3pi_CLEO",0.37,"beta3Phase_3pi_CLEO",-0.15*M_PI);
  m_beta4 = ReadComplexParam(params.p_model,"beta4Mag_3pi_CLEO",0.87,"beta4Phase_3pi_CLEO",0.53*M_PI);
  m_beta5 = ReadComplexParam(params.p_model,"beta5Mag_3pi_CLEO",0.71,"beta5Phase_3pi_CLEO",0.56*M_PI);
  m_beta6 = ReadComplexParam(params.p_model,"beta6Mag_3pi_CLEO",2.10,"beta6Phase_3pi_CLEO",0.23*M_PI);
  m_beta7 = ReadComplexParam(params.p_model,"beta7Mag_3pi_CLEO",0.77,"beta7Phase_3pi_CLEO",-0.54*M_PI);
  msg_Out()<<"### Propagator structure for \"F1_0_PiPlusPiZeroPiZero (CLEO), "
	   <<"mode="<<int(m_mode)<<"\" (FORM_FACTOR = "<<int(m_ffmodel)<<"):\n"
	   <<"###   rho(CLEO): M = "<<m_MCLEOrho<<" GeV, Gamma = "<<m_GCLEOrho<<" GeV\n"
	   <<"###   rho'(CLEO): M = "<<m_MCLEOrhop<<" GeV, Gamma = "<<m_GCLEOrhop<<" GeV\n"
	   <<"###   sigma: M = "<<m_MCLEOsigma<<" GeV, Gamma = "<<m_GCLEOsigma<<" GeV\n"
	   <<"###   f0: M = "<<m_MCLEOf0<<" GeV, Gamma = "<<m_GCLEOf0<<" GeV\n"
	   <<"###   f2: M = "<<m_MCLEOf2<<" GeV, Gamma = "<<m_GCLEOf2<<" GeV\n"
	   <<"###   a1(CLEO): M = "<<m_MCLEOa1<<" GeV, Gamma = "<<m_GCLEOa1<<" GeV\n";
}

double F1_0_PiPlusPiZeroPiZero::
CLEO_P(const double & S,const double & m1,const double & m2) {
  double arg = (S-sqr(m1+m2))*(S-sqr(m1-m2));
  return (arg>0. && S>0. ? sqrt(arg)/sqrt(S) : 0.);
}

double F1_0_PiPlusPiZeroPiZero::
CLEO_GammaL(const double & S,const double & M,const double & Gamma,
	    const double & m1,const double & m2,int L) {
  double PS = CLEO_P(S,m1,m2), PM = CLEO_P(sqr(M),m1,m2);
  if (PM<=0. || S<=0.) return 0.;
  return Gamma*(M/sqrt(S))*pow(PS/PM, 2*L+1);
}

Complex F1_0_PiPlusPiZeroPiZero::
CLEO_BWL(const double & S,const double & M,const double & Gamma,
	 const double & m1,const double & m2,int L) {
  double G = CLEO_GammaL(S,M,Gamma,m1,m2,L);
  return sqr(M)/Complex(S-sqr(M), -M*G);
}

double F1_0_PiPlusPiZeroPiZero::CLEO_w1(const double & x) const {
  if (x<0.1753) return 0.;
  double dx = x-0.1753;
  if (x<0.823) return 5.809*pow(dx,3)*(1.-3.0098*dx+4.5792*pow(dx,3));
  return -13.914+27.679*x-13.393*sqr(x)+3.1924*pow(x,3)-0.10487*pow(x,4);
}

double F1_0_PiPlusPiZeroPiZero::CLEO_w2(const double & x) const {
  if (x<0.1676) return 0.;
  double dx = x-0.1676;
  if (x<0.823) return 6.2845*pow(dx,3)*(1.-2.9595*dx+4.3355*pow(dx,3));
  return -15.411+32.088*x-17.666*sqr(x)+4.9355*pow(x,3)-0.37498*pow(x,4);
}

double F1_0_PiPlusPiZeroPiZero::CLEO_wK(const double & x) const {
  if (x<=sqr(m_MCLEOKstar+m_mCLEOK)) return 0.;
  // p_{K*K}(x) = lambda^{1/2}(x,MK*^2,mK^2)/(2x), Eq.(eq:CLEO-pK).
  double MKst2 = sqr(m_MCLEOKstar), mK2 = sqr(m_mCLEOK);
  double lambda = sqr(x-MKst2-mK2)-4.*MKst2*mK2;
  return (lambda>0. ? sqrt(lambda)/(2.*x) : 0.);
}

double F1_0_PiPlusPiZeroPiZero::CLEO_WGA(const double & Q2) const {
  double C3pi = sqr(0.2384), CKst = sqr(4.7621)*C3pi;
  return C3pi*(CLEO_w1(Q2)+CLEO_w2(Q2)) + CKst*CLEO_wK(Q2);
}

Complex F1_0_PiPlusPiZeroPiZero::CLEO_Aa1(const double & Q2) const {
  // Eq.(eq:CLEO-Aa1). The 1.3281*0.806 normalization constant in the
  // denominator is transcribed literally from the note - it is an
  // internal TAUOLA normalization of WGA relative to the pole width,
  // not independently re-derived here.
  Complex D1 = Complex(Q2-sqr(m_MCLEOa1),
		       -m_MCLEOa1*m_GCLEOa1/(1.3281*0.806)*CLEO_WGA(Q2));
  Complex term1 = sqr(m_MCLEOa1)/D1;
  Complex term2(0.,0.);
  if (m_betaCLEOa1p!=Complex(0.,0.)) {
    Complex D2 = Complex(Q2-sqr(m_MCLEOa1p),
			 -m_MCLEOa1p*m_GCLEOa1p/(1.3281*0.806)*CLEO_WGA(Q2));
    term2 = m_betaCLEOa1p*sqr(m_MCLEOa1p)/D2;
  }
  return term1+term2;
}

Complex F1_0_PiPlusPiZeroPiZero::
CLEO_F1term(const double & Q2,const double & s1,const double & s2,
	    const double & s3,const double & m1_2,const double & m2_2,
	    const double & m3_2) const {
  // Eq.(eq:CLEO-F1). FIX: the BW_L "daughter masses appropriate to the
  // subchannel" (per the note's own BW_L(s;R) shorthand) are NOT the
  // same pion mass for every term - s1=(p2+p3)^2 is formed by pions
  // 2,3 (so its BW_L needs masses m2,m3), s2=(p1+p3)^2 by pions 1,3,
  // s3=(p1+p2)^2 by pions 1,2. An earlier version used sqrt(m1_2)
  // uniformly everywhere - harmless for pi-pi-pi+ (all three masses
  // equal) but WRONG for pi0pi0pi- (m1=m2=mpi0 differs from m3=mpi+-),
  // silently mixing up which daughter-mass pair enters each BW_L call.
  double m1 = sqrt(m1_2), m2 = sqrt(m2_2), m3 = sqrt(m3_2);
  Complex Aa1 = CLEO_Aa1(Q2);
  Complex bracket =
    m_beta1*CLEO_BWL(s1,m_MCLEOrho, m_GCLEOrho, m2,m3,1) +
    m_beta2*CLEO_BWL(s1,m_MCLEOrhop,m_GCLEOrhop,m2,m3,1)
    - m_beta3*((s3-m3_2)-(s1-m1_2))/3. * CLEO_BWL(s2,m_MCLEOrho, m_GCLEOrho, m1,m3,1)
    - m_beta4*((s3-m3_2)-(s1-m1_2))/3. * CLEO_BWL(s2,m_MCLEOrhop,m_GCLEOrhop,m1,m3,1)
    + m_beta5*(Q2+s3-m2_2)*(2.*m3_2+2.*m1_2-s3)/(18.*s3) *
      CLEO_BWL(s3,m_MCLEOf2,m_GCLEOf2,m1,m2,2)
    + (2./3.)*m_beta6*CLEO_BWL(s3,m_MCLEOsigma,m_GCLEOsigma,m1,m2,0)
    + (2./3.)*m_beta7*CLEO_BWL(s3,m_MCLEOf0,m_GCLEOf0,m1,m2,0);
  return Aa1*bracket;
}

Complex F1_0_PiPlusPiZeroPiZero::
CLEO_F3term(const double & Q2,const double & s1,const double & s2,
	    const double & s3,const double & m1_2,const double & m2_2,
	    const double & m3_2) const {
  // Eq.(eq:CLEO-F3) - CLEO's own "F3" (third transverse axial basis
  // coefficient), NOT the anomalous vector form factor. See the
  // class-level comment for how this combines with F1/F2 here. Same
  // per-subchannel daughter-mass fix as CLEO_F1term above.
  double m1 = sqrt(m1_2), m2 = sqrt(m2_2), m3 = sqrt(m3_2);
  Complex Aa1 = CLEO_Aa1(Q2);
  Complex bracket =
    (m_beta3/3.)*((s2-m2_2)-(s3-m3_2))*CLEO_BWL(s1,m_MCLEOrho, m_GCLEOrho, m2,m3,1)
    + (m_beta3/3.)*((s3-m3_2)-(s1-m1_2))*CLEO_BWL(s2,m_MCLEOrho, m_GCLEOrho, m1,m3,1)
    + (m_beta4/3.)*((s2-m2_2)-(s3-m3_2))*CLEO_BWL(s1,m_MCLEOrhop,m_GCLEOrhop,m2,m3,1)
    + (m_beta4/3.)*((s3-m3_2)-(s1-m1_2))*CLEO_BWL(s2,m_MCLEOrhop,m_GCLEOrhop,m1,m3,1)
    - (m_beta5/2.)*((s1-m1_2)-(s2-m2_2))*CLEO_BWL(s3,m_MCLEOf2,m_GCLEOf2,m1,m2,2);
  return Aa1*bracket;
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_KS(const double & s123,const double & s12,const double & s13) {
  // Base-class convention (see FF_0_PPP.H): s12=(p1+p2)^2=paper's s3
  // (the non-resonant like-sign pair), s13=(p1+p3)^2=paper's s2. The
  // paper's genuine s1=(p2+p3)^2 - needed for F2=F(s2,s1,Q^2)=
  // N*BWa1(Q^2)*Brho(s1), Eq.(3.5) - is NOT one of the two base-class
  // arguments and must be reconstructed via momentum conservation
  // (s1+s2+s3=Q^2+m1^2+m2^2+m3^2). An earlier version of this
  // function used s12 (paper's s3) directly for F2, which puts the
  // rho propagator on the never-resonant like-sign pion pair instead
  // of the correct opposite-sign combination - confirmed wrong by
  // explicit numerical check against Eq.(3.3)-(3.5).
  if (m_ffmodel==ff_model::KS_CLEO &&
      (m_mode==FF_0_PPP_mode::piP_pi0_pi0 || m_mode==FF_0_PPP_mode::piM_piP_piP)) {
    double m1_2 = m_masses2[m_pi[0]], m2_2 = m_masses2[m_pi[1]], m3_2 = m_masses2[m_pi[2]];
    double paper_s2 = s13, paper_s3 = s12;
    double paper_s1 = s123 - s12 - s13 + m1_2 + m2_2 + m3_2;
    // R_3pi (+1 charged/-1 neutral per the note) is an overall sign
    // common to every form factor of a given channel - physically
    // immaterial for a standalone channel (cancels in every rate
    // bilinear; the note itself says so explicitly for an analogous
    // overall-sign choice) and omitted here, consistent with the
    // existing KS90 branch above, which also does not apply it.
    // The earlier CAVEAT above is now RESOLVED by the note, in two ways.
    //
    // (1) No F3 folding. The note's "Active transverse basis in the
    // default TAUOLA-CLEO current" gives the active coefficients
    // explicitly as c1=+2sqrt2/3, c2=-2sqrt2/3, c3=0 (Eq.(eq:CLEO-
    // active-coefficients)), so the two-slot map is simply
    // F1_eff=F1^CLEO, F2_eff=F2^CLEO (Eq.(eq:CLEO-two-slot-map)). The
    // generic reduction F_i - (c3/c_i)F3 (Eq.(eq:CLEO-generic-fold))
    // has a vanishing fold term here. The note's "Implementation
    // warning" states outright that c1=c2=c3=1 is NOT the default
    // TAUOLA-CLEO convention, which is what the previous version
    // assumed. CLEO_F3term is deliberately kept below: it is correct
    // as written and is needed to compare against a redundant
    // three-vector basis, it just has coefficient zero here.
    //
    // (2) Slot assignment. The note's F1^CLEO leads with BW(s1;rho),
    // s1=(p2+p3)^2, and multiplies (p2-p3) - which is v2 in
    // VA_0_PiPiPi::Calc, not v1. Since this class's F1 slot multiplies
    // v1=T(p1-p3), that slot must carry the (1,3) subsystem, i.e. the
    // note's F2=F1(1<->2). This is exactly the pairing the KS90 branch
    // below already uses (F1 -> rhos(s13), F2 -> rhos(paper_s1)); the
    // previous version had the two orderings the other way round, so
    // the CLEO and KS branches of the same class disagreed and only KS
    // matched the basis vectors. No extra sign is needed: c2=-c1 is
    // already absorbed by building v1 as (p1-p3) rather than (p3-p1).
    if (!m_isF2) {
      return CLEO_F1term(s123, paper_s2, paper_s1, paper_s3, m2_2, m1_2, m3_2);
    }
    else {
      return CLEO_F1term(s123, paper_s1, paper_s2, paper_s3, m1_2, m2_2, m3_2);
    }
  }
  if (m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // FM95 Tab.I row 6: G1(Q^2,s2,s3)=T_K1a(Q^2)T_K*^(1)(s2),
    // G2(Q^2,s1,s3)=T_K1a(Q^2)T_K*^(1)(s1) - same s1-reconstruction
    // as the pion modes above (s2=s13, s1 needs momentum conservation).
    if (p_TKstar1_pi0K==NULL) return Complex(0.,0.);
    Complex TK1a = (p_TK1a!=NULL ? (*p_TK1a)(s123) : Complex(0.,0.));
    if (m_isF2) {
      double s1 = s123 - s12 - s13 +
	m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
      return TK1a * (*p_TKstar1_pi0K)(s1);
    }
    return TK1a * (*p_TKstar1_pi0K)(s13);
  }
  if (p_a1s==NULL || p_rhos==NULL) return Complex(0.,0.);
  // NOTE: do NOT multiply by m_norm here. FF_0_PPP_Base::operator()
  // already applies it ("case ff_model::KS: return m_norm*FF_KS(...)"),
  // so every FF_XXX() override must return the DIMENSIONLESS shape only -
  // as the base-class fallbacks state ("F=1 times whatever is already in
  // m_norm") and as every other derived class in this file does.
  // This function used to apply m_norm a second time, squaring it in the
  // rate: for 3pi, |m_norm| = 2sqrt2*Vud/(3 fpi) = 9.93, so the partial
  // width came out |m_norm|^2 = 98.7 times too large (measured 96.6/94.1
  // times the PDG value for pi-pi0pi0 / pi-pi-pi+ before this fix).
  // The pi0_pi0_KM branch above had the same double-application, with
  // |m_norm| = 0.573 there, so that channel was instead a factor
  // |m_norm|^2 = 0.328 too SMALL - its apparent agreement with the PDG
  // rate was an artefact of this bug, not evidence the mode was right.
  if (m_isF2) {
    double s1 = s123 - s12 - s13 +
      m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
    return (*p_a1s)(s123) * (*p_rhos)(s1);
  }
  return (*p_a1s)(s123) * (*p_rhos)(s13);
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_RChiPT(const double & s123,const double & s1,const double & s2) {
  // RChiPT was never derived/implemented for this channel (a genuine
  // gap, not a physics-forced zero) - fall back to the constant
  // (isospin/CKM-only) form factor per the general policy.
  return Complex(1.,0.);
}

///////////////////////////////////////////////////////////////////////////
//
// RChL 3pi axial form factor, Eqs.(4),(6),(8) of 1203.3955, plus the
// sigma-meson extension Eqs.(3)-(7) of 1310.1053.
// Base-class arg convention: s123=q^2, s1_arg=(p1+p2)^2=s3(paper),
// s2_arg=(p1+p3)^2=s2(paper) - see FF_0_PPP.H comment. s1(paper) is
// reconstructed below. F2 follows from Bose symmetry, Eq below (10) of
// 1203.3955: F2(q2,s2,s1)=F1(q2,s1,s2) - i.e. F2 is obtained by calling
// this same code with the paper's s1,s2 swapped.
//
///////////////////////////////////////////////////////////////////////////

Complex F1_0_PiPlusPiZeroPiZero::
FF_RChL2012(const double & s123,const double & s1_arg,const double & s2_arg) {
  double q2 = s123, s3 = s1_arg, s2 = s2_arg;
  double s1 = q2 - s2 - s3 + 3.*m_mpi2_iso; // m1^2=m2^2=m3^2=mpi2 here
  if (m_isF2) std::swap(s1,s2); // Bose symmetry
  double R3pi = (m_mode==FF_0_PPP_mode::piM_piP_piP ? 1. : -1.);

  if (p_rhos==NULL) return Complex(0.,0.);
  Complex propS1 = (*p_rhos)(s1), propS2 = (*p_rhos)(s2);

  // Fchi_1
  Complex Fchi1(-2.*sqrt(2.)/3.,0.);

  // FR_1, Eq.(6)
  Complex FR1 = (sqrt(2.)*m_FV_rchl*m_GV_rchl)/(3.*sqr(m_F_rchl)) *
    ( 3.*s1*propS1 -
      (2.*m_GV_rchl/m_FV_rchl-1.) *
      ( (2.*q2-2.*s1-s3)*propS1 + (s3-s1)*propS2 ) );

  // FRR_1, Eq.(6)
  double Ha = RChL::H(s1/q2,m_mpi2_iso/q2,m_lambda_p,m_lambda_pp);
  double Hb = RChL::H(s2/q2,m_mpi2_iso/q2,m_lambda_p,m_lambda_pp);
  Complex a1prop = q2/Complex(q2-sqr(m_Ma1_rchl),
			      -m_Ma1_rchl*RChL::Gamma_a1_PionFit(q2,m_Mrho_rchl,m_mpi2_iso));
  Complex FRR1 = (4.*m_FA_rchl*m_GV_rchl)/(3.*sqr(m_F_rchl)) * a1prop *
    ( ( -(m_lambda_p+m_lambda_pp)*3.*s1 + Ha*(2.*q2+s1-s3) ) * propS1 +
      Hb*(s3-s1)*propS2 );

  // sigma extension, Eqs.(3)-(4) [pi-pi-pi+] or Eq.(15) [pi0pi0pi-].
  double GsigS1 = RChL::GammaSigma(s1,sqr(m_Msigma),m_Gsigma0,m_mpi2_iso);
  double GsigS2 = RChL::GammaSigma(s2,sqr(m_Msigma),m_Gsigma0,m_mpi2_iso);
  double GsigS3 = RChL::GammaSigma(s3,sqr(m_Msigma),m_Gsigma0,m_mpi2_iso);
  Complex sigmaR, sigmaRR;
  if (m_mode==FF_0_PPP_mode::piM_piP_piP) {
    sigmaR  = m_alphasigma*RChL::BWsigma(s1,sqr(m_Msigma),m_Msigma,GsigS1)*
                          RChL::Fsigma(q2,s1,m_mpi2_iso,m_Rsigma2)
            + m_betasigma*RChL::BWsigma(s2,sqr(m_Msigma),m_Msigma,GsigS2)*
                          RChL::Fsigma(q2,s2,m_mpi2_iso,m_Rsigma2);
    sigmaRR = m_gammasigma*RChL::BWsigma(s1,sqr(m_Msigma),m_Msigma,GsigS1)*
                          RChL::Fsigma(q2,s1,m_mpi2_iso,m_Rsigma2)
            + m_deltasigma*RChL::BWsigma(s2,sqr(m_Msigma),m_Msigma,GsigS2)*
                          RChL::Fsigma(q2,s2,m_mpi2_iso,m_Rsigma2);
  }
  else {
    sigmaR  = m_alphasigma*RChL::BWsigma(s3,sqr(m_Msigma),m_Msigma,GsigS3)*
                          RChL::Fsigma(q2,s3,m_mpi2_iso,m_Rsigma2);
    sigmaRR = m_gammasigma*RChL::BWsigma(s3,sqr(m_Msigma),m_Msigma,GsigS3)*
                          RChL::Fsigma(q2,s3,m_mpi2_iso,m_Rsigma2);
  }
  FR1  += (sqrt(2.)*m_FV_rchl*m_GV_rchl)/(3.*sqr(m_F_rchl)) * sigmaR;
  FRR1 += (4.*m_FA_rchl*m_GV_rchl)/(3.*sqr(m_F_rchl)) * a1prop * sigmaRR;

  return (Fchi1+FR1+FRR1) * R3pi;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class F3_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  Complex m_gamma, m_delta;

  // RChL2012: this class computes the paper's F4 (pseudoscalar form
  // factor, Eq.9-10 of 1203.3955) - the "F3_0_PPP" name here is
  // legacy KS naming (cf. VA_0_PiPiPi.H's p_f3 comment), unrelated to
  // the paper's own F3 (which is fixed to 0 by construction, Sec.2.1).
  double  m_F_rchl, m_FV_rchl, m_GV_rchl, m_Mrho_rchl, m_mpi2_iso;
  double  m_Mrhop_rchl, m_betarhop;
  // Combined rho+rho' propagator (Eq.11), same RChL_BW+Summed_Propagator
  // construction as F1_0_PiPlusPiZeroPiZero's p_rhos - see that class's
  // Construct() comment. alpha_2 (Eq.10) gets the same Eq.(11)
  // substitution as F1^R/F1^RR per the paper's "insert Eq.(11) into
  // our Eqs.(6) and (10)" - this was a plain (rho-only, wrong-sign)
  // propagator in an earlier version of this class; fixed here.
  Propagator_Base * p_rho_combo;
  
  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChL2012(const double & s123,const double & s1,const double & s2);
public :
  F3_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
    p_rho_combo(NULL)
  {
    Construct(params);
  }
  ~F3_0_PiPlusPiZeroPiZero() {
    if (p_rho_combo) { delete p_rho_combo; p_rho_combo = NULL; }
  }
};

void F3_0_PiPlusPiZeroPiZero::Construct(const FF_Parameters & params) {
  if (m_ffmodel==ff_model::RChL2012 &&
      (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ||
       m_mode==FF_0_PPP_mode::piM_piP_piP)) {
    m_norm       = (*params.p_model)("Vud",Tools::Vud) /
                   (*params.p_model)("F_rchl3pi",0.091337);
    m_F_rchl     = (*params.p_model)("F_rchl3pi",0.091337);
    m_FV_rchl    = (*params.p_model)("FV_rchl3pi",0.168652);
    m_GV_rchl    = sqr(m_F_rchl)/m_FV_rchl;
    m_Mrho_rchl  = (*params.p_model)("Mrho_rchl3pi",0.771849);
    m_Mrhop_rchl = (*params.p_model)("Mrhop_rchl3pi",1.350000);
    m_betarhop   = (*params.p_model)("betarhop_rchl3pi",-0.318551);
    m_mpi2_iso   = sqr((Flavour(kf_pi).HadMass()+
			2.*Flavour(kf_pi_plus).HadMass())/3.);
    // Same charged/neutral rho fix as F1_0_PiPlusPiZeroPiZero::Construct()
    // above - see that comment for the physics reasoning.
    Flavour rhoFlav = (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ?
		      Flavour(kf_rho_770_plus) : Flavour(kf_rho_770));
    Propagator_Base * rho  = new RChL_BW(LineShapes->Get(rhoFlav));
    Propagator_Base * rhop = new RChL_BW(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    Summed_Propagator * combo = new Summed_Propagator();
    combo->Add(rho,  Complex(1.,0.));
    combo->Add(rhop, Complex(m_betarhop,0.));
    p_rho_combo = combo;
  }
  DumpPropagatorStructure(std::string("F3_0_PiPlusPiZeroPiZero, mode=")+
			   std::to_string(int(m_mode)), int(m_ffmodel), p_rho_combo);
}

///////////////////////////////////////////////////////////////////////////
//
// F4 of Eq.(9)-(10) of 1203.3955 (pseudoscalar form factor). Same arg
// convention / s1-reconstruction as F1_0_PiPlusPiZeroPiZero::FF_RChL2012
// above - see that function's comment.
//
///////////////////////////////////////////////////////////////////////////
Complex F3_0_PiPlusPiZeroPiZero::
FF_RChL2012(const double & s123,const double & s1_arg,const double & s2_arg) {
  if (p_rho_combo==NULL) return Complex(0.,0.);
  double q2 = s123, s3 = s1_arg, s2 = s2_arg;
  double s1 = q2 - s2 - s3 + 3.*m_mpi2_iso;
  double R3pi = (m_mode==FF_0_PPP_mode::piM_piP_piP ? 1. : -1.);
  double kappa = (m_mode==FF_0_PPP_mode::piM_piP_piP ? 1. : 0.5);

  double Fchi4 = (2.*sqrt(2.)/3.) * m_mpi2_iso *
    ( 3.*(s3-m_mpi2_iso) - q2*(1.+2.*kappa*R3pi) ) /
    ( 2.*q2*(q2-m_mpi2_iso) );

  Complex propS1 = (*p_rho_combo)(s1);
  Complex propS2 = (*p_rho_combo)(s2);
  // alpha_2(q2,s1,s2) = (3GV/FV)(s1/q2)(mpi2/(q2-mpi2)) (s3-s2)/propS1^-1
  Complex alpha2_21 = (3.*m_GV_rchl/m_FV_rchl)*(s2/q2)*(m_mpi2_iso/(q2-m_mpi2_iso))*
                      (s3-s1) * propS2;
  Complex alpha2_12 = (3.*m_GV_rchl/m_FV_rchl)*(s1/q2)*(m_mpi2_iso/(q2-m_mpi2_iso))*
                      (s3-s2) * propS1;
  Complex FR4 = -(sqrt(2.)*m_FV_rchl*m_GV_rchl)/(3.*sqr(m_F_rchl)) *
                (alpha2_21+alpha2_12);

  return (Complex(Fchi4,0.)+FR4) * R3pi;
}

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  Complex m_gamma, m_delta;

  // pi0_pi0_KM only (FM95 Tab.II row 6): vector ("anomaly") form
  // factor is IDENTICALLY ZERO for piP_pi0_pi0/piM_piP_piP (G-parity
  // of the pion system - see the original Kuhn-Santamaria discussion),
  // but non-zero once a kaon replaces one pion. Needs T_K*^(2), the
  // THREE-resonance K*(892)+lambda K*(1410)+mu K*(1714) mix, Eq.(42),
  // lambda=6.5/(-26)=-0.25, mu=1/(-26)=-0.038 (same lambda,mu as the
  // rho sector, Sec.V). K*(1680)/K*'' now built and registered
  // (Kstar_Decays.H/.C) - see Construct() for the remaining
  // mass-convention caveat (1.700/0.235 GeV here vs the FM95 Kpi
  // CLEO-tune's number, both attached to the same kf-code).
  Summed_Propagator * p_TKstar1_v, * p_TKstar2;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s1,const double & s2);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    // piP_pi0_pi0/piM_piP_piP: vector current forced to vanish by
    // G-parity - genuine physics, stays hard zero regardless of model.
    // pi0_pi0_KM: no RChiPT implementation exists (a genuine gap, not
    // a physics-forced zero) - falls back to the constant form factor.
    if (m_mode==FF_0_PPP_mode::pi0_pi0_KM) return Complex(1.,0.);
    return Complex(0.,0.);
  }
public :
  FS_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.)),
    p_TKstar1_v(NULL), p_TKstar2(NULL)
  {
    Construct(params);
  }
  ~FS_0_PiPlusPiZeroPiZero() {
    if (p_TKstar1_v) { delete p_TKstar1_v; p_TKstar1_v = NULL; }
    if (p_TKstar2)   { delete p_TKstar2;   p_TKstar2   = NULL; }
  }
};

void FS_0_PiPlusPiZeroPiZero::Construct(const FF_Parameters & params)  {
  if (m_ffmodel==ff_model::KS && m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // FM95 Tab.II row 6: A^(abc)=sin(theta_c), F3^(abc)=A/(2sqrt2 pi^2
    // fpi^3) G3, Eq.(25). Cabibbo-suppressed (Vus), |Delta S|=1.
    m_norm = Complex((*params.p_model)("Vus",Tools::Vus) /
		      (2.*sqrt(2.)*sqr(M_PI)*pow(m_fpi,3)), 0.);

    double betaKst = (*params.p_model)("betaKstar_pi0pi0K",-0.135);
    Propagator_Base * Kstar892_1 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
    Propagator_Base * Kstar1410_1 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    p_TKstar1_v = new Summed_Propagator();
    p_TKstar1_v->Add(Kstar892_1,  Complex(1., 0.));
    p_TKstar1_v->Add(Kstar1410_1, Complex(betaKst,0.));
    // T_K*^(2), Eq.(42): K*(892) + lambda*K*(1410) + mu*K*(1714)/K*(1680).
    // K*(1680) lineshape registered (Kstar_Decays.H/.C, kf 30313/30323).
    // Weights now overridable (unified nomenclature) rather than
    // hardcoded literals.
    Complex lambda = ReadComplexParam(params.p_model,
				      "gammaMag_pi0pi0K",-0.25,"gammaPhase_pi0pi0K");
    Complex mu     = ReadComplexParam(params.p_model,
				      "deltaMag_pi0pi0K",-0.038,"deltaPhase_pi0pi0K");
    Total_Width_Base * wKstarpp = LineShapes->Get(Flavour(kf_K_star_1680_plus));
    p_TKstar2 = new Summed_Propagator();
    Propagator_Base * Kstar892_2 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
    Propagator_Base * Kstar1410_2 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    p_TKstar2->Add(Kstar892_2,  Complex(1.,0.));
    p_TKstar2->Add(Kstar1410_2, lambda);
    if (wKstarpp==NULL) {
      msg_Error()<<"Error in "<<METHOD<<": missing K*(1714)/K*'' lineshape "
		 <<"for pi0_pi0_KM vector form factor (FM95 Eq.42) - "
		 <<"dropping the mu*BW_K*'' term (mu=-0.038, small but not "
		 <<"zero) until this resonance is added to Line_Shapes.\n";
    }
    else {
      Propagator_Base * Kstarpp = new BreitWigner(wKstarpp);
      p_TKstar2->Add(Kstarpp, mu);
    }
  }
  std::string label = std::string("FS_0_PiPlusPiZeroPiZero, mode=")+
                       std::to_string(int(m_mode));
  DumpPropagatorStructure(label+" [T_K*^(1)_v]", int(m_ffmodel), p_TKstar1_v);
  DumpPropagatorStructure(label+" [T_K*^(2)]",   int(m_ffmodel), p_TKstar2);
}

Complex FS_0_PiPlusPiZeroPiZero::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (m_mode!=FF_0_PPP_mode::pi0_pi0_KM) return Complex(0.,0.);
  if (p_TKstar1_v==NULL || p_TKstar2==NULL) return Complex(0.,0.);
  // G3(Q^2,s1,s2,s3) = (1/4) T_K*^(2)(Q^2) [T_K*^(1)(s1) - T_K*^(1)(s2)]
  // s2=s13 directly available; s1 reconstructed via momentum conservation.
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  return 0.25 * (*p_TKstar2)(s123) *
    ( (*p_TKstar1_v)(s1) - (*p_TKstar1_v)(s13) );
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// K^+ K^- pi^- form factors - STUBS ONLY, all return 0.
//
// Per arXiv:1509.09140 Sec.3.3, two parametrizations are eventually
// wanted here:
//  - CPC version [Jadach/Was/Decker/Kuhn CPC 76 (1993) 361]: dominant
//    production via a1 -> K* K and a1 -> rho pi for the axial-vector
//    form factors (F1, F2 here), rho' -> (rho pi; K* K) for the vector
//    form factor (F5, i.e. FS here); product of Breit-Wigner
//    amplitudes for each resonance.
//  - RChL: vector current from the Wess-Zumino term / odd-intrinsic-
//    parity amplitude [Dumm, Roig, Pich, Portoles, PRD 81 (2010)
//    034031]; direct-vertex + single-resonance + double-resonance
//    contributions.
// Neither is implemented yet - these classes only exist so VA_0_KKPi
// (see VA_0_KKPi.C) has something to link against while we fill in
// the actual physics.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F_0_KKPi_Stub : public FF_0_PPP_Base {
  void    Construct(const FF_Parameters & params) {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F_0_KKPi_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// K^+ K^- pi^- axial form factors F1, F2, RChL2012 (Eqs.12-14 of
// 1203.3955, Sec.2.2, "K- pi- K+" channel).
//
// INDEX MAPPING CAVEAT: the paper attributes momenta for this channel
// as (p1,p2,p3) = (K-, pi-, K+) (Table 1), whereas VA_0_KKPi.H fixes
// (p1,p2,p3) = (K+, K-, pi-). Consequently the base class's arguments
// (s123=q^2, s1_arg=(p1+p2)^2, s2_arg=(p1+p3)^2 in OUR ordering) relate
// to the PAPER's (s1,s2) as: paper_s2 = s1_arg, paper_s1 = s2_arg (the
// two are swapped relative to the 3pi case - see the derivation this
// class's FF_RChL2012 comment references). paper_s3 is unaffected by
// the swap since s1+s2 is invariant either way.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_KKPi : public FF_0_PPP_Base {
  bool    m_isF2;
  double  m_F_rchl, m_FV_rchl, m_FA_rchl, m_GV_rchl;
  double  m_Mrho_rchl, m_MKst_rchl, m_Ma1_rchl, m_mK2, m_mpi2;
  Total_Width_Base * p_GRho, * p_GKst;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
  Complex FF_RChL2012(const double & s123,const double & s1,const double & s2);
public:
  F1_0_KKPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), m_isF2(false), p_GRho(NULL), p_GKst(NULL)
  {
    if (params.m_name=="F2_0_KKPi") m_isF2 = true;
    Construct(params);
  }
};

void F1_0_KKPi::Construct(const FF_Parameters & params) {
  m_norm      = (*params.p_model)("Vud",Tools::Vud) /
                (*params.p_model)("F_rchlKKpi",0.0924); // even # of kaons: cos(theta_C)/F
  m_F_rchl    = (*params.p_model)("F_rchlKKpi",0.0924);
  m_FV_rchl   = (*params.p_model)("FV_rchlKKpi",0.18);
  m_FA_rchl   = (*params.p_model)("FA_rchlKKpi",0.149);
  m_GV_rchl   = sqr(m_F_rchl)/m_FV_rchl;
  m_Mrho_rchl = (*params.p_model)("Mrho_rchlKKpi",0.775);
  m_MKst_rchl = (*params.p_model)("MKst_rchlKKpi",0.8953);
  m_Ma1_rchl  = (*params.p_model)("Ma1_rchlKKpi",1.12);
  m_mK2       = sqr(Flavour(kf_K_plus).HadMass());
  m_mpi2      = sqr(Flavour(kf_pi_plus).HadMass());
  p_GRho = LineShapes->Get(Flavour(kf_rho_770_plus));
  p_GKst = LineShapes->Get(Flavour(kf_K_star_892_plus));
  // This class uses bare Total_Width_Base pointers directly in its
  // formula (not wrapped in a Summed_Propagator/Propagator_Base), so
  // DumpPropagatorStructure (which expects the latter) doesn't apply -
  // print the raw masses/on-shell widths directly instead.
  msg_Out()<<"### Propagator structure for \"F1_0_KKPi, mode="<<int(m_mode)
	   <<(m_isF2 ? " (F2)" : " (F1)")<<"\" (FORM_FACTOR = "<<int(m_ffmodel)<<"):\n";
  if (p_GRho!=NULL) msg_Out()<<"###   rho(770): M = "<<m_Mrho_rchl<<" GeV,  "
			     <<"Gamma(M^2) = "<<(*p_GRho)(sqr(m_Mrho_rchl))<<" GeV\n";
  if (p_GKst!=NULL) msg_Out()<<"###   K*(892): M = "<<m_MKst_rchl<<" GeV,  "
			     <<"Gamma(M^2) = "<<(*p_GKst)(sqr(m_MKst_rchl))<<" GeV\n";
  msg_Out()<<"###   a1(1260) (inline, RChL::Gamma_a1_PionFit): M = "
	   <<m_Ma1_rchl<<" GeV\n";
}

Complex F1_0_KKPi::FF_RChL2012(const double & s123,const double & s1_arg,
			       const double & s2_arg) {
  if (p_GRho==NULL || p_GKst==NULL) return Complex(0.,0.);
  double q2 = s123;
  double paper_s2 = s1_arg, paper_s1 = s2_arg; // see class-level comment
  double paper_s3 = q2 - paper_s1 - paper_s2 + 2.*m_mK2 + m_mpi2;

  Complex propRho(1./Complex(sqr(m_Mrho_rchl)-paper_s2,
			     -m_Mrho_rchl*(*p_GRho)(paper_s2)));
  Complex propKst(1./Complex(sqr(m_MKst_rchl)-paper_s1,
			     -m_MKst_rchl*(*p_GKst)(paper_s1)));
  // NB: this a1 propagator sign convention (Ma1^2-q^2, not q^2-Ma1^2)
  // is exactly as printed in Eq.(13)/(14) of 1203.3955 - see the file
  // header note on the corresponding sign in Eq.(6)/(8) for 3pi, which
  // uses the opposite convention. Transcribed literally, not "fixed".
  Complex a1prop = q2/Complex(sqr(m_Ma1_rchl)-q2,
			      -m_Ma1_rchl*RChL::Gamma_a1_PionFit(q2,m_Mrho_rchl,m_mpi2));

  if (!m_isF2) {
    Complex Fchi1(-sqrt(2.)/3.,0.);
    Complex FR1 = -(sqrt(2.)/6.)*(m_FV_rchl*m_GV_rchl)/sqr(m_F_rchl) *
      ( RChL::BR(paper_s1,paper_s3,m_mK2,m_mK2,m_GV_rchl,m_FV_rchl) * propRho +
	RChL::AR(q2,paper_s1,paper_s3,m_mK2,m_mK2,m_mpi2,m_GV_rchl,m_FV_rchl) * propKst );
    Complex FRR1 = (2./3.)*(m_FA_rchl*m_GV_rchl)/sqr(m_F_rchl) * a1prop *
      ( RChL::BRR(q2,paper_s1,paper_s3,paper_s2,m_mK2,m_mK2,m_mpi2,0.,0.) * propRho +
	RChL::ARR(q2,paper_s1,paper_s3,m_mK2,m_mK2,m_mpi2,0.,0.) * propKst );
    return Fchi1+FR1+FRR1;
  }
  else {
    Complex Fchi2(-sqrt(2.)/3.,0.); // Fchi2 = Fchi1
    Complex FR2 = -(sqrt(2.)/6.)*(m_FV_rchl*m_GV_rchl)/sqr(m_F_rchl) *
      ( RChL::AR(q2,paper_s2,paper_s3,m_mK2,m_mpi2,m_mK2,m_GV_rchl,m_FV_rchl) * propRho +
	RChL::BR(paper_s2,paper_s3,m_mK2,m_mpi2,m_GV_rchl,m_FV_rchl) * propKst );
    Complex FRR2 = (2./3.)*(m_FA_rchl*m_GV_rchl)/sqr(m_F_rchl) * a1prop *
      ( RChL::ARR(q2,paper_s2,paper_s3,m_mK2,m_mpi2,m_mK2,0.,0.) * propRho +
	RChL::BRR(q2,paper_s2,paper_s3,paper_s1,m_mK2,m_mpi2,m_mK2,0.,0.) * propKst );
    return Fchi2+FR2+FRR2;
  }
}

//////////////////////////////////////////////////////////////////////////////////
//
// K^+ K^- pi^- vector form factor F5, RChL2012 (Eq.15 of 1203.3955).
// Same index-mapping caveat as F1_0_KKPi above.
//
//////////////////////////////////////////////////////////////////////////////////

class F5_0_KKPi : public FF_0_PPP_Base {
  double  m_F_rchl, m_FV_rchl, m_GV_rchl, m_MV_rchl, m_MKst_rchl;
  double  m_mK2, m_mpi2, m_thetaV;
  double  m_c125, m_c1256, m_c1238, m_c4, m_d3, m_d123, m_g13, m_g2, m_g4, m_g5;
  Total_Width_Base * p_GRho, * p_GKst, * p_GOmega, * p_GPhi;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
  Complex FF_RChL2012(const double & s123,const double & s1,const double & s2);
public:
  F5_0_KKPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_GRho(NULL), p_GKst(NULL), p_GOmega(NULL), p_GPhi(NULL)
  {
    Construct(params);
  }
};

void F5_0_KKPi::Construct(const FF_Parameters & params) {
  m_norm    = (*params.p_model)("Vud",Tools::Vud) /
              (*params.p_model)("F_rchlKKpi",0.0924);
  m_F_rchl  = (*params.p_model)("F_rchlKKpi",0.0924);
  m_FV_rchl = (*params.p_model)("FV_rchlKKpi",0.18);
  m_GV_rchl = sqr(m_F_rchl)/m_FV_rchl;
  m_MV_rchl = (*params.p_model)("Mrho_rchlKKpi",0.775);
  m_MKst_rchl = (*params.p_model)("MKst_rchlKKpi",0.8953);
  m_mK2     = sqr(Flavour(kf_K_plus).HadMass());
  m_mpi2    = sqr(Flavour(kf_pi_plus).HadMass());
  m_thetaV  = (*params.p_model)("thetaV_rchlKKpi",35.26*M_PI/180.);
  // Table 5 defaults (channels 7-9). Only the combinations Table 5
  // actually tabulates are used - see RChL_Functions.H's note.
  m_c125  = (*params.p_model)("c125_rchlKKpi",0.);
  m_c1256 = (*params.p_model)("c1256_rchlKKpi",
			      -3.*m_FV_rchl*m_MV_rchl/(96.*sqr(M_PI)*sqrt(2.)*sqr(m_F_rchl)));
  m_c1238 = (*params.p_model)("c1238_rchlKKpi",0.);
  m_c4    = (*params.p_model)("c4_rchlKKpi",-0.07);
  m_d3    = (*params.p_model)("d3_rchlKKpi",-sqr(m_MV_rchl)/(64.*sqr(M_PI)*sqr(m_F_rchl)));
  m_d123  = (*params.p_model)("d123_rchlKKpi",0.05);
  m_g13   = (*params.p_model)("g13_rchlKKpi",
			      -2.*m_MV_rchl/(192.*sqr(M_PI)*sqrt(2.)*m_FV_rchl));
  m_g2    = (*params.p_model)("g2_rchlKKpi",
			      m_MV_rchl/(192.*sqr(M_PI)*sqrt(2.)*m_FV_rchl));
  m_g4    = (*params.p_model)("g4_rchlKKpi",-0.72);
  m_g5    = (*params.p_model)("g5_rchlKKpi",-0.6-2.*m_g4); // 2g4+g5=-0.6

  p_GRho   = LineShapes->Get(Flavour(kf_rho_770_plus));
  p_GKst   = LineShapes->Get(Flavour(kf_K_star_892_plus));
  p_GOmega = LineShapes->Get(Flavour(kf_omega_782));
  // phi(1020) lineshape now registered (Omega_Decays.H/.C,
  // kf_phi_1020=333) - this comment previously said none existed;
  // fixed to actually use it. At the default IDEAL mixing angle
  // (thetaV=35.26deg) its contribution still vanishes identically in
  // Eq.(15) (cos^2(thetaV)(1-sqrt(2)tan(thetaV))=0), so this only
  // matters if thetaV_rchlKKpi is overridden away from ideal mixing -
  // previously that case was silently missing the phi piece entirely;
  // now it is included (see FF_RChL2012 below, omegaPhiTerm).
  p_GPhi   = LineShapes->Get(Flavour(kf_phi_1020));
  if (p_GPhi==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": phi(1020) lineshape unexpectedly "
	       <<"unavailable for F5_0_KKPi - the phi piece of Eq.(15) will "
	       <<"be dropped (harmless at the default ideal mixing angle, "
	       <<"where it vanishes identically anyway, but not otherwise).\n";
  }
  msg_Out()<<"### Propagator structure for \"F5_0_KKPi, mode="<<int(m_mode)
	   <<"\" (FORM_FACTOR = "<<int(m_ffmodel)<<"):\n";
  if (p_GRho!=NULL)   msg_Out()<<"###   rho(770): M = "<<m_MV_rchl<<" GeV,  "
			       <<"Gamma(M^2) = "<<(*p_GRho)(sqr(m_MV_rchl))<<" GeV\n";
  if (p_GKst!=NULL)   msg_Out()<<"###   K*(892): M = "<<m_MKst_rchl<<" GeV,  "
			       <<"Gamma(M^2) = "<<(*p_GKst)(sqr(m_MKst_rchl))<<" GeV\n";
  if (p_GOmega!=NULL) msg_Out()<<"###   omega(782): Gamma(0.78194^2) = "
			       <<(*p_GOmega)(sqr(0.78194))<<" GeV\n";
  if (p_GPhi!=NULL)   msg_Out()<<"###   phi(1020): Gamma(1.020^2) = "
			       <<(*p_GPhi)(sqr(1.020))<<" GeV,  mixing angle = "
			       <<m_thetaV*180./M_PI<<" deg\n";
}

Complex F5_0_KKPi::FF_RChL2012(const double & s123,const double & s1_arg,
			       const double & s2_arg) {
  if (p_GRho==NULL || p_GKst==NULL) return Complex(0.,0.);
  double q2 = s123;
  double paper_s2 = s1_arg, paper_s1 = s2_arg;

  double Mrho2 = sqr(m_MV_rchl);
  Complex propRhoQ2(1./Complex(Mrho2-q2,-m_MV_rchl*(*p_GRho)(q2)));
  Complex propKst(1./Complex(sqr(m_MKst_rchl)-paper_s1,
			     -m_MKst_rchl*(*p_GKst)(paper_s1)));
  double  sin2V=sqr(sin(m_thetaV)), cos2V=sqr(cos(m_thetaV));
  double  Momega=0.78194, Gomega=0.00843;
  Complex propOmega(1./Complex(sqr(Momega)-paper_s2,
			       -Momega*(p_GOmega!=NULL?(*p_GOmega)(paper_s2):Gomega)));
  double  mixOmega = 1.+sqrt(2.)/tan(m_thetaV);
  double  mixPhi   = 1.-sqrt(2.)*tan(m_thetaV); // ->0 at ideal mixing
  // phi(1020) piece: previously always dropped (no lineshape existed).
  // Now included whenever p_GPhi is available - vanishes identically
  // at the default ideal mixing angle (mixPhi=0 there) regardless, so
  // this only changes anything if thetaV_rchlKKpi is overridden.
  double  Mphi=1.020, Gphi=0.00426;
  Complex propPhi(1./Complex(sqr(Mphi)-paper_s2,
			     -Mphi*(p_GPhi!=NULL?(*p_GPhi)(paper_s2):Gphi)));
  Complex omegaPhiTerm = sin2V*mixOmega*propOmega + cos2V*mixPhi*propPhi;

  double CR_s2 = RChL::CR(q2,paper_s2,m_mK2,m_mK2,m_mpi2,
			  m_c125,m_c1256,m_c1238,m_c4);
  double CR_s1 = RChL::CR(q2,paper_s1,m_mK2,m_mpi2,m_mK2,
			  m_c125,m_c1256,m_c1238,m_c4);
  double DR_val = RChL::DR(q2,paper_s2,paper_s1,m_mK2,m_mpi2,m_g13,m_g2,m_g4,m_g5);

  Complex Fchi5(sqrt(2.),0.);
  Complex FR5 = (16.*sqr(M_PI)*m_GV_rchl/m_MV_rchl) *
    ( CR_s2*omegaPhiTerm + CR_s1*propKst -
      (2.*m_FV_rchl/m_GV_rchl)*DR_val*propRhoQ2 );

  double CRR_s1 = RChL::CRR(q2,paper_s1,m_mK2,m_d3,m_d123);
  double CRR_s2 = RChL::CRR(q2,paper_s2,m_mpi2,m_d3,m_d123);
  // NOTE: this second occurrence of the omega piece is left as the
  // omega-only "sin2V*mixOmega*propOmega" combination, NOT the
  // omegaPhiTerm used in FR5 above - unlike that one, there is no
  // comment/derivation here confirming the same phi structure applies,
  // so this is not extended by assumption. Flag if you know this
  // should also include the phi piece.
  Complex FRR5 = -16.*sqrt(2.)*sqr(M_PI)*m_FV_rchl*m_GV_rchl * propRhoQ2 *
    ( CRR_s1*propKst + CRR_s2*sin2V*mixOmega*propOmega );

  return Fchi5+FR5+FRR5;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// "KS"-mode (Kuehn-Santamaria-style) form factors from Finkemeier &
// Mirkes, hep-ph/9503474 (FM95), for the KpiK isospin family:
// K^-pi^-K^+, K^0pi^-K0bar, K_Spi^-K_S, K_Spi^-K_L, K_Lpi^-K_L,
// K^-pi^0K^0 (Tabs.I/II, rows 1-5, plus Eq.(49)/(50)/(52) for the
// K_S/K_L final states, Sec.VI).
//
// Momentum convention: FixMode() (above) assigns p1=q1,p2=q2,p3=q3 in
// FM95's own per-channel order (see the enum comments in FF_0_PPP.H
// and the per-channel notes below), i.e. the SAME base-class
// convention as the pion modes: s12=(p1+p2)^2=paper's s3,
// s13=(p1+p3)^2=paper's s2, and paper's s1=(p2+p3)^2 is reconstructed
// via momentum conservation. This current is wired through the new
// VA_0_KPiK Current class (see VA_0_KPiK.H/.C), which uses the SAME
// (corrected) v1=(p1-p3)_T, v2=(p2-p3)_T basis as VA_0_PiPiPi.
//
// Shared resonances used throughout this section:
//  - BW_A1(Q^2): a1(1260), Eq.(30)/(31) - already available
//    (kf_a_1_1260_plus), reused as-is.
//  - T_rho^(1)(s) = [BW_rho(770) + beta_rho BW_rho(1370)]/(1+beta_rho),
//    beta_rho=-0.145, Eq.(8)/(9). NOTE: this is a SIMPLER 2-resonance
//    mix than the 3-resonance (rho+rho'+rho'', gamma=-0.167,delta=0.05)
//    KS parametrization already used for the 2pi/KK channels in
//    FF_0_PP.C - that mismatch was already flagged there and is
//    reproduced/reinforced here (FM95 independently corroborates
//    beta_rho=-0.145 with NO third resonance). Frank: please advise
//    whether FF_0_PP.C's pipi_plus/KK_plus gamma/delta should be
//    changed to match this simpler form, or whether these are
//    deliberately different fits for different purposes.
//  - T_K*^(1)(s) = [BW_K*(892) + beta_K* BW_K*(1410)]/(1+beta_K*),
//    beta_K*=-0.135, Eq.(10) - already available (used by FF_0_PP.C's
//    Kpi_plus mode), reused as-is.
//  - T_rho^(2)(s), T_K*^(2)(s): THREE-resonance vector-current mixes,
//    Eq.(42), lambda=-0.25, mu=-0.038. K*(1680)/K*'' now built and
//    registered (Kstar_Decays.H/.C, kf 30313/30323); reused here.
//  - T_omega(s) = [BW_omega(782) + eps BW_phi(1020)]/(1+eps), eps=0.05,
//    Eq.(36)/(37)/(38). kf_phi_1020=333 confirmed against the real
//    ATOOLS/Phys/Flavour_Tags.H (Omega_Decays.H/.C).
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_KPiK : public FF_0_PPP_Base {
  bool    m_isF2;
  Propagator_Base   * p_BWA1;
  Summed_Propagator * p_Trho1, * p_TKstar1;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  F1_0_KPiK(const FF_Parameters & params) :
    FF_0_PPP_Base(params), m_isF2(false),
    p_BWA1(NULL), p_Trho1(NULL), p_TKstar1(NULL)
  {
    if (params.m_name=="F2_0_KPiK") m_isF2 = true;
    Construct(params);
  }
  ~F1_0_KPiK() {
    if (p_BWA1)    { delete p_BWA1;    p_BWA1    = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
  }
};

void F1_0_KPiK::Construct(const FF_Parameters & params) {
  if (m_ffmodel!=ff_model::KS) return;
  double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vud = (*params.p_model)("Vud",Tools::Vud);
  double pref = 2.*sqrt(2.)/(3.*fpi);
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_KP:
  case FF_0_PPP_mode::K0_piM_K0bar:
    m_norm = Complex(-pref*Vud/2., 0.); break;                 // Tab.I rows 1,2
  case FF_0_PPP_mode::KS_piM_KS:
    m_norm = Complex(-pref*Vud/4., 0.); break;                 // row 3
  case FF_0_PPP_mode::KS_piM_KL:
    m_norm = Complex(-pref*Vud/4., 0.); break;                 // row 4
  case FF_0_PPP_mode::KL_piM_KL:
    m_norm = Complex( pref*Vud/4., 0.); break;                 // Eq.(50): = -KSpiKS
  case FF_0_PPP_mode::KM_pi0_K0:
    m_norm = Complex( pref*3.*Vud/(2.*sqrt(2.)), 0.);          // row 5
    // K_S/K_L projection - same reasoning as FF_0_PP.C's Kpi_plus/
    // KK_plus and this file's piM_K0bar_pi0 (m_isKSKL in FF_0_PPP.H).
    if (m_isKSKL) m_norm /= sqrt(2.);
    break;
  default: m_norm = Complex(0.,0.); break;
  }

  p_BWA1 = new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
  // T_rho^(1)/T_K*^(1): reuse the SAME unified KS parameter names as
  // the two-meson pipi_plus/Kpi_plus channels (FF_0_PP.C) - these are
  // the identical physical rho(770)+beta_rho*rho(1450) and
  // K*(892)+beta_K* K*(1410) mixes, just appearing again inside this
  // 3-body current; tuning the 2-body fit now automatically keeps
  // this 3-body current consistent with it, rather than needing two
  // independently-maintained copies of the same numbers.
  Complex gammaRho  = ReadComplexParam(params.p_model,
				       "gammaMag_pipi3",-0.145,"gammaPhase_pipi3");
  Complex gammaKst  = ReadComplexParam(params.p_model,
				       "gammaMag_Kpi_100",-0.135,"gammaPhase_Kpi_100");
  Propagator_Base * rho770  =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
  Propagator_Base * rho1450 =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, gammaRho); // beta_rho, Eq.(9)
  Propagator_Base * Kstar892  =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, gammaKst); // beta_K*, Eq.(10)

  std::string label = std::string("F1_0_KPiK, mode=")+std::to_string(int(m_mode))+
                       (m_isF2 ? " (F2)" : " (F1)");
  DumpPropagatorStructure(label+" [a1]",       int(m_ffmodel), p_BWA1);
  DumpPropagatorStructure(label+" [T_rho^1]",  int(m_ffmodel), p_Trho1);
  DumpPropagatorStructure(label+" [T_K*^1]",   int(m_ffmodel), p_TKstar1);
}

Complex F1_0_KPiK::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_BWA1==NULL || p_Trho1==NULL || p_TKstar1==NULL) return Complex(0.,0.);
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  Complex BWA1 = (*p_BWA1)(s123);
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_KP:
  case FF_0_PPP_mode::K0_piM_K0bar:
    // G1=BWA1*T_rho1(s2=s13); G2=BWA1*T_K*1(s1)
    return (m_isF2 ? BWA1*(*p_TKstar1)(s1) : BWA1*(*p_Trho1)(s13));
  case FF_0_PPP_mode::KS_piM_KS:
  case FF_0_PPP_mode::KL_piM_KL:
    // G1=BWA1*T_K*1(s3=s12); G2=-BWA1*[T_K*1(s1)+T_K*1(s3=s12)]
    return ( m_isF2 ?
	     -BWA1*((*p_TKstar1)(s1)+(*p_TKstar1)(s12)) :
	     BWA1*(*p_TKstar1)(s12) );
  case FF_0_PPP_mode::KS_piM_KL:
    // G1=BWA1*[2T_rho1(s2=s13)+T_K*1(s3=s12)]; G2=BWA1*[T_K*1(s1)-T_K*1(s3=s12)]
    return ( m_isF2 ?
	     BWA1*((*p_TKstar1)(s1)-(*p_TKstar1)(s12)) :
	     BWA1*(2.*(*p_Trho1)(s13)+(*p_TKstar1)(s12)) );
  case FF_0_PPP_mode::KM_pi0_K0:
    // G1=BWA1*[2/3 T_rho1(s2=s13)+1/3 T_K*1(s3=s12)];
    // G2=1/3 BWA1*[T_K*1(s1)-T_K*1(s3=s12)]
    return ( m_isF2 ?
	     (1./3.)*BWA1*((*p_TKstar1)(s1)-(*p_TKstar1)(s12)) :
	     BWA1*((2./3.)*(*p_Trho1)(s13)+(1./3.)*(*p_TKstar1)(s12)) );
  default: return Complex(0.,0.);
  }
}

//////////////////////////////////////////////////////////////////////////////////
//
// Vector ("anomaly") form factor F3(=G3, code name FS) for the KpiK
// family, Tab.II rows 1-5 + Eq.(52) for K_Lpi^-K_L.
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_KPiK : public FF_0_PPP_Base {
  Summed_Propagator * p_Trho2, * p_TKstar1, * p_Tomega;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  FS_0_KPiK(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_Trho2(NULL), p_TKstar1(NULL), p_Tomega(NULL)
  {
    Construct(params);
  }
  ~FS_0_KPiK() {
    if (p_Trho2)   { delete p_Trho2;   p_Trho2   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
    if (p_Tomega)  { delete p_Tomega;  p_Tomega  = NULL; }
  }
};

void FS_0_KPiK::Construct(const FF_Parameters & params) {
  if (m_ffmodel!=ff_model::KS) return;
  double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vud = (*params.p_model)("Vud",Tools::Vud);
  double pref = 1./(2.*sqrt(2.)*sqr(M_PI)*pow(fpi,3));
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_KP:    m_norm = Complex(-pref*Vud,0.);        break;
  case FF_0_PPP_mode::K0_piM_K0bar: m_norm = Complex( pref*Vud,0.);        break;
  case FF_0_PPP_mode::KS_piM_KS:    m_norm = Complex(-pref*Vud/2.,0.);     break;
  case FF_0_PPP_mode::KS_piM_KL:    m_norm = Complex( pref*Vud/2.,0.);     break;
  case FF_0_PPP_mode::KL_piM_KL:    m_norm = Complex( pref*Vud/2.,0.);     break; // Eq.(52): =-KSpiKS
  case FF_0_PPP_mode::KM_pi0_K0:
    m_norm = Complex(-pref*Vud/sqrt(2.),0.);
    if (m_isKSKL) m_norm /= sqrt(2.); // same K_S/K_L projection as F1_0_KPiK
    break;
  default: m_norm = Complex(0.,0.); break;
  }

  // T_rho^(2), Eq.(42): rho(770)+lambda rho(1500)+mu rho(1750). Reuses
  // kf_rho_1450_plus/kf_rho_1700_plus as the closest available
  // resonances - FLAG: FM95 wants DIFFERENT masses/widths for these
  // (1.500/0.220 and 1.750/0.120 GeV) than the axial T_rho^(1) above
  // (1.370/0.510) or the original Kuhn-Santamaria pi-pi fit already in
  // FF_0_PP.C. Since BreitWigner's pole mass is tied to the Flavour's
  // own HadMass(), the SAME kf-coded object cannot simultaneously carry
  // three different pole masses for three different papers. Please
  // advise: either accept whichever mass is in Sherpa's particle table
  // for kf_rho_1450_plus/kf_rho_1700_plus everywhere, or provide new,
  // distinguishable kf-codes/parametrised Flavours for this specific
  // 3-body vector-current usage.
  Complex lambda = ReadComplexParam(params.p_model,
				    "gammaMag_KpiK_vector",-0.25,"gammaPhase_KpiK_vector");
  Complex mu     = ReadComplexParam(params.p_model,
				    "deltaMag_KpiK_vector",-0.038,"deltaPhase_KpiK_vector");
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
  Propagator_Base * rho1500 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho2 = new Summed_Propagator();
  p_Trho2->Add(rho770,  Complex(1.,0.));
  p_Trho2->Add(rho1500, lambda);
  Total_Width_Base * wRhopp = LineShapes->Get(Flavour(kf_rho_1700_plus));
  if (wRhopp!=NULL) p_Trho2->Add(new BreitWigner(wRhopp), mu);

  // T_K*^(1): same unified KS parameter name as the two-meson Kpi_plus
  // channel (see F1_0_KPiK's Construct() comment above for rationale).
  Complex gammaKst = ReadComplexParam(params.p_model,
				      "gammaMag_Kpi_100",-0.135,"gammaPhase_Kpi_100");
  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, gammaKst);

  // T_omega, Eq.(36): omega(782)+eps phi(1020). phi(1020) confirmed
  // and registered (Omega_Decays.H/.C, kf_phi_1020=333).
  Complex eps = ReadComplexParam(params.p_model,
				 "gammaMag_omegaphi",0.05,"gammaPhase_omegaphi");
  Propagator_Base * omega782 = new BreitWigner(LineShapes->Get(Flavour(kf_omega_782)));
  p_Tomega = new Summed_Propagator();
  p_Tomega->Add(omega782, Complex(1.,0.));
  Total_Width_Base * wPhi = LineShapes->Get(Flavour(kf_phi_1020));
  if (wPhi!=NULL) p_Tomega->Add(new BreitWigner(wPhi), eps);
  else msg_Error()<<"Error in "<<METHOD<<": phi(1020) lineshape "
		  <<"unexpectedly unavailable for T_omega (FM95 Eq.36) - "
		  <<"falling back to plain BW_omega(782), a ~5% approximation.\n";

  std::string label = std::string("FS_0_KPiK, mode=")+std::to_string(int(m_mode));
  DumpPropagatorStructure(label+" [T_rho^2]",  int(m_ffmodel), p_Trho2);
  DumpPropagatorStructure(label+" [T_K*^1]",   int(m_ffmodel), p_TKstar1);
  DumpPropagatorStructure(label+" [T_omega]",  int(m_ffmodel), p_Tomega);
}

Complex FS_0_KPiK::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_Trho2==NULL || p_TKstar1==NULL || p_Tomega==NULL) return Complex(0.,0.);
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  Complex Trho2 = (*p_Trho2)(s123);
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_KP:
  case FF_0_PPP_mode::K0_piM_K0bar:
    // (sqrt2-1)[sqrt2 T_omega(s2=s13) + T_K*1(s1)]
    return (sqrt(2.)-1.)*Trho2*( sqrt(2.)*(*p_Tomega)(s13) + (*p_TKstar1)(s1) );
  case FF_0_PPP_mode::KS_piM_KS:
  case FF_0_PPP_mode::KL_piM_KL:
    // (sqrt2-1)[T_K*1(s1) - T_K*1(s3=s12)]
    return (sqrt(2.)-1.)*Trho2*( (*p_TKstar1)(s1) - (*p_TKstar1)(s12) );
  case FF_0_PPP_mode::KS_piM_KL:
    // (sqrt2-1)[2sqrt2 T_omega(s2=s13) + T_K*1(s1) + T_K*1(s3=s12)]
    return (sqrt(2.)-1.)*Trho2*
      ( 2.*sqrt(2.)*(*p_Tomega)(s13) + (*p_TKstar1)(s1) + (*p_TKstar1)(s12) );
  case FF_0_PPP_mode::KM_pi0_K0:
    // (sqrt2-1)[T_K*1(s3=s12) - T_K*1(s1)]
    return (sqrt(2.)-1.)*Trho2*( (*p_TKstar1)(s12) - (*p_TKstar1)(s1) );
  default: return Complex(0.,0.);
  }
}

// F4 stub (always 0 - paper sets the scalar form factor to 0 for all
// channels, Sec.III below Eq.19).
class F3_0_KPiK_Stub : public FF_0_PPP_Base {
  void    Construct(const FF_Parameters & params) {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F3_0_KPiK_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// KpiPi isospin family: K^-pi^-pi^+, pi^-K0bar pi^0 (Tabs.I/II rows
// 7,8). Same momentum convention as the KpiK family above
// (p1=q1,p2=q2,p3=q3 in FM95's own order); wired through the new
// VA_0_KPiPi Current class.
//
// T_K1^(a) (K1(1400)+xi K1(1270)) and T_K1^(b) (pure K1(1270)) both
// use the K1(1270)/K1(1400) lineshapes built and registered in
// K1_Decays.H/.C (kf 10313/10323, 20313/20323, confirmed against the
// real ATOOLS/Phys/Flavour_Tags.H).
//
//////////////////////////////////////////////////////////////////////////////////

class F1_0_KPiPi : public FF_0_PPP_Base {
  bool    m_isF2;
  Complex m_xiK1;
  Summed_Propagator * p_TK1a, * p_Trho1, * p_TKstar1;
  Propagator_Base    * p_TK1b;

  // --- CLEO K1 data-driven alternative (KS_CLEO=102) ---
  // tau_two_meson_currents_KS_RChiT.tex Sec."CLEO K1 data-driven
  // alternative" - ONLY given there for the charged K^-pi^-pi^+ mode
  // (KM_piM_piP); no analogous combination is given for piM_K0bar_pi0,
  // so KS_CLEO falls back to the constant form factor for that mode.
  // Uses the SAME s_Kpi=s13/s_pipi=reconstructed-s1 assignment as the
  // FM95 KS branch (F1 uses the K*pi subchannel, F2 the Krho one), but
  // with CLEO's OWN K1(1270)/K1(1400) masses/widths (genuinely
  // broader than FM95's) and real weights A,B,C,D rather than FM95's
  // xi-mixed T_K1^(a)/T_K1^(b). The K* propagator here is a BARE
  // K*(892) (no K*(1410) admixture, unlike FM95's own T_K*^(1)) -
  // reuses Sherpa's registered running width.
  double  m_ACLEO, m_BCLEO, m_CCLEO, m_DCLEO;
  double  m_MK1_1270_CLEO, m_GK1_1270_CLEO, m_MK1_1400_CLEO, m_GK1_1400_CLEO;
  Propagator_Base * p_KstCLEO;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  F1_0_KPiPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), m_isF2(false), m_xiK1(0.33,0.),
    p_TK1a(NULL), p_Trho1(NULL), p_TKstar1(NULL), p_TK1b(NULL),
    p_KstCLEO(NULL)
  {
    if (params.m_name=="F2_0_KPiPi") m_isF2 = true;
    Construct(params);
  }
  ~F1_0_KPiPi() {
    if (p_TK1a)    { delete p_TK1a;    p_TK1a    = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
    if (p_TK1b)    { delete p_TK1b;    p_TK1b    = NULL; }
    if (p_KstCLEO) { delete p_KstCLEO; p_KstCLEO = NULL; }
  }
};

void F1_0_KPiPi::Construct(const FF_Parameters & params) {
  if (m_ffmodel==ff_model::KS || m_ffmodel==ff_model::KS_f0 ||
      m_ffmodel==ff_model::KS_flatte) {
  double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vus = (*params.p_model)("Vus",Tools::Vus);
  double pref = 2.*sqrt(2.)/(3.*fpi);
  // xi_K1 weights K1(1270) against K1(1400) in T_K1^(a), Eq.(9):
  //   T_K1^(a) = BW_K1(1400) + xi_K1 BW_K1(1270).
  // FM95 quote a real 0.33, but that was fitted alongside THEIR K1
  // line shapes.  Exposed as a complex here for two reasons: the two
  // physical K1 states are K1A/K1B (3P1/1P1) mixtures, so a relative
  // PHASE between them is physically expected and a real coefficient
  // cannot produce it; and the magnitude needs re-deriving whenever the
  // K1 line shapes change - as they did when K1(1270) was given its
  // PDG channels (see K1_Decays.C).
  //
  // The default is MODEL-DEPENDENT.  FM95's 0.33 is not a fit: Eq.(33)
  // DERIVES |xi| = 0.33 from Gamma(K1(1270)->K*pi)/Gamma(K1(1400)->K*pi) after
  // phase-space correction, against a fixed-width K1(1270).  Model 106 replaces
  // that K1(1270) with a Flatte, so the ratio the derivation rests on no longer
  // holds and xi has to be re-determined; 0.70 is the measured optimum against
  // BaBar 2007 / Belle 2010 / CLEO 2000.  Carrying FM95's 0.33 into 106 would
  // silently give a much worse model (chi2/ndf 40.9 vs 10.6 at 150k).
  m_xiK1 = ReadComplexParam(params.p_model,"xiK1Mag",
                            m_ffmodel==ff_model::KS_flatte ? 0.70 : 0.33,
                            "xiK1Phase");
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_piP:
    m_norm = Complex(-pref*Vus/2., 0.); break;                   // row 7
  case FF_0_PPP_mode::piM_K0bar_pi0:
    m_norm = Complex( pref*3.*Vus/(2.*sqrt(2.)), 0.);
    // K_S/K_L projection - see the m_isKSKL comment in FF_0_PPP.H and
    // the identical reasoning in FF_0_PP.C's Construct_KK/Construct_Kpi.
    if (m_isKSKL) m_norm /= sqrt(2.);
    break;                                                       // row 8
  default: m_norm = Complex(0.,0.); break;
  }

  Total_Width_Base * wK11400 = LineShapes->Get(Flavour(kf_K_1_1400_plus));
  Total_Width_Base * wK11270 = LineShapes->Get(Flavour(kf_K_1_1270_plus));
  if (wK11400==NULL || wK11270==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": missing K1(1270)/K1(1400) "
	       <<"lineshape(s) - T_K1^(a) and T_K1^(b) will be treated as "
	       <<"identically zero.\n";
  }
  else {
    // 106 swaps the K1(1270) Breit-Wigner for the two-channel Flatte in
    // BOTH T_K1^(a) and T_K1^(b), so each side sees the same threshold
    // dynamics.  K1(1400) stays a Breit-Wigner either way.
    const bool flatte = (m_ffmodel==ff_model::KS_flatte);
    // FM95 Eq.(34) is explicit that the K1 envelopes are CONSTANT-width,
    // normalized Breit-Wigners, with Eq.(35) mK1(1400)=1.402, G=0.174.  Model
    // 100 keeps the running-width registry propagator it was validated with;
    // 106 follows the paper, because a running Gamma(s) grows with s and
    // suppresses exactly the 1.45-1.7 GeV tail this resonance has to supply.
    // The values are overridable: they are parameters of FM95's fit, not the
    // flavour's PDG pole, so they are NOT read from the shared registry.
    // Named for the MODEL, not for a source: the same slot legitimately
    // carries either fit, and naming it after one while defaulting it to the
    // other is how the ff-ppp-work branch's "FM95_m_K1_1400 = 1.463" came to
    // be read as an FM95 number for years.  Both fits, for the record:
    //   FM95 Eq.(35), hep-ph/9503474   1.402 / 0.174
    //   CLEO (Asner et al.)            1.463 / 0.300   <- default here
    // FF_0_PPP.C's own CLEO branch carries the identical CLEO pair as
    // MK1_1400_CLEO / GK1_1400_CLEO.  106 defaults to CLEO's because the data
    // prefer it decisively (chi2/ndf 20.8 -> 10.6 at 150k, BaBar d03 12.7 ->
    // 4.0); set these to 1.402/0.174 to recover FM95's own K1(1400).
    const double mK1400 = (*params.p_model)("MK1_1400_106", 1.463);
    const double GK1400 = (*params.p_model)("GK1_1400_106", 0.300);
    p_TK1a = new Summed_Propagator();
    p_TK1a->Add(flatte ? (Propagator_Base*)new FM95_Fixed_BW(mK1400,GK1400)
                       : (Propagator_Base*)new BreitWigner(wK11400),
                Complex(1.,0.));
    p_TK1a->Add(flatte ? MakeK1_1270_Flatte()
                       : (Propagator_Base*)new BreitWigner(wK11270), m_xiK1);
    p_TK1b = flatte ? MakeK1_1270_Flatte()
                    : (Propagator_Base*)new BreitWigner(wK11270);
  }
  // Same unified KS parameter names as pipi_plus/Kpi_plus (FF_0_PP.C)
  // and F1_0_KPiK above - same physical rho(770)+rho(1450) and
  // K*(892)+K*(1410) mixes reused yet again.
  Complex gammaRho = ReadComplexParam(params.p_model,
				      "gammaMag_pipi3",-0.145,"gammaPhase_pipi3");
  Complex gammaKst = ReadComplexParam(params.p_model,
				      "gammaMag_Kpi_100",-0.135,"gammaPhase_Kpi_100");
  // 106 puts rho(770) on a Gounaris-Sakurai shape.  GS carries the dispersive
  // real part of the self-energy that a running-width BW drops, which is the
  // correct analytic description of a broad resonance sitting on its own
  // threshold (Gounaris & Sakurai, PRL 21 (1968) 244).  Measured on the ff-ppp
  // branch as -29%/-45% on the 3pi spectra and a wash on Kpipi; it is taken
  // here on correctness grounds, not because it moves this channel.
  // rho(1450) deliberately stays a Breit-Wigner: GS assumes an ELASTIC pi pi
  // width, true of rho(770) at 99.9% but not of rho(1450) (19.7% pi pi).
  // NOTE this is a departure from FM95, which uses a running-width BW
  // throughout - same category as the Flatte K1(1270).
  const resonance_type rho_type = (m_ffmodel==ff_model::KS_flatte) ?
                                  resonance_type::GS : resonance_type::running;
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)),
                                              rho_type);
  Propagator_Base * rho1450 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, gammaRho);
  // KS_f0 (103): add an f0(500)/sigma admixture to the pipi tower -
  // NOT present in FM95's original current, requested addition. No
  // literature source for this specific weight - default chosen to
  // be a modest, same-order-of-magnitude admixture to rho(1450)'s own
  // weight, reusing Sherpa's own registered f0(600) running width
  // (kf_f_0_600) per the established "reuse lineshape machinery"
  // preference. Please retune deltaMag_pipi_f0 against data.
  if (m_ffmodel==ff_model::KS_f0) {
    Complex deltaF0 = ReadComplexParam(params.p_model,
				       "deltaMag_pipi_f0",-0.2,"deltaPhase_pipi_f0");
    p_Trho1->Add(new BreitWigner(LineShapes->Get(Flavour(kf_f_0_600))), deltaF0);
  }
  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, gammaKst);

  std::string label = std::string("F1_0_KPiPi, mode=")+std::to_string(int(m_mode))+
                       (m_isF2 ? " (F2)" : " (F1)");
  DumpPropagatorStructure(label+" [T_K1^a]",  int(m_ffmodel), p_TK1a);
  DumpPropagatorStructure(label+" [T_K1^b]",  int(m_ffmodel), p_TK1b);
  DumpPropagatorStructure(label+" [T_rho^1]", int(m_ffmodel), p_Trho1);
  DumpPropagatorStructure(label+" [T_K*^1]",  int(m_ffmodel), p_TKstar1);
  }
  else if ((m_ffmodel==ff_model::KS_CLEO || m_ffmodel==ff_model::CLEO_f0) &&
	   m_mode==FF_0_PPP_mode::KM_piM_piP) {
    double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
    double Vus = (*params.p_model)("Vus",Tools::Vus);
    // FIXME: the note gives the "-2N/(3Fpi)"/"-N/(sqrt3 Fpi)" prefactor
    // STRUCTURE for F1/F2 respectively, but no explicit numeric value
    // for the overall N (unlike FM95's own A^{K-pi-pi+}=-sin(thetaC)/2,
    // which gives a concrete CKM coefficient). N_Kpipi_CLEO below is a
    // placeholder (default 1) for whatever residual normalization
    // convention N represents - please supply if known.
    double NCLEO = (*params.p_model)("N_Kpipi_CLEO", 1.);
    if (!m_isF2) m_norm = Complex(-2.*NCLEO*Vus/(3.*fpi), 0.);
    else         m_norm = Complex(   -NCLEO*Vus/(sqrt(3.)*fpi), 0.);

    m_CCLEO = (*params.p_model)("CCLEO_Kpipi", 0.20);
    m_DCLEO = (*params.p_model)("DCLEO_Kpipi", 0.27);
    m_ACLEO = (*params.p_model)("ACLEO_Kpipi", 0.94);
    m_BCLEO = (*params.p_model)("BCLEO_Kpipi", 0.);
    // CLEO-fitted K1(1270)/K1(1400) - genuinely broader than FM95's
    // own values (0.090/0.174 GeV) - own masses/widths, not reused
    // from K1_Decays.H's registered lineshapes (built for FM95's own,
    // much narrower, fit). Running-width convention for these two is
    // NOT specified in the note for this channel - a plain FIXED
    // width is used here (most conservative reading), unlike the
    // explicit BW_L prescription given for the 3pi CLEO current;
    // flag if a running width is wanted instead.
    m_MK1_1270_CLEO = (*params.p_model)("MK1_1270_CLEO", 1.254);
    m_GK1_1270_CLEO = (*params.p_model)("GK1_1270_CLEO", 0.26);
    m_MK1_1400_CLEO = (*params.p_model)("MK1_1400_CLEO", 1.463);
    m_GK1_1400_CLEO = (*params.p_model)("GK1_1400_CLEO", 0.30);

    Complex gammaRho = ReadComplexParam(params.p_model,
					"gammaMag_pipi3",-0.145,"gammaPhase_pipi3");
    Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
    Propagator_Base * rho1450 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    p_Trho1 = new Summed_Propagator();
    p_Trho1->Add(rho770,  Complex(1.,0.));
    p_Trho1->Add(rho1450, gammaRho);
    // CLEO_f0 (104): same f0(500)/sigma addition as KS_f0 above -
    // see that branch's comment for the caveat on this weight.
    if (m_ffmodel==ff_model::CLEO_f0) {
      Complex deltaF0 = ReadComplexParam(params.p_model,
					 "deltaMag_pipi_f0",-0.2,"deltaPhase_pipi_f0");
      p_Trho1->Add(new BreitWigner(LineShapes->Get(Flavour(kf_f_0_600))), deltaF0);
    }
    p_KstCLEO = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));

    std::string label = std::string("F1_0_KPiPi (CLEO), mode=")+std::to_string(int(m_mode))+
                         (m_isF2 ? " (F2)" : " (F1)");
    DumpPropagatorStructure(label+" [T_rho^1]", int(m_ffmodel), p_Trho1);
    DumpPropagatorStructure(label+" [K*(892)]", int(m_ffmodel), p_KstCLEO);
    msg_Out()<<"###   K1(1270) (CLEO): M = "<<m_MK1_1270_CLEO<<" GeV, Gamma = "
	     <<m_GK1_1270_CLEO<<" GeV (fixed width)\n"
	     <<"###   K1(1400) (CLEO): M = "<<m_MK1_1400_CLEO<<" GeV, Gamma = "
	     <<m_GK1_1400_CLEO<<" GeV (fixed width)\n";
  }
}

Complex F1_0_KPiPi::
FF_KS(const double & s123,const double & s12,const double & s13) {
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  if ((m_ffmodel==ff_model::KS_CLEO || m_ffmodel==ff_model::CLEO_f0) &&
      m_mode==FF_0_PPP_mode::KM_piM_piP) {
    if (p_Trho1==NULL || p_KstCLEO==NULL) return Complex(0.,0.);
    // s_Kpi=s13, s_pipi=s1 (reconstructed) - same assignment as the
    // FM95 KS branch below (F1 uses the K*pi subchannel, F2 the Krho
    // one). Fixed-width K1(1270)/K1(1400) BWs (see Construct()).
    Complex BWK1_1270 = sqr(m_MK1_1270_CLEO) /
      Complex(sqr(m_MK1_1270_CLEO)-s123,-m_MK1_1270_CLEO*m_GK1_1270_CLEO);
    Complex BWK1_1400 = sqr(m_MK1_1400_CLEO) /
      Complex(sqr(m_MK1_1400_CLEO)-s123,-m_MK1_1400_CLEO*m_GK1_1400_CLEO);
    if (!m_isF2) return (m_CCLEO*BWK1_1270 + m_DCLEO*BWK1_1400) * (*p_KstCLEO)(s13);
    return (m_ACLEO*BWK1_1270 + m_BCLEO*BWK1_1400) * (*p_Trho1)(s1);
  }
  if (p_Trho1==NULL || p_TKstar1==NULL) return Complex(0.,0.);
  Complex TK1a = (p_TK1a!=NULL ? (*p_TK1a)(s123) : Complex(0.,0.));
  Complex TK1b = (p_TK1b!=NULL ? (*p_TK1b)(s123) : Complex(0.,0.));
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_piP:
    // G1=T_K1a*T_K*1(s2=s13); G2=T_K1b*T_rho1(s1)
    return (m_isF2 ? TK1b*(*p_Trho1)(s1) : TK1a*(*p_TKstar1)(s13));
  case FF_0_PPP_mode::piM_K0bar_pi0:
    // G1=2/3 T_K1b*T_rho1(s2=s13) + 1/3 T_K1a*T_K*1(s3=s12)
    // G2=1/3 T_K1a*[T_K*1(s1) - T_K*1(s3=s12)]
    return ( m_isF2 ?
	     (1./3.)*TK1a*((*p_TKstar1)(s1)-(*p_TKstar1)(s12)) :
	     (2./3.)*TK1b*(*p_Trho1)(s13) + (1./3.)*TK1a*(*p_TKstar1)(s12) );
  default: return Complex(0.,0.);
  }
}

//////////////////////////////////////////////////////////////////////////////////
//
// Vector form factor for the KpiPi family, Tab.II rows 7,8. Uses
// T_K*^(2) - K*(1680)/K*'' built and registered, same as
// FS_0_PiPlusPiZeroPiZero above.
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_KPiPi : public FF_0_PPP_Base {
  Summed_Propagator * p_TKstar2, * p_Trho1, * p_TKstar1;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  FS_0_KPiPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_TKstar2(NULL), p_Trho1(NULL), p_TKstar1(NULL)
  {
    Construct(params);
  }
  ~FS_0_KPiPi() {
    if (p_TKstar2) { delete p_TKstar2; p_TKstar2 = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
  }
};

void FS_0_KPiPi::Construct(const FF_Parameters & params) {
  // CLEO (102/104) set F3 = 0 in their own analysis and folded the
  // uncertainty into their systematics, so the anomalous piece is built only
  // for the FM95-family models.  Kiers et al., arXiv:0808.1707, restore F3 on
  // top of the CLEO F1/F2; that was tried here and is recorded in FF_KS below.
  if (m_ffmodel!=ff_model::KS && m_ffmodel!=ff_model::KS_f0 &&
      m_ffmodel!=ff_model::KS_flatte) return;
  double fpi  = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vus  = (*params.p_model)("Vus",Tools::Vus);
  double pref = 1./(2.*sqrt(2.)*sqr(M_PI)*pow(fpi,3));
  // NOTE on the Clebsch: these are FM95 Tab.II (vector) values, which are
  // deliberately NOT the Tab.I (axial) ones used by F1_0_KPiPi::Construct.
  // Eq.(eq:FM-3P-master) writes a single A^(abc) for F1,F2,F3, but that is
  // a simplification of the source: the note's own KKpi section gives
  // A_KKpi = -cos(theta_C)/2 for the axial and A_{V,KKpi} = -cos(theta_C)
  // for the vector, a factor 2 apart. So axial and anomalous Clebsches
  // must be read from their respective tables, not derived from each
  // other - do not ""unify"" them. (Checked 2026-08-24 while chasing the
  // Kpipi rate excess; deriving these from the axial values made the
  // agreement worse and contradicts the KKpi counter-example.)
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_piP:    m_norm = Complex(pref*Vus,0.);          break;
  case FF_0_PPP_mode::piM_K0bar_pi0:
    m_norm = Complex(pref*sqrt(2.)*Vus,0.);
    // K_S/K_L projection - see F1_0_KPiPi::Construct above.
    if (m_isKSKL) m_norm /= sqrt(2.);
    break;
  default: m_norm = Complex(0.,0.); break;
  }

  Complex lambda = ReadComplexParam(params.p_model,
				    "gammaMag_KpiK_vector",-0.25,"gammaPhase_KpiK_vector");
  Complex mu     = ReadComplexParam(params.p_model,
				    "deltaMag_KpiK_vector",-0.038,"deltaPhase_KpiK_vector");
  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar2 = new Summed_Propagator();
  p_TKstar2->Add(Kstar892,  Complex(1.,0.));
  p_TKstar2->Add(Kstar1410, lambda);
  Total_Width_Base * wKstarpp = LineShapes->Get(Flavour(kf_K_star_1680_plus));
  if (wKstarpp!=NULL) p_TKstar2->Add(new BreitWigner(wKstarpp), mu);
  else msg_Error()<<"Error in "<<METHOD<<": missing K*(1714)/K*'' lineshape "
		  <<"(FM95 Eq.42) - dropping the mu-term, see the identical "
		  <<"flag in FS_0_PiPlusPiZeroPiZero above.\n";

  Complex gammaRho = ReadComplexParam(params.p_model,
				      "gammaMag_pipi3",-0.145,"gammaPhase_pipi3");
  // 106 puts rho(770) on a Gounaris-Sakurai shape.  GS carries the dispersive
  // real part of the self-energy that a running-width BW drops, which is the
  // correct analytic description of a broad resonance sitting on its own
  // threshold (Gounaris & Sakurai, PRL 21 (1968) 244).  Measured on the ff-ppp
  // branch as -29%/-45% on the 3pi spectra and a wash on Kpipi; it is taken
  // here on correctness grounds, not because it moves this channel.
  // rho(1450) deliberately stays a Breit-Wigner: GS assumes an ELASTIC pi pi
  // width, true of rho(770) at 99.9% but not of rho(1450) (19.7% pi pi).
  // NOTE this is a departure from FM95, which uses a running-width BW
  // throughout - same category as the Flatte K1(1270).
  const resonance_type rho_type = (m_ffmodel==ff_model::KS_flatte) ?
                                  resonance_type::GS : resonance_type::running;
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),
                                              rho_type);
  Propagator_Base * rho1450 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, gammaRho);
  // KS_f0 (103): same f0(500)/sigma addition as F1_0_KPiPi - see that
  // class's Construct() for the caveat on this weight.
  if (m_ffmodel==ff_model::KS_f0) {
    Complex deltaF0 = ReadComplexParam(params.p_model,
				       "deltaMag_pipi_f0",-0.2,"deltaPhase_pipi_f0");
    p_Trho1->Add(new BreitWigner(LineShapes->Get(Flavour(kf_f_0_600))), deltaF0);
  }

  Complex gammaKst = ReadComplexParam(params.p_model,
				      "gammaMag_Kpi_100",-0.135,"gammaPhase_Kpi_100");
  Propagator_Base * Kstar892_1  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410_1 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892_1,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410_1, gammaKst);

  std::string label = std::string("FS_0_KPiPi, mode=")+std::to_string(int(m_mode));
  DumpPropagatorStructure(label+" [T_K*^2]",  int(m_ffmodel), p_TKstar2);
  DumpPropagatorStructure(label+" [T_rho^1]", int(m_ffmodel), p_Trho1);
  DumpPropagatorStructure(label+" [T_K*^1]",  int(m_ffmodel), p_TKstar1);
}

Complex FS_0_KPiPi::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_TKstar2==NULL || p_Trho1==NULL || p_TKstar1==NULL) return Complex(0.,0.);
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  Complex TKstar2 = (*p_TKstar2)(s123);
  switch (m_mode) {
  // NOTE, from chasing FM95 Tab.II back to its own primary source:
  // Decker, Mirkes, Sauer and Was, Z.Phys. C58 (1993) 445, Eq.(29) write this
  // bracket as [T_rho^(1)(s1) + alpha T_K*(s2)]/(1+alpha) and FIT alpha =
  // -0.20..-0.25 from tau -> K- pi- K+, because the K* is not seen in
  // e+e- -> K Kbar pi.  FM95's 1/2 [T_rho + T_K*] below is that expression at
  // alpha = 1, i.e. FM95 silently dropped the fitted suppression.  Decker's
  // alpha was TRIED here (2026-08-28) and made the Kpipi shape WORSE - chi2/ndf
  // 57.2 -> 63.4, degrading BABAR d04 (the K- pi+ mass, where the K* term
  // lives) 6.0 -> 8.1.  Not adopted.  Note this is a weak test of alpha: F3
  // carries ~13% of the width here, whereas Decker fitted it on K- pi- K+
  // where the anomalous piece is proportionally much larger.  If the question
  // is reopened, reopen it on KKpi.
  case FF_0_PPP_mode::KM_piM_piP:
    // G3 = 1/2 T_K*2 [T_rho1(s1) + T_K*1(s2=s13)]
    return 0.5*TKstar2*( (*p_Trho1)(s1) + (*p_TKstar1)(s13) );
  case FF_0_PPP_mode::piM_K0bar_pi0:
    // G3 = 1/4 T_K*2 [2 T_rho1(s2=s13) + T_K*1(s1) + T_K*1(s3=s12)]
    return 0.25*TKstar2*( 2.*(*p_Trho1)(s13) + (*p_TKstar1)(s1) + (*p_TKstar1)(s12) );
  default: return Complex(0.,0.);
  }
}

class F3_0_KPiPi_Stub : public FF_0_PPP_Base {
  void    Construct(const FF_Parameters & params) {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F3_0_KPiPi_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// eta pi^- pi^0 / eta' pi^- pi^0 (tau_two_meson_currents_KS_RChiT.tex,
// Sec."eta pi pi and eta' pi pi"). G-parity kills the axial (F1,F2)
// AND scalar (F4) form factors identically in the isospin limit - the
// ONLY nonzero piece is the anomalous vector current F5(=FS here),
// Eq.(eq:eta-pipi-old) (old TAUOLA/VMD current):
//   F5 = N_WZW^eta(') * T_rho^(2)(Q^2) * T_rho^(1)(s1),
//   s1 = (p2+p3)^2 = s_pipi (paper's own s1, reconstructed via
//        momentum conservation, same convention as every other
//        3-meson channel in this file).
// T_rho^(1) (beta_rho=-0.145) and T_rho^(2) (lambda=-0.25, mu=-0.038)
// reuse the SAME unified parameter names as the KpiK/KpiPi/3pi
// channels above (gammaMag_pipi_KS, gammaMag_KpiK_vector_KS/
// deltaMag_KpiK_vector_KS) - same physical rho towers appearing yet
// again. A dedicated R$\chi$T current (Gomez Dumm & Roig, WZW+VJP+
// VVP+VPPP operators) also exists in the note but needs several
// further couplings (c_i-type odd-intrinsic-parity constants) not
// given numerically there - NOT implemented; this class is the
// simpler "old VMD/TAUOLA" current only.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// F1/F2 stub (axial vanishes by G-parity - genuine physics, not a gap).
class F1_0_EtaPiPi_Stub : public FF_0_PPP_Base {
  void    Construct(const FF_Parameters & params) {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F1_0_EtaPiPi_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};

// F4 stub (always 0, same convention as every other channel here).
class F3_0_EtaPiPi_Stub : public FF_0_PPP_Base {
  void    Construct(const FF_Parameters & params) {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F3_0_EtaPiPi_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};

class FS_0_EtaPiPi : public FF_0_PPP_Base {
  Summed_Propagator * p_Trho1, * p_Trho2;

  void    Construct(const FF_Parameters & params);
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // dedicated RChiPT current not implemented - constant fallback
  }
public:
  FS_0_EtaPiPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_Trho1(NULL), p_Trho2(NULL)
  {
    Construct(params);
  }
  ~FS_0_EtaPiPi() {
    if (p_Trho1) { delete p_Trho1; p_Trho1 = NULL; }
    if (p_Trho2) { delete p_Trho2; p_Trho2 = NULL; }
  }
};

void FS_0_EtaPiPi::Construct(const FF_Parameters & params) {
  if (m_ffmodel!=ff_model::KS) return;
  bool isPrime = (m_mode==FF_0_PPP_mode::EtaprimePiPi_pi0);
  double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  // The note's Eq.(PPPcurrent) writes the anomalous term as
  // -i c5/(4 pi^2 F^2) eps^{mu nu rho sigma} p1 p2 p3 F5, INSIDE an
  // overall normalization N.  This code has no separate N - everything
  // is folded into m_norm - so transcribing that 1/(4 pi^2 F^2)
  // literally drops the 1/F that N carries, which is what happened
  // here: the prefactor used sqr(fpi) where it needs pow(fpi,3).
  //
  // The dimensional argument is internal to this codebase and does not
  // depend on resolving the note's convention.  Every VA_0_XXX 3-meson
  // current assembles m_norm*(F1*v1 + F2*v2 + F3*q + FS*v4) with
  // v4 = cross(p1,p2,p3), so every FS_* feeding that same slot must
  // share dimensions.  FS_0_PiPlusPiZeroPiZero uses
  // Vus/(2 sqrt2 pi^2 fpi^3) (FM95 Eq.(25)); this one used fpi^2, two
  // powers, and was therefore 1/fpi = 10.8x too small in amplitude.
  //
  // Effect: Gamma(tau- -> pi- pi0 eta nu) 2.917e-17 -> 3.415e-15 GeV,
  // i.e. 0.0095x -> 1.12x the Belle measurement
  // B = (1.35 +- 0.03 +- 0.07)e-3 [arXiv:0811.0088, Inami et al].
  // N_etapipi stays a genuine O(1) mixing coefficient, default 1, as
  // it should be - it was previously absorbing a factor of 10 that was
  // never a mixing effect.
  double Nmix = (*params.p_model)(isPrime ? "N_etaprimepipi" : "N_etapipi", 1.);
  m_norm = Complex((*params.p_model)("Vud",Tools::Vud) * Nmix /
		    (4.*sqr(M_PI)*pow(fpi,3)), 0.);

  Complex gammaRho = ReadComplexParam(params.p_model,
				      "gammaMag_pipi3",-0.145,"gammaPhase_pipi3");
  Propagator_Base * rho770_1  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
  Propagator_Base * rho1450_1 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770_1,  Complex(1.,0.));
  p_Trho1->Add(rho1450_1, gammaRho);

  Complex lambda = ReadComplexParam(params.p_model,
				    "gammaMag_KpiK_vector",-0.25,"gammaPhase_KpiK_vector");
  Complex mu     = ReadComplexParam(params.p_model,
				    "deltaMag_KpiK_vector",-0.038,"deltaPhase_KpiK_vector");
  Propagator_Base * rho770_2  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
  Propagator_Base * rho1500   = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho2 = new Summed_Propagator();
  p_Trho2->Add(rho770_2, Complex(1.,0.));
  p_Trho2->Add(rho1500,  lambda);
  Total_Width_Base * wRhopp = LineShapes->Get(Flavour(kf_rho_1700_plus));
  if (wRhopp!=NULL) p_Trho2->Add(new BreitWigner(wRhopp), mu);

  std::string label = std::string("FS_0_EtaPiPi, mode=")+std::to_string(int(m_mode));
  DumpPropagatorStructure(label+" [T_rho^1]", int(m_ffmodel), p_Trho1);
  DumpPropagatorStructure(label+" [T_rho^2]", int(m_ffmodel), p_Trho2);
}

Complex FS_0_EtaPiPi::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_Trho1==NULL || p_Trho2==NULL) return Complex(0.,0.);
  // Paper's s1=(p2+p3)^2=s_pipi is reconstructed via momentum
  // conservation (same convention as every other 3-meson channel).
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  return (*p_Trho2)(s123) * (*p_Trho1)(s1);
}


DECLARE_FF_GETTER(FF_0_PPP_Base,"FF_0_PPP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PPP_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<"("<<params.m_name<<", N_f = "<<params.m_flavs.size()<<"):\n";
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    // msg_Out()<<"    *  i = "<<i<<": "<<params.m_pi[i]<<"\n";
    if (params.m_flavs[params.m_pi[i]].IsMeson()) Nmesons++;
  }
  if (Nmesons!=3) return NULL;
  // Below a first round of decays/currents for which we have both
  // Kuehn-Santamaria and RChiPT parametrizations
  msg_Out()<<"Flavours: "
	   <<params.m_flavs[params.m_pi[0]]<<" "
	   <<params.m_flavs[params.m_pi[1]]<<" "
	   <<params.m_flavs[params.m_pi[2]]<<".\n";
  if (//   pi^+ pi^0 pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi_plus) ||
      //   pi^- pi^+ pi^+
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi_plus)      
      ) {
    msg_Out()<<METHOD<<" for "<<params.m_name<<"\n";
    if (params.m_name=="F1_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F2_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="FS_0_PPP") return new FS_0_PiPlusPiZeroPiZero(params);
  }
  //   pi^0 pi^0 K^- (FM95 hep-ph/9503474) - reuses the same classes
  //   as above; F1_0_PiPlusPiZeroPiZero/FS_0_PiPlusPiZeroPiZero switch
  //   on m_mode internally to pick the K1/K* resonances for this case.
  if (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi &&
      params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi &&
      params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_plus) {
    if (params.m_name=="F1_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F2_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="FS_0_PPP") return new FS_0_PiPlusPiZeroPiZero(params);
  }
  // K^+ K^- pi^- (cf. VA_0_KKPi.H for the ordering convention).
  // F1, F2, FS(=F5) have real RChL2012 implementations (see F1_0_KKPi,
  // F5_0_KKPi above); F3 is identically 0 for this channel regardless
  // of model (c3=0 in Table 1 of 1203.3955 / Table 1 of 1509.09140),
  // so it stays a stub for every ff_model.
  if (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
      params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_plus &&
      params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi_plus) {
    if (params.m_name=="F1_0_KKPi") return new F1_0_KKPi(params);
    if (params.m_name=="F2_0_KKPi") return new F1_0_KKPi(params);
    if (params.m_name=="F3_0_KKPi") return new F_0_KKPi_Stub(params);
    if (params.m_name=="FS_0_KKPi") return new F5_0_KKPi(params);
  }
  //////////////////////////////////////////////////////////////////
  // FM95 (hep-ph/9503474) "KS"-mode KpiK family: K^-pi^-K^+,
  // K^0pi^-K0bar, K_Spi^-K_S, K_Spi^-K_L, K_Lpi^-K_L, K^-pi^0K^0.
  // Wired to the new VA_0_KPiK Current class (params.m_name="*_0_KPiK"),
  // NOT VA_0_KKPi (different momentum-index convention - see the
  // KM_piM_KP comment at FixMode() above).
  //////////////////////////////////////////////////////////////////
  if ((params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_plus) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_S) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_S &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_L) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_S) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_L &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_L) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi &&
       (params.m_flavs[params.m_pi[2]].Kfcode()==kf_K ||
	params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_S ||
	params.m_flavs[params.m_pi[2]].Kfcode()==kf_K_L))) {
    if (params.m_name=="F1_0_KPiK") return new F1_0_KPiK(params);
    if (params.m_name=="F2_0_KPiK") return new F1_0_KPiK(params);
    if (params.m_name=="F3_0_KPiK") return new F3_0_KPiK_Stub(params);
    if (params.m_name=="FS_0_KPiK") return new FS_0_KPiK(params);
  }
  //////////////////////////////////////////////////////////////////
  // FM95 KpiPi family: K^-pi^-pi^+, pi^-K0bar pi^0. Wired to the new
  // VA_0_KPiPi Current class.
  //////////////////////////////////////////////////////////////////
  if ((params.m_flavs[params.m_pi[0]].Kfcode()==kf_K_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi_plus) ||
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       (params.m_flavs[params.m_pi[1]].Kfcode()==kf_K ||
	params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_S ||
	params.m_flavs[params.m_pi[1]].Kfcode()==kf_K_L) &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi)) {
    if (params.m_name=="F1_0_KPiPi") return new F1_0_KPiPi(params);
    if (params.m_name=="F2_0_KPiPi") return new F1_0_KPiPi(params);
    if (params.m_name=="F3_0_KPiPi") return new F3_0_KPiPi_Stub(params);
    if (params.m_name=="FS_0_KPiPi") return new FS_0_KPiPi(params);
  }
  //////////////////////////////////////////////////////////////////
  // eta pi^- pi^0 / eta' pi^- pi^0. Exact order required (eta(') is
  // NOT interchangeable with either pion - see the FixMode()/enum
  // comments). Wired to a new VA_0_EtaPiPi Current class.
  //////////////////////////////////////////////////////////////////
  if ((params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta ||
       params.m_flavs[params.m_pi[0]].Kfcode()==kf_eta_prime_958) &&
      params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
      params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi) {
    if (params.m_name=="F1_0_EtaPiPi") return new F1_0_EtaPiPi_Stub(params);
    if (params.m_name=="F2_0_EtaPiPi") return new F1_0_EtaPiPi_Stub(params);
    if (params.m_name=="F3_0_EtaPiPi") return new F3_0_EtaPiPi_Stub(params);
    if (params.m_name=="FS_0_EtaPiPi") return new FS_0_EtaPiPi(params);
  }
  return NULL;
}




