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
  m_norm(Complex(0.,0.)) {
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
  // K^-(q1) pi^0(q2) K^0(q3) - row 5.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi &&
	   m_flavs[m_pi[2]].Kfcode()==kf_K)          m_mode = FF_0_PPP_mode::KM_pi0_K0;
  // K^-(q1) pi^-(q2) pi^+(q3) - row 7.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_K_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi_plus)    m_mode = FF_0_PPP_mode::KM_piM_piP;
  // pi^-(q1) K0bar(q2) pi^0(q3) - row 8.
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_K &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi)         m_mode = FF_0_PPP_mode::piM_K0bar_pi0;
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
  case ff_model::KS:      return m_norm * FF_KS(s123,s12,s13);
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
  double  m_xiK1;

  // --- RChL2012 (1203.3955 Sec.2.1 + 1310.1053 sigma extension) ---
  // Defaults below are the BaBar-fitted values of Table I of
  // 1310.1053 (pi- pi- pi+ channel), used for BOTH 3pi charge modes
  // per that paper's explicit recommendation (Sec. "Case of pi0pi0pi-
  // mode"): only alpha_sigma, gamma_sigma differ (scaled by 0.63 for
  // the neutral mode), everything else stays as in Table I.
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
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s1,const double & s2);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2);
  Complex FF_RChL2012(const double & s123,const double & s1,const double & s2);
  /*
    Complex A(const double & m2,const double & s,const double & mu2);
    double  Gamma_V(const double & s);
    double  Gamma_Vp(const double & s);
    double  Gamma_Vpp(const double & s);
  */
public :
  F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params);
  ~F1_0_PiPlusPiZeroPiZero();
};

F1_0_PiPlusPiZeroPiZero::F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
  FF_0_PPP_Base(params),
  p_a1s(NULL), p_rhos(NULL),
  p_TK1a(NULL), p_TKstar1_pi0K(NULL), m_xiK1(0.33),
  m_isF2(false), 
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  FixParameters(params);
  Construct();
}

F1_0_PiPlusPiZeroPiZero::~F1_0_PiPlusPiZeroPiZero() {
  if (p_a1s)  { delete p_a1s;  p_a1s  = NULL; }
  if (p_rhos) { delete p_rhos; p_rhos = NULL; }
  if (p_TK1a)         { delete p_TK1a;         p_TK1a         = NULL; }
  if (p_TKstar1_pi0K) { delete p_TKstar1_pi0K; p_TKstar1_pi0K = NULL; }
}

void F1_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {
  if (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ||
      m_mode==FF_0_PPP_mode::piM_piP_piP) {
    m_norm    = Complex(0., -((2.*sqrt(2)*(*params.p_model)("Vud", Tools::Vud)) /
			      (3.*m_fpi) ));
    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.14500,0.0000);
      m_delta = Complex( 0.00000,0.0000);
      m_alpha = Complex( 0.00185,0.0000);
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
      m_lambda_pp  = RChL::Lambda_pp(m_FV_rchl,m_GV_rchl,m_lambda_p);
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
  }
  else if (m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // Finkemeier-Mirkes hep-ph/9503474 (FM95), Tab.I row 6:
    // A^(abc) = sin(theta_c)/4, G1=T_K1a(Q^2)T_K*^(1)(s2),
    // G2=T_K1a(Q^2)T_K*^(1)(s1). Eq.(23)/(24): F_i = (2sqrt2 A/3fpi) G_i,
    // Cabibbo-SUPPRESSED (sin, not cos) since this is a |Delta S|=1
    // (one-kaon) channel - use Vus, not Vud.
    if (m_ffmodel==ff_model::KS) {
      m_norm = Complex((2.*sqrt(2.)*(*params.p_model)("Vus",Tools::Vus)*
			(*params.p_model)("sinThetaC_over_4_num",1.)/4.) /
			(3.*m_fpi), 0.);
      // xi: relative K1(1270) admixture in T_K1^(a), Eq.(32)-(33) of
      // FM95: |xi|=0.33, sign preferred by data is xi=+0.33 (Sec.VII).
      m_xiK1 = (*params.p_model)("xiK1",0.33);
    }
  }
}

void F1_0_PiPlusPiZeroPiZero::Construct() {
  if (m_ffmodel==ff_model::KS && m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // T_K1^(a) = [BW_K1(1400) + xi*BW_K1(1270)] / (1+xi), Eq.(32).
    // K1(1270)/K1(1400) lineshapes now built and registered
    // (K1_Decays.H/.C, kf 10313/10323 and 20313/20323, confirmed
    // against the real ATOOLS/Phys/Flavour_Tags.H), using the
    // proper running widths from their rho-K/K*-pi decay channels
    // rather than FM95's own simple fixed-width prescription -
    // a reasonable, arguably better-motivated substitution.
    Total_Width_Base * wK11270 = LineShapes->Get(Flavour(kf_K_1_1270_plus));
    Total_Width_Base * wK11400 = LineShapes->Get(Flavour(kf_K_1_1400_plus));
    if (wK11270==NULL || wK11400==NULL) {
      msg_Error()<<"Error in "<<METHOD<<": missing K1(1270)/K1(1400) "
		 <<"lineshape(s) for pi0_pi0_KM (FM95) - T_K1^(a) will be "
		 <<"treated as identically zero until these are added to "
		 <<"Line_Shapes::Init(). See the FixParameters()/Construct() "
		 <<"comments in FF_0_PPP.C for what is needed.\n";
    }
    else {
      Propagator_Base * K11400 = new BreitWigner(wK11400);
      Propagator_Base * K11270 = new BreitWigner(wK11270);
      p_TK1a = new Summed_Propagator();
      p_TK1a->Add(K11400, Complex(1.,0.));
      p_TK1a->Add(K11270, Complex(m_xiK1,0.));
    }
    // T_K*^(1) = [BW_K*(892) + beta_K* BW_K*(1410)]/(1+beta_K*),
    // Eq.(10), beta_K*=-0.135 (Sec.II). Both lineshapes already exist
    // (used by FF_0_PP.C's Kpi_plus mode) - reused here directly.
    Propagator_Base * Kstar892  =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
    Propagator_Base * Kstar1410 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    p_TKstar1_pi0K = new Summed_Propagator();
    p_TKstar1_pi0K->Add(Kstar892,  Complex(1., 0.));
    p_TKstar1_pi0K->Add(Kstar1410, Complex(-0.135,0.));
  }
  else if (m_ffmodel==ff_model::KS) {
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
    // a1 is NOT built into p_a1s here (stays NULL for this model): its
    // propagator needs a q^2 numerator and the Eq.(17)-of-1310.1053
    // polynomial width fit, neither of which fit Propagator_Base's
    // existing hierarchy - it stays inline in FF_RChL2012 below.
    //
    // Charged/neutral rho choice matches the KS branch above: for
    // pi0 pi0 pi-, the rho in each pi0-pi- sub-system is CHARGED
    // (rho- -> pi0 pi-); for pi- pi- pi+, the rho in the pi-pi+
    // sub-system is NEUTRAL (rho0 -> pi+ pi-, hence also mixing with
    // omega there in the KS branch). This was backwards in an earlier
    // version of this branch - fixed here.
    Flavour rhoFlav = (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ?
		      Flavour(kf_rho_770_plus) : Flavour(kf_rho_770));
    Propagator_Base * rho  = new RChL_BW(LineShapes->Get(rhoFlav));
    Propagator_Base * rhop = new RChL_BW(LineShapes->Get(Flavour(kf_rho_1450_plus)));
    p_rhos = new Summed_Propagator();
    p_rhos->Add(rho,  Complex(1.,0.));
    p_rhos->Add(rhop, Complex(m_betarhop,0.));
  }
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
  if (m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // FM95 Tab.I row 6: G1(Q^2,s2,s3)=T_K1a(Q^2)T_K*^(1)(s2),
    // G2(Q^2,s1,s3)=T_K1a(Q^2)T_K*^(1)(s1) - same s1-reconstruction
    // as the pion modes above (s2=s13, s1 needs momentum conservation).
    if (p_TKstar1_pi0K==NULL) return Complex(0.,0.);
    Complex TK1a = (p_TK1a!=NULL ? (*p_TK1a)(s123) : Complex(0.,0.));
    if (m_isF2) {
      double s1 = s123 - s12 - s13 +
	m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
      return m_norm * TK1a * (*p_TKstar1_pi0K)(s1);
    }
    return m_norm * TK1a * (*p_TKstar1_pi0K)(s13);
  }
  if (p_a1s==NULL || p_rhos==NULL) return Complex(0.,0.);
  if (m_isF2) {
    double s1 = s123 - s12 - s13 +
      m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
    return m_norm * (*p_a1s)(s123) * (*p_rhos)(s1);
  }
  return m_norm * (*p_a1s)(s123) * (*p_rhos)(s13);
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
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
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
    FixParameters(params);
    Construct();
  }
  ~F3_0_PiPlusPiZeroPiZero() {
    if (p_rho_combo) { delete p_rho_combo; p_rho_combo = NULL; }
  }
};

void F3_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {
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
  }
}

void F3_0_PiPlusPiZeroPiZero::Construct() {
  if (m_ffmodel==ff_model::RChL2012) {
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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
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
    FixParameters(params);
    Construct();
  }
  ~FS_0_PiPlusPiZeroPiZero() {
    if (p_TKstar1_v) { delete p_TKstar1_v; p_TKstar1_v = NULL; }
    if (p_TKstar2)   { delete p_TKstar2;   p_TKstar2   = NULL; }
  }
};

void FS_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {
  if (m_ffmodel==ff_model::KS && m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    // FM95 Tab.II row 6: A^(abc)=sin(theta_c), F3^(abc)=A/(2sqrt2 pi^2
    // fpi^3) G3, Eq.(25). Cabibbo-suppressed (Vus), |Delta S|=1.
    m_norm = Complex((*params.p_model)("Vus",Tools::Vus) /
		      (2.*sqrt(2.)*sqr(M_PI)*pow(m_fpi,3)), 0.);
  }
}

void FS_0_PiPlusPiZeroPiZero::Construct()  {
  if (m_ffmodel==ff_model::KS && m_mode==FF_0_PPP_mode::pi0_pi0_KM) {
    Propagator_Base * Kstar892_1 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
    Propagator_Base * Kstar1410_1 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    p_TKstar1_v = new Summed_Propagator();
    p_TKstar1_v->Add(Kstar892_1,  Complex(1., 0.));
    p_TKstar1_v->Add(Kstar1410_1, Complex(-0.135,0.));
    // T_K*^(2), Eq.(42): K*(892) + lambda*K*(1410) + mu*K*(1714)/K*(1680).
    // RESOLVED: K*(1680) lineshape now added (Kstar_Decays.H/.C, kf
    // 30313/30323), registered in Line_Shapes::Init(). Its pole mass
    // comes from whatever Sherpa's own particle table has for that
    // kf-code, not necessarily either paper's quoted number (1.714 GeV
    // per an earlier guess here vs 1.700 GeV per
    // tau_two_meson_currents_KS_RChiT.tex's CLEO-tune usage) - the
    // NULL-check below is now just a defensive fallback (in case the
    // particle isn't registered in some build), not an expected path.
    Total_Width_Base * wKstarpp = LineShapes->Get(Flavour(kf_K_star_1680_plus));
    p_TKstar2 = new Summed_Propagator();
    Propagator_Base * Kstar892_2 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
    Propagator_Base * Kstar1410_2 =
      new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
    p_TKstar2->Add(Kstar892_2,  Complex(1.,0.));
    p_TKstar2->Add(Kstar1410_2, Complex(-0.25,0.));
    if (wKstarpp==NULL) {
      msg_Error()<<"Error in "<<METHOD<<": missing K*(1714)/K*'' lineshape "
		 <<"for pi0_pi0_KM vector form factor (FM95 Eq.42) - "
		 <<"dropping the mu*BW_K*'' term (mu=-0.038, small but not "
		 <<"zero) until this resonance is added to Line_Shapes.\n";
    }
    else {
      Propagator_Base * Kstarpp = new BreitWigner(wKstarpp);
      p_TKstar2->Add(Kstarpp, Complex(-0.038,0.));
    }
  }
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
  void    FixParameters(const FF_Parameters & params) {}
  void    Construct() {}
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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
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
    FixParameters(params);
    Construct();
  }
};

void F1_0_KKPi::FixParameters(const FF_Parameters & params) {
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
}

void F1_0_KKPi::Construct() {
  p_GRho = LineShapes->Get(Flavour(kf_rho_770_plus));
  p_GKst = LineShapes->Get(Flavour(kf_K_star_892_plus));
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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
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
    FixParameters(params);
    Construct();
  }
};

void F5_0_KKPi::FixParameters(const FF_Parameters & params) {
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
}

void F5_0_KKPi::Construct() {
  p_GRho   = LineShapes->Get(Flavour(kf_rho_770_plus));
  p_GKst   = LineShapes->Get(Flavour(kf_K_star_892_plus));
  p_GOmega = LineShapes->Get(Flavour(kf_omega_782));
  // No phi(1020) lineshape exists in Line_Shapes.C yet; at the default
  // ideal mixing angle its contribution vanishes identically in Eq.(15)
  // anyway (cos^2(thetaV)(1-sqrt(2)tan(thetaV))=0 at thetaV=35.26deg),
  // so p_GPhi staying NULL is harmless UNLESS thetaV is changed away
  // from ideal mixing via thetaV_rchlKKpi - in that case this term is
  // silently dropped rather than contributing. Flagged here rather than
  // silently wrong: add a phi(1020) Total_Width_Base if non-ideal
  // mixing is ever needed.
  p_GPhi   = NULL;
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
  Complex omegaPhiTerm = sin2V*mixOmega*propOmega; // phi term dropped, see Construct()

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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
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
    FixParameters(params);
    Construct();
  }
  ~F1_0_KPiK() {
    if (p_BWA1)    { delete p_BWA1;    p_BWA1    = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
  }
};

void F1_0_KPiK::FixParameters(const FF_Parameters & params) {
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
    m_norm = Complex( pref*3.*Vud/(2.*sqrt(2.)), 0.); break;   // row 5
  default: m_norm = Complex(0.,0.); break;
  }
}

void F1_0_KPiK::Construct() {
  if (m_ffmodel!=ff_model::KS) return;
  p_BWA1 = new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
  // T_rho^(1): only needed for the rows that actually reference it
  // (1,2,4,5) - build unconditionally, it is cheap and harmless for
  // rows 3/KL that don't use it.
  Propagator_Base * rho770  =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
  Propagator_Base * rho1450 =
    new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, Complex(-0.145,0.)); // beta_rho, Eq.(9)
  Propagator_Base * Kstar892  =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 =
    new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, Complex(-0.135,0.)); // beta_K*, Eq.(10)
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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  FS_0_KPiK(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_Trho2(NULL), p_TKstar1(NULL), p_Tomega(NULL)
  {
    FixParameters(params);
    Construct();
  }
  ~FS_0_KPiK() {
    if (p_Trho2)   { delete p_Trho2;   p_Trho2   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
    if (p_Tomega)  { delete p_Tomega;  p_Tomega  = NULL; }
  }
};

void FS_0_KPiK::FixParameters(const FF_Parameters & params) {
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
  case FF_0_PPP_mode::KM_pi0_K0:    m_norm = Complex(-pref*Vud/sqrt(2.),0.); break;
  default: m_norm = Complex(0.,0.); break;
  }
}

void FS_0_KPiK::Construct() {
  if (m_ffmodel!=ff_model::KS) return;
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
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
  Propagator_Base * rho1500 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho2 = new Summed_Propagator();
  p_Trho2->Add(rho770,  Complex(1.,0.));
  p_Trho2->Add(rho1500, Complex(-0.25,0.)); // lambda
  Total_Width_Base * wRhopp = LineShapes->Get(Flavour(kf_rho_1700_plus));
  if (wRhopp!=NULL) p_Trho2->Add(new BreitWigner(wRhopp), Complex(-0.038,0.)); // mu

  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, Complex(-0.135,0.));

  // T_omega, Eq.(36): omega(782)+eps phi(1020). phi(1020) confirmed
  // and registered (Omega_Decays.H/.C, kf_phi_1020=333).
  Propagator_Base * omega782 = new BreitWigner(LineShapes->Get(Flavour(kf_omega_782)));
  p_Tomega = new Summed_Propagator();
  p_Tomega->Add(omega782, Complex(1.,0.));
  Total_Width_Base * wPhi = LineShapes->Get(Flavour(kf_phi_1020));
  if (wPhi!=NULL) p_Tomega->Add(new BreitWigner(wPhi), Complex(0.05,0.)); // eps
  else msg_Error()<<"Error in "<<METHOD<<": phi(1020) lineshape "
		  <<"unexpectedly unavailable for T_omega (FM95 Eq.36) - "
		  <<"falling back to plain BW_omega(782), a ~5% approximation.\n";
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
  void    FixParameters(const FF_Parameters & params) {}
  void    Construct() {}
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
  double  m_xiK1;
  Summed_Propagator * p_TK1a, * p_Trho1, * p_TKstar1;
  Propagator_Base    * p_TK1b;

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  F1_0_KPiPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), m_isF2(false), m_xiK1(0.33),
    p_TK1a(NULL), p_Trho1(NULL), p_TKstar1(NULL), p_TK1b(NULL)
  {
    if (params.m_name=="F2_0_KPiPi") m_isF2 = true;
    FixParameters(params);
    Construct();
  }
  ~F1_0_KPiPi() {
    if (p_TK1a)    { delete p_TK1a;    p_TK1a    = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
    if (p_TK1b)    { delete p_TK1b;    p_TK1b    = NULL; }
  }
};

void F1_0_KPiPi::FixParameters(const FF_Parameters & params) {
  if (m_ffmodel!=ff_model::KS) return;
  double fpi = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vus = (*params.p_model)("Vus",Tools::Vus);
  double pref = 2.*sqrt(2.)/(3.*fpi);
  m_xiK1 = (*params.p_model)("xiK1",0.33);
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_piP:
    m_norm = Complex(-pref*Vus/2., 0.); break;                   // row 7
  case FF_0_PPP_mode::piM_K0bar_pi0:
    m_norm = Complex( pref*3.*Vus/(2.*sqrt(2.)), 0.); break;     // row 8
  default: m_norm = Complex(0.,0.); break;
  }
}

void F1_0_KPiPi::Construct() {
  if (m_ffmodel!=ff_model::KS) return;
  Total_Width_Base * wK11400 = LineShapes->Get(Flavour(kf_K_1_1400_plus));
  Total_Width_Base * wK11270 = LineShapes->Get(Flavour(kf_K_1_1270_plus));
  if (wK11400==NULL || wK11270==NULL) {
    msg_Error()<<"Error in "<<METHOD<<": missing K1(1270)/K1(1400) "
	       <<"lineshape(s) - T_K1^(a) and T_K1^(b) will be treated as "
	       <<"identically zero (see F1_0_PiPlusPiZeroPiZero's Construct() "
	       <<"for what is needed to fix this).\n";
  }
  else {
    p_TK1a = new Summed_Propagator();
    p_TK1a->Add(new BreitWigner(wK11400), Complex(1.,0.));
    p_TK1a->Add(new BreitWigner(wK11270), Complex(m_xiK1,0.));
    p_TK1b = new BreitWigner(wK11270);
  }
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
  Propagator_Base * rho1450 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, Complex(-0.145,0.));
  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410, Complex(-0.135,0.));
}

Complex F1_0_KPiPi::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_Trho1==NULL || p_TKstar1==NULL) return Complex(0.,0.);
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
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

  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s12,const double & s13);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(1.,0.); // RChiPT not implemented for this channel - constant fallback
  }
public:
  FS_0_KPiPi(const FF_Parameters & params) :
    FF_0_PPP_Base(params), p_TKstar2(NULL), p_Trho1(NULL), p_TKstar1(NULL)
  {
    FixParameters(params);
    Construct();
  }
  ~FS_0_KPiPi() {
    if (p_TKstar2) { delete p_TKstar2; p_TKstar2 = NULL; }
    if (p_Trho1)   { delete p_Trho1;   p_Trho1   = NULL; }
    if (p_TKstar1) { delete p_TKstar1; p_TKstar1 = NULL; }
  }
};

void FS_0_KPiPi::FixParameters(const FF_Parameters & params) {
  if (m_ffmodel!=ff_model::KS) return;
  double fpi  = (*params.p_model)("fpi",0.1307)/sqrt(2.);
  double Vus  = (*params.p_model)("Vus",Tools::Vus);
  double pref = 1./(2.*sqrt(2.)*sqr(M_PI)*pow(fpi,3));
  switch (m_mode) {
  case FF_0_PPP_mode::KM_piM_piP:    m_norm = Complex(pref*Vus,0.);          break;
  case FF_0_PPP_mode::piM_K0bar_pi0: m_norm = Complex(pref*sqrt(2.)*Vus,0.); break;
  default: m_norm = Complex(0.,0.); break;
  }
}

void FS_0_KPiPi::Construct() {
  if (m_ffmodel!=ff_model::KS) return;
  Propagator_Base * Kstar892  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar2 = new Summed_Propagator();
  p_TKstar2->Add(Kstar892,  Complex(1.,0.));
  p_TKstar2->Add(Kstar1410, Complex(-0.25,0.));
  Total_Width_Base * wKstarpp = LineShapes->Get(Flavour(kf_K_star_1680_plus));
  if (wKstarpp!=NULL) p_TKstar2->Add(new BreitWigner(wKstarpp), Complex(-0.038,0.));
  else msg_Error()<<"Error in "<<METHOD<<": missing K*(1714)/K*'' lineshape "
		  <<"(FM95 Eq.42) - dropping the mu-term, see the identical "
		  <<"flag in FS_0_PiPlusPiZeroPiZero above.\n";
  Propagator_Base * rho770  = new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
  Propagator_Base * rho1450 = new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
  p_Trho1 = new Summed_Propagator();
  p_Trho1->Add(rho770,  Complex(1.,0.));
  p_Trho1->Add(rho1450, Complex(-0.145,0.));
  Propagator_Base * Kstar892_1  = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_892_plus)));
  Propagator_Base * Kstar1410_1 = new BreitWigner(LineShapes->Get(Flavour(kf_K_star_1410_plus)));
  p_TKstar1 = new Summed_Propagator();
  p_TKstar1->Add(Kstar892_1,  Complex(1., 0.));
  p_TKstar1->Add(Kstar1410_1, Complex(-0.135,0.));
}

Complex FS_0_KPiPi::
FF_KS(const double & s123,const double & s12,const double & s13) {
  if (p_TKstar2==NULL || p_Trho1==NULL || p_TKstar1==NULL) return Complex(0.,0.);
  double s1 = s123 - s12 - s13 +
    m_masses2[m_pi[0]] + m_masses2[m_pi[1]] + m_masses2[m_pi[2]];
  Complex TKstar2 = (*p_TKstar2)(s123);
  switch (m_mode) {
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
  void    FixParameters(const FF_Parameters & params) {}
  void    Construct() {}
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return Complex(0.,0.);
  }
public:
  F3_0_KPiPi_Stub(const FF_Parameters & params) : FF_0_PPP_Base(params) {}
};


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
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_K)) {
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
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_K &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi)) {
    if (params.m_name=="F1_0_KPiPi") return new F1_0_KPiPi(params);
    if (params.m_name=="F2_0_KPiPi") return new F1_0_KPiPi(params);
    if (params.m_name=="F3_0_KPiPi") return new F3_0_KPiPi_Stub(params);
    if (params.m_name=="FS_0_KPiPi") return new FS_0_KPiPi(params);
  }
  return NULL;
}




