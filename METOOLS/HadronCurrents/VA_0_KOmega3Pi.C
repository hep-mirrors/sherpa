#include "METOOLS/HadronCurrents/VA_0_KOmega3Pi.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// K-omega -> K pi+ pi- pi0 (tau_two_meson_currents_KS_RChiT.tex,
// Sec."K omega -> K pi+pi-pi0: four-pseudoscalar channel"). See
// VA_0_KOmega3Pi.H for the class-level overview. Full derivation of
// the implemented formula:
//
// Generic PV production current (Eq.PV-current), with V=omega,
// P=K, p1=p_omega (=k, the omega's own momentum = p_++p_-+p_0 by
// momentum conservation), p2=p_K, Q=p1+p2:
//   H_mu = -g(Q^2) eps_{mu,nu,alpha,beta} eps_omega*^nu p1^alpha p2^beta
//          - i f(Q^2) eps_omega*_mu
//          - i{a1(Q^2) p1_mu + a2(Q^2) p2_mu}(eps_omega*.Q)
// Stripping off eps_omega*^nu (replaced by an open index rho, since
// it will be contracted through the omega propagator instead of a
// final-state polarisation vector) gives the two-index object
//   H^{mu,rho} = -g(Q^2) eps^{mu,rho,alpha,beta} p1_alpha p2_beta
//                - i f(Q^2) g^{mu,rho}
//                - i{a1(Q^2) p1^mu + a2(Q^2) p2^mu} Q^rho.
//
// The omega->3pi anomalous vertex (Eq.omega3pi) contributes an
// "omega-decay vector"
//   W^rho = i G_omega3pi(s_+-,s_+0,s_-0) eps^{rho,beta,gamma,delta}
//           p_+_beta p_-_gamma p_0_delta,
// which is IDENTICALLY TRANSVERSE to k=p1=p_++p_-+p_0 (k.W=0 always,
// since eps is antisymmetric and any term in k.W repeats one of
// p_+,p_-,p_0 in two Levi-Civita slots). Contracting W through the
// massive-vector propagator (Eq.PV-vector-prop),
//   P^{rho,alpha}_omega(k) = [-g^{rho,alpha}+k^rho k^alpha/M_omega^2]
//                            / D_omega(k^2),
// the k^rho k^alpha/M_omega^2 piece drops out identically because
// k.W=0, leaving simply
//   P_{omega,rho,alpha}(k) W^alpha = -W_rho / D_omega(k^2)
//                                  =: Vtilde_rho.
// So the full current reduces to ordinary 4-vector arithmetic:
//   J^mu = H^{mu,rho}(Q^2) Vtilde_rho
//        = -g(Q^2) eps^{mu,rho,alpha,beta} p1_alpha p2_beta Vtilde_rho
//          - i f(Q^2) Vtilde^mu
//          - i{a1(Q^2)p1^mu+a2(Q^2)p2^mu}(Q.Vtilde).
// The first term is itself a (real p1, real p2, complex Vtilde)
// triple-vector Levi-Civita contraction; since Vtilde is a complex
// SCALAR prefactor times the real vector X=cross(p_+,p_-,p_0) (see
// below), linearity of the Levi-Civita contraction means this can be
// computed as (complex prefactor) * cross(p1,p2,X) using the same
// real-valued cross() helper twice - no genuine complex Levi-Civita
// or rank-2 tensor code is needed anywhere in this implementation.
//
///////////////////////////////////////////////////////////////////////////

VA_0_KOmega3Pi::VA_0_KOmega3Pi(const ATOOLS::Flavour_Vector& flavs,
			       const std::vector<int>& indices,
			       const std::string& name) :
  Current_Base(flavs, indices, name),
  m_MK(0.), m_MK1a(0.), m_GK1a(0.), m_MK1b(0.), m_GK1b(0.),
  m_MKsta(0.), m_GKsta(0.), m_MKstb(0.), m_GKstb(0.),
  m_Momega(0.), p_omega_width(NULL), p_rho_width(NULL),
  m_hKtKwK(0.,0.), m_hA1tA1wK(0.,0.), m_hA2tA2wK(0.,0.),
  m_hV1tV1wK(0.,0.), m_hV2tV2wK(0.,0.), m_Nomega3pi(1.,0.),
  m_norm(1.)
{
  msg_Out()<<METHOD<<"(N_f = "<<m_flavs.size()<<"):\n";
  for (size_t i=0;i<p_i.size();i++) {
    msg_Out()<<"    *  i = "<<i<<": "<<p_i[i]<<"  --> "
	     <<m_flavs[p_i[i]]<<".\n";
  }
}

VA_0_KOmega3Pi::~VA_0_KOmega3Pi() {}

Complex VA_0_KOmega3Pi::
DRfixed(const double & Q2,const double & M,const double & G) const {
  // D_R(Q^2) = Q^2-M_R^2+iM_R Gamma_R, returning 1/D_R.
  //
  // FIXED width, as arXiv:0709.4039 Eq.(8) writes it (no Q^2
  // dependence on Gamma_R), and as its alpha_omegaK fit to
  // B(tau->omega K nu) was performed - so the fitted coupling only
  // means what it says in this convention.
  //
  // Registry running widths were tried here and are WORSE: every pole
  // (K1 at 1.27/1.40, K* at 0.89/1.41) lies below the kinematic range,
  // which runs to m_tau^2 = 3.16 GeV^2, so the form factors are probed
  // on the high-s side where a running Gamma(s) grows and suppresses
  // 1/|D|^2 in the region carrying the rate. Measured 2026-08-27:
  // B(tau->K-omega nu)/CLEO fell 0.536 -> 0.193 for the default
  // configuration. Do not swap these for registry shapes without
  // refitting alpha_omegaK alongside.
  return 1./Complex(Q2-sqr(M), M*G);
}

Complex VA_0_KOmega3Pi::FF_f(const double & Q2) const {
  // Eq.(PV-f): f(Q^2) = -1/2(Q^2+M_omega^2-m_K^2) sum_i hA_i tA_i/D_Ai(Q^2)
  Complex SA = m_hA1tA1wK*DRfixed(Q2,m_MK1a,m_GK1a)
             + m_hA2tA2wK*DRfixed(Q2,m_MK1b,m_GK1b);
  return -0.5*(Q2+sqr(m_Momega)-sqr(m_MK))*SA;
}

Complex VA_0_KOmega3Pi::FF_g(const double & Q2) const {
  // g(Q^2) = sum_i hV_i tV_i/D_Vi(Q^2) - NOTE the 1/2 printed in
  // arXiv:1408.0086 Eq.(14) is dropped deliberately.
  //
  // That 1/2 is inherited from the original source, Flores-Tlalpa and
  // Lopez Castro, PRD 77 (2008) 113011 [arXiv:0709.4039], whose Eq.(5)
  // writes the vector term in the q_+/q_- basis (q_pm = p_V +- p_P):
  //   i g eps^{abmn} eps*_b q_{+m} q_{-n},   g = 1/2 sum_j f_Vj g_VjVP/D_Vj.
  // Since eps^{abmn} q_{+m} q_{-n} = -2 eps^{abmn} p_1m p_2n, rewriting
  // in the p_1,p_2 basis used by 1408.0086 Eq.(13) - the basis this
  // class implements - carries a factor 2, so g_here = 2 g_0709.4039
  // = sum_i, not 1/2 sum_i. 1408.0086 changed the tensor basis but kept
  // the 1/2, which is inconsistent with its OWN Eq.(25): x_0 =
  // -1/4 g^2 lambda - f^2 reproduces 0709.4039's Eq.(7) alpha =
  // |f|^2 + lambda |g|^2 only for g_here = 2 g_0709.4039.
  // The a_1, a_2 coefficients cross-check the same mapping exactly:
  // 0709.4039 Eq.(8) has a_+ = 1/2 S_A + 2 S_P and a_- = 1/2 S_P, and
  // with q_+ = Q, q_- = p_1 - p_2 the bracket is (a_+ + a_-)p_1 +
  // (a_+ - a_-)p_2, giving 5/2 and 3/2 - exactly Eq.(14)'s a_1, a_2.
  Complex SV = m_hV1tV1wK*DRfixed(Q2,m_MKsta,m_GKsta)
             + m_hV2tV2wK*DRfixed(Q2,m_MKstb,m_GKstb);
  return SV;
}

Complex VA_0_KOmega3Pi::FF_a1(const double & Q2) const {
  // Eq.(PV-a1): a1 = 5/2 hK tK/D_K(Q^2) + 1/2 sum_i hA_i tA_i/D_Ai(Q^2)
  Complex SP = m_hKtKwK*DRfixed(Q2,m_MK,0.); // Gamma_K=0, stable
  Complex SA = m_hA1tA1wK*DRfixed(Q2,m_MK1a,m_GK1a)
             + m_hA2tA2wK*DRfixed(Q2,m_MK1b,m_GK1b);
  return 2.5*SP + 0.5*SA;
}

Complex VA_0_KOmega3Pi::FF_a2(const double & Q2) const {
  // Eq.(PV-a2): a2 = 3/2 hK tK/D_K(Q^2) + 1/2 sum_i hA_i tA_i/D_Ai(Q^2)
  Complex SP = m_hKtKwK*DRfixed(Q2,m_MK,0.);
  Complex SA = m_hA1tA1wK*DRfixed(Q2,m_MK1a,m_GK1a)
             + m_hA2tA2wK*DRfixed(Q2,m_MK1b,m_GK1b);
  return 1.5*SP + 0.5*SA;
}

Complex VA_0_KOmega3Pi::
Gomega3pi(const double & spm,const double & sp0,const double & sm0) const {
  // Eq.(omega3pi-vmd), pure rho-exchange (contact term g_omega3pi^ct
  // set to 0, "the pure rho-exchange implementation" per the note).
  // BW_rho(s) here uses Sherpa's own registered rho(770) RUNNING
  // width (bare propagator form, matching this equation's own D_rho-
  // style convention - NOT the M^2-normalised BreitWigner class used
  // elsewhere in this codebase, to avoid an unwanted extra M^2
  // numerator factor creeping in).
  if (p_rho_width==NULL) return Complex(0.,0.);
  double Mrho = p_rho_width->Flav().HadMass();
  Complex BWpm = 1./Complex(spm-sqr(Mrho), Mrho*(*p_rho_width)(spm));
  Complex BWp0 = 1./Complex(sp0-sqr(Mrho), Mrho*(*p_rho_width)(sp0));
  Complex BWm0 = 1./Complex(sm0-sqr(Mrho), Mrho*(*p_rho_width)(sm0));
  return m_Nomega3pi*(BWpm+BWp0+BWm0);
}

void VA_0_KOmega3Pi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D pK   = moms[p_i[0]];
  Vec4D pPiP = moms[p_i[1]];
  Vec4D pPiM = moms[p_i[2]];
  Vec4D pPi0 = moms[p_i[3]];
  Vec4D p1   = pPiP+pPiM+pPi0;   // omega momentum, k
  Vec4D p2   = pK;
  Vec4D Q    = p1+p2;
  double Q2  = Q.Abs2();
  double k2  = p1.Abs2();

  double spm = (pPiP+pPiM).Abs2();
  double sp0 = (pPiP+pPi0).Abs2();
  double sm0 = (pPiM+pPi0).Abs2();

  Complex G3pi = Gomega3pi(spm,sp0,sm0);

  // omega running-width propagator, Eq.(PV-vector-prop)'s D_omega(k^2)
  // with the explicit k^2-dependent Gamma_omega(k^2) - bare 1/D form,
  // NOT the M^2-normalised BreitWigner class (same reasoning as
  // Gomega3pi above).
  Complex Domega_inv(0.,0.);
  if (p_omega_width!=NULL) {
    double Gamma_k2 = (*p_omega_width)(k2);
    Domega_inv = 1./Complex(k2-sqr(m_Momega), m_Momega*Gamma_k2);
  }

  // Vtilde_rho = -W_rho/D_omega(k^2), W_rho = i*G3pi*eps_{rho,b,g,d}
  // p_+^b p_-^g p_0^d = i*G3pi*cross(pPiP,pPiM,pPi0)_rho. Combine the
  // "i" and the "-1/D_omega" into one complex prefactor multiplying
  // the real vector X - see the class-level derivation comment.
  Vec4D  X = cross(pPiP,pPiM,pPi0);
  Complex prefactorV = -Complex(0.,1.)*G3pi*Domega_inv;
  Vec4C  Vtilde = prefactorV * Vec4C(X);

  double Q2_forFF = Q2; // the K-omega system's own Q^2 (production current argument)
  Complex fFF  = FF_f(Q2_forFF);
  Complex gFF  = FF_g(Q2_forFF);
  Complex a1FF = FF_a1(Q2_forFF);
  Complex a2FF = FF_a2(Q2_forFF);

  // First term: -g(Q^2)*cross(p1,p2,Vtilde) = -g(Q^2)*prefactorV*cross(p1,p2,X)
  Vec4C term1 = -gFF*prefactorV*Vec4C(cross(p1,p2,X));
  // Second term: -i*f(Q^2)*Vtilde
  Vec4C term2 = -Complex(0.,1.)*fFF*Vtilde;
  // Third term: -i*{a1*p1+a2*p2}*(Q.Vtilde), Q.Vtilde=prefactorV*(Q*X)
  Complex QdotVtilde = prefactorV*(Q*X);
  Vec4C term3 = -Complex(0.,1.)*(a1FF*Vec4C(p1)+a2FF*Vec4C(p2))*QdotVtilde;

  Insert( m_norm*(term1+term2+term3), 0);
}

void VA_0_KOmega3Pi::SetModelParameters(struct GeneralModel model) {
  // Masses/widths: "may use the values in Table Krho-PV when
  // reproducing the same meson-dominance setup" - per the note, this
  // IS allowed (unlike the couplings below, which must NOT be copied
  // from Krho). Defaults taken from that table.
  m_MK    = model("MK_KOmega",    0.494);
  m_MK1a  = model("MK1_1270_KOmega", 1.272);
  m_GK1a  = model("GK1_1270_KOmega", 0.090);
  m_MK1b  = model("MK1_1400_KOmega", 1.403);
  m_GK1b  = model("GK1_1400_KOmega", 0.174);
  m_MKsta = model("MKstar_892_KOmega",  0.892);
  m_GKsta = model("GKstar_892_KOmega",  0.0508);
  // Second vector is K*(1410), NOT K*(1680). arXiv:0709.4039 Sec.IV is
  // explicit: "g will be assumed to be dominated by the K* = K*(892)
  // and the K'* = K*(1410) intermediate vector mesons (the coupling of
  // the K''* = K*(1680) to the VK^- system is more suppressed)".
  // arXiv:1408.0086 Tab.I labels this column K*(1680) with m=1717,
  // Gamma=322, and fills its weak coupling with 242e3 MeV^2 - which is
  // 0709.4039's f^cqm_{K1'}, the covariant-quark-model coupling of the
  // K1(1400) AXIAL meson, not a vector coupling at all. PDG K*(1410)
  // values are used here instead.
  m_MKstb = model("MKstar_1410_KOmega", 1.414);
  m_GKstb = model("GKstar_1410_KOmega", 0.232);
  m_Momega = Flavour(kf_omega_782).HadMass();
  p_omega_width = LineShapes->Get(Flavour(kf_omega_782));
  p_rho_width   = LineShapes->Get(Flavour(kf_rho_770_plus));


  // h_R*t_{R omega K} coupling products, in GeV.
  //
  // Wang, Guo, Liu and Li, Eur. Phys. J. C74 (2014) 3140
  // [arXiv:1408.0086], Sec.III.A Tab.I, gives this meson-dominance
  // table for K rho.  The same paper states directly, in the paragraph
  // following its Eq.(21), that the LEADING-ORDER production form
  // factors are shared between the two channels: "we adopt the same
  // form factors in the K rho and K omega annihilation processes.
  // According to Eq.(13), we have H_mu^{0 rho} = H_mu^{0 omega}".
  // Their Eq.(13) is exactly the PV current implemented here, so these
  // products come from the source itself and are not an SU(3)
  // assumption of ours.  The relative MINUS sign appearing in the same
  // paragraph, H^{1 rho} = -H^{1 omega}, belongs to the second-order
  // weak diagram Fig.1(b), which exists only to supply a CP-violating
  // phase and has no counterpart in this current.
  //
  // UNITS.  Tab.I heads its rows h_R (10^3 MeV^2) and t_{R V K}
  // (10^-3 MeV^-1) but overrides both inline for the K^- column
  // (0.159 MeV^-1, -3170 MeV).  NEITHER is right.  0709.4039 gives the
  // primary values: f_K = (159.8 +- 1.5) MeV from K->mu nu (dimension
  // MeV), and g_{K+ omega K-} = 3.17 +- 0.03 from phi->KK plus SU(3)
  // with ideal omega-phi mixing (DIMENSIONLESS, its Eq.14).  So the
  // product is 159.8 * -3.17 = -504 MeV = -0.504 GeV, carrying mass^1
  // exactly like the axial and vector products - which is what makes
  // Eq.(14)'s a_1 dimensionally homogeneous, its K-pole and K_1 terms
  // then agreeing at mass^-1.  Reading Tab.I with its row-header units
  // happens to give the same figure through two compensating errors;
  // reading the inline units literally does not, and leaves a_1 short
  // by one power of mass.
  // The resonance masses and widths defaulted above are Tab.I's too.
  m_hKtKwK   = ReadComplexParam(&model,"hKtK_KOmegaMag",  -0.50403,"hKtK_KOmegaPhase");
  m_hA1tA1wK = ReadComplexParam(&model,"hA1tA1_KOmegaMag",-0.41710,"hA1tA1_KOmegaPhase");
  m_hA2tA2wK = ReadComplexParam(&model,"hA2tA2_KOmegaMag", 0.08160,"hA2tA2_KOmegaPhase");
  // Vector sector taken from 0709.4039 directly:
  //   f_K* g_K*wK = 188.9e3 MeV^2 * 8.71e-3 MeV^-1 = 1.6453 GeV  (Eqs.13,23)
  //   f_K'* g_K'*wK = alpha_omegaK * (the above),  alpha_omegaK fitted
  //   by 0709.4039 Sec.V to reproduce B(tau->omega K nu) = 4.1e-4,
  //   giving the two-fold solution alpha = +0.54+-0.38 or -0.77+-0.40
  //   (+0.55 / -0.78 in its cqm variant). The positive branch is used:
  //   0709.4039 Fig.3 notes the phi-K mass distribution favours it.
  //   The tau rate is nearly insensitive to the branch - the two give
  //   widths 0.06% apart - so this choice is about consistency with the
  //   source, not about the normalisation.
  m_hV1tV1wK = ReadComplexParam(&model,"hV1tV1_KOmegaMag", 1.64532,"hV1tV1_KOmegaPhase");
  m_hV2tV2wK = ReadComplexParam(&model,"hV2tV2_KOmegaMag", 0.90493,"hV2tV2_KOmegaPhase");
  // N = g_omegarhopi*g_rhopipi for the omega->3pi vertex, in GeV^-1.
  // Not a free normalisation: a product of two measurable VMD
  // couplings, multiplying the whole current, so the rate goes as
  // |N|^2.  Fixed by requiring the same pure-rho-exchange vertex
  // reproduce Gamma(omega->pi+pi-pi0) = 0.892*8.68 MeV = 7.74 MeV,
  // integrated over the omega Dalitz plot with
  //   sum_pol |M|^2 = |G|^2 * (-X.X),   X^a = eps^{abcd} p+_b p-_c p0_d,
  // and the exact rest-frame identity -X.X = m_omega^2 |p+ x p-|^2
  // (verified numerically to machine precision; the m_omega^2 is easy
  // to drop and costs a factor 1.63 in the rate if you do).
  // That gives N = 179.89 GeV^-1, i.e. g_omegarhopi = 29.98 GeV^-1 for
  // g_rhopipi = 6.0, against the WZW/VMD value 3g^2/(4pi^2 f_pi) =
  // 29.6 GeV^-1 - agreement to 1.3%, an independent check that does
  // not involve the tau rate at all.
  m_Nomega3pi = ReadComplexParam(&model,"Nomega3piMag",179.89,"Nomega3piPhase");

  double Vus = model("Vus", Tools::Vus);
  m_norm = Vus; // K-omega is a |dS|=1 (Cabibbo-suppressed) channel

  msg_Out()<<"### VA_0_KOmega3Pi parameters:\n"
	   <<"###   K: M = "<<m_MK<<" GeV (stable, Gamma=0)\n"
	   <<"###   K1(1270): M = "<<m_MK1a<<" GeV, Gamma = "<<m_GK1a<<" GeV\n"
	   <<"###   K1(1400): M = "<<m_MK1b<<" GeV, Gamma = "<<m_GK1b<<" GeV\n"
	   <<"###   K*(892): M = "<<m_MKsta<<" GeV, Gamma = "<<m_GKsta<<" GeV\n"
	   <<"###   K*(1410): M = "<<m_MKstb<<" GeV, Gamma = "<<m_GKstb<<" GeV\n"
	   <<"###   omega(782): M = "<<m_Momega
	   <<" GeV (running width via LineShapes)\n"
	   <<"###   h_K*t_{K omega K} = "<<m_hKtKwK<<" GeV [1408.0086 Tab.I]\n"
	   <<"###   h_A1*t_{A1 omega K} = "<<m_hA1tA1wK<<" GeV [1408.0086 Tab.I]\n"
	   <<"###   h_A2*t_{A2 omega K} = "<<m_hA2tA2wK<<" GeV [1408.0086 Tab.I]\n"
	   <<"###   h_V1*t_{V1 omega K} = "<<m_hV1tV1wK<<" GeV [1408.0086 Tab.I]\n"
	   <<"###   h_V2*t_{V2 omega K} = "<<m_hV2tV2wK<<" GeV [1408.0086 Tab.I]\n"
	   <<"###   N_omega3pi = "<<m_Nomega3pi<<" GeV^-1 [omega->3pi calibrated]\n";
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_KOmega3Pi,"VA_0_KOmega3Pi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_KOmega3Pi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow K \\omega \\rightarrow K \\pi^+\\pi^-\\pi^0 $ \n\n"
    <<"Order: 0 = K, 1 = pi^+, 2 = pi^-, 3 = pi^0 - all four daughters \n"
    <<"are physically distinct, so this order is required but not a \n"
    <<"special convention choice (no Bose symmetrization needed). \n\n"
    <<"A genuinely NEW four-pseudoscalar current: quasi-two-body K-omega \n"
    <<"weak production (Wang-Guo PV current) with the omega polarization \n"
    <<"eliminated through its own propagator, feeding the anomalous \n"
    <<"omega->3pi vertex (pure rho-exchange VMD realisation). \n\n"
    <<"The h_R*t_{R omega K} strong-coupling products default to the \n"
    <<"meson-dominance values of Wang, Guo, Liu and Li, EPJC 74 (2014) \n"
    <<"3140 [arXiv:1408.0086] Tab.I, which that paper applies to K omega \n"
    <<"as well as K rho (its Eq.(13) H^{0 rho} = H^{0 omega}). Override \n"
    <<"with hKtK_KOmega{Mag,Phase}, hA1tA1_KOmega{Mag,Phase}, \n"
    <<"hA2tA2_KOmega{Mag,Phase}, hV1tV1_KOmega{Mag,Phase}, \n"
    <<"hV2tV2_KOmega{Mag,Phase}. N_omega3pi = g_omegarhopi*g_rhopipi \n"
    <<"is fixed at 140.79 GeV^-1 by calibrating the same omega->3pi \n"
    <<"vertex against the measured Gamma(omega->pi+pi-pi0); it scales \n"
    <<"the whole current, so the rate goes as its square. \n\n"
    <<"Reference: tau_two_meson_currents_KS_RChiT.tex, Sec.'K omega -> \n"
    <<"K pi+pi-pi0: four-pseudoscalar channel' \n"
    <<std::endl;
}
