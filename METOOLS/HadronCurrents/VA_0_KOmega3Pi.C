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
  // Eq.(PV-DR): D_R(Q^2) = Q^2-M_R^2+iM_R Gamma_R (FIXED width, as
  // literally written - no (Q^2) dependence on Gamma_R shown in this
  // equation, unlike the omega propagator below). Returns 1/D_R.
  return 1./Complex(Q2-sqr(M), M*G);
}

Complex VA_0_KOmega3Pi::FF_f(const double & Q2) const {
  // Eq.(PV-f): f(Q^2) = -1/2(Q^2+M_omega^2-m_K^2) sum_i hA_i tA_i/D_Ai(Q^2)
  Complex SA = m_hA1tA1wK*DRfixed(Q2,m_MK1a,m_GK1a)
             + m_hA2tA2wK*DRfixed(Q2,m_MK1b,m_GK1b);
  return -0.5*(Q2+sqr(m_Momega)-sqr(m_MK))*SA;
}

Complex VA_0_KOmega3Pi::FF_g(const double & Q2) const {
  // Eq.(PV-g): g(Q^2) = 1/2 sum_i hV_i tV_i/D_Vi(Q^2)
  Complex SV = m_hV1tV1wK*DRfixed(Q2,m_MKsta,m_GKsta)
             + m_hV2tV2wK*DRfixed(Q2,m_MKstb,m_GKstb);
  return 0.5*SV;
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
  m_MKstb = model("MKstar_1680_KOmega", 1.717);
  m_GKstb = model("GKstar_1680_KOmega", 0.322);
  m_Momega = Flavour(kf_omega_782).HadMass();
  p_omega_width = LineShapes->Get(Flavour(kf_omega_782));
  p_rho_width   = LineShapes->Get(Flavour(kf_rho_770_plus));

  // FIXME: h_R*t_{R omega K} coupling products. The note is explicit
  // ("this catalogue does not silently substitute the K rho numbers
  // ... A code should expose h_R t_{R omega K} as an explicit
  // parameter block") that these must NOT default to the Krho table's
  // own h_R*t_{R rho K} values - defaulting to 0 here instead, which
  // makes this current IDENTICALLY ZERO until real Komega numbers are
  // supplied. This is a deliberate, conservative choice (matching the
  // note's own caution) rather than a silent wrong-physics default.
  m_hKtKwK   = ReadComplexParam(&model,"hKtK_KOmegaMag",0.,"hKtK_KOmegaPhase");
  m_hA1tA1wK = ReadComplexParam(&model,"hA1tA1_KOmegaMag",0.,"hA1tA1_KOmegaPhase");
  m_hA2tA2wK = ReadComplexParam(&model,"hA2tA2_KOmegaMag",0.,"hA2tA2_KOmegaPhase");
  m_hV1tV1wK = ReadComplexParam(&model,"hV1tV1_KOmegaMag",0.,"hV1tV1_KOmegaPhase");
  m_hV2tV2wK = ReadComplexParam(&model,"hV2tV2_KOmegaMag",0.,"hV2tV2_KOmegaPhase");
  // g_omegarhopi*g_rhopipi normalisation for the omega->3pi vertex -
  // also not given a specific number in the note; placeholder default
  // 1, like N_Keta/N_etapipi elsewhere in this codebase.
  m_Nomega3pi = ReadComplexParam(&model,"Nomega3piMag",1.,"Nomega3piPhase");

  double Vus = model("Vus", Tools::Vus);
  m_norm = Vus; // K-omega is a |dS|=1 (Cabibbo-suppressed) channel

  msg_Out()<<"### VA_0_KOmega3Pi parameters:\n"
	   <<"###   K: M = "<<m_MK<<" GeV (stable, Gamma=0)\n"
	   <<"###   K1(1270): M = "<<m_MK1a<<" GeV, Gamma = "<<m_GK1a<<" GeV\n"
	   <<"###   K1(1400): M = "<<m_MK1b<<" GeV, Gamma = "<<m_GK1b<<" GeV\n"
	   <<"###   K*(892): M = "<<m_MKsta<<" GeV, Gamma = "<<m_GKsta<<" GeV\n"
	   <<"###   K*(1680): M = "<<m_MKstb<<" GeV, Gamma = "<<m_GKstb<<" GeV\n"
	   <<"###   omega(782): M = "<<m_Momega
	   <<" GeV (running width via LineShapes)\n"
	   <<"###   h_K*t_{K omega K} = "<<m_hKtKwK<<" (0 = FIXME, not yet supplied)\n"
	   <<"###   h_A1*t_{A1 omega K} = "<<m_hA1tA1wK<<" (0 = FIXME)\n"
	   <<"###   h_A2*t_{A2 omega K} = "<<m_hA2tA2wK<<" (0 = FIXME)\n"
	   <<"###   h_V1*t_{V1 omega K} = "<<m_hV1tV1wK<<" (0 = FIXME)\n"
	   <<"###   h_V2*t_{V2 omega K} = "<<m_hV2tV2wK<<" (0 = FIXME)\n"
	   <<"###   N_omega3pi = "<<m_Nomega3pi<<" (placeholder, default 1)\n";
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
    <<"Status: the h_R*t_{R omega K} strong-coupling products (5 complex \n"
    <<"parameters: hKtK_KOmega{Mag,Phase}, hA1tA1_KOmega{Mag,Phase}, \n"
    <<"hA2tA2_KOmega{Mag,Phase}, hV1tV1_KOmega{Mag,Phase}, \n"
    <<"hV2tV2_KOmega{Mag,Phase}) all default to 0 - this current is \n"
    <<"IDENTICALLY ZERO until real Komega coupling values are supplied; \n"
    <<"the source note is explicit that these must NOT be copied from \n"
    <<"the analogous Krho values (SU(3) is a model input, not an \n"
    <<"identity). N_omega3pi (default 1) is a placeholder normalisation \n"
    <<"for the omega->3pi vertex (g_omegarhopi*g_rhopipi), also not \n"
    <<"given a specific number in the source. \n\n"
    <<"Reference: tau_two_meson_currents_KS_RChiT.tex, Sec.'K omega -> \n"
    <<"K pi+pi-pi0: four-pseudoscalar channel' \n"
    <<std::endl;
}
