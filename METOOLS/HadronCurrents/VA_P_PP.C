#include "METOOLS/HadronCurrents/VA_P_PP.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_PP::VA_P_PP(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& indices,
                 const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_fplus(NULL), p_fzero(NULL),
  p_V(NULL), p_A0(NULL), p_A1(NULL), p_A2(NULL),
  p_VT(NULL), p_A0T(NULL), p_A1T(NULL), p_A2T(NULL),
  p_S(NULL), p_P(NULL), p_D(NULL),
  m_cS(1.,0.), m_cP(1.,0.), m_cD(0.,0.),
  m_useS(true), m_useP(true), m_useD(false),
  m_phase(0.5*M_PI), m_m1(0.), m_m2(0.),
  m_isoS(1.), m_isoP(1.), m_isoD(1.),
  m_gS(1.,0.), m_gP(1.,0.), m_gD(1.,0.)
{}

//////////////////////////////////////////////////////////////////////////////
// Wave projectors and resonance-decay vertex factors
//
// In the rest frame of the two-body system L^mu = (0, p_1-p_2), with
// |L| = 2p. The covariant form guaranteeing L.Q = 0 in any frame is
//
//     L^mu = (p_1-p_2)^mu - ((m_1^2-m_2^2)/s) Q^mu,      L.L = -4p^2.
//
// L IS DELIBERATELY NOT NORMALISED. An earlier version divided by 2p
// so that Lhat.Lhat = -1 and Lhat could be dropped in for a
// polarization vector. That reasoning was wrong. The correct
// replacement in the narrow-width limit is
//
//     eps*^mu  ->  g_V BW(s) L^mu,
//
// where g_V is the V -> P P coupling: the decay amplitude is
// g_V eps.(p_1-p_2) = g_V eps.L, so the vertex carries one power of L,
// i.e. one power of the breakup momentum. Normalising L away discarded
// both the coupling and that momentum dependence, and the D -> K pi
// widths came out a factor 14.3 too small, against a phase-space
// weighted <(g 2p)^2> = 12.1 for the K*(892) - i.e. essentially all of
// the deficit. It also distorted the m_Kpi SHAPE, since 2p(s) grows
// with s and its omission suppressed the upper half of the spectrum.
//
// The couplings are derived from each propagator's own pole mass and
// on-shell width rather than being fitted:
//     V -> PP :  Gamma = g^2 p^3 / (6 pi M^2)
//                  ->  g^2 = 6 pi M^2 Gamma / p^3
//     T -> PP :  Gamma = 4 g^2 p^5 / (15 pi M^2)
//                  ->  g^2 = 15 pi M^2 Gamma / (4 p^5)
// using Sum_lambda |eps.L|^2 = -L.L = 4p^2 and
// Sum_lambda |eps_{mu nu}L^mu L^nu|^2 = (2/3)(L.L)^2 = (32/3)p^4.
//
// This also puts the centrifugal p-dependence where it belongs - at
// the production vertex - rather than double-counting it against the
// running width inside the propagator.
//////////////////////////////////////////////////////////////////////////////

Vec4C VA_P_PP::WaveP(const Vec4D & p1,const Vec4D & p2,
                     const Vec4D & Q,const double & s) const
{
  const double dm2 = p1.Abs2()-p2.Abs2();
  return Vec4C((p1-p2) - (dm2/s)*Q);
}

Vec4C VA_P_PP::WaveD(const Vec4D & p1,const Vec4D & p2,
                     const Vec4D & Q,const double & s,
                     const Vec4D & pD) const
{
  // L is real by construction, so work with a Vec4D and convert once.
  const double dm2 = p1.Abs2()-p2.Abs2();
  const Vec4D  L   = (p1-p2) - (dm2/s)*Q;
  // T^{mu nu} = L^mu L^nu - (1/3)(L.L)(g^{mu nu} - Q^mu Q^nu/s),
  // symmetric, traceless and transverse to Q. Contracting one index
  // with p_D and dividing by M gives eps_T of Eq.(effectiveT), so the
  // D-wave term goes through the very same kernel as the vector.
  const double LdotL = L.Abs2();          // = -4 p^2
  const double LdotP = L*pD;
  const double QdotP = Q*pD;
  const Vec4D  T = LdotP*L - (LdotL/3.)*(pD-(QdotP/s)*Q);
  return Vec4C((1./m_M)*T);
}

double VA_P_PP::Breakup(const double & s) const {
  const double n = (s-sqr(m_m1+m_m2))*(s-sqr(m_m1-m_m2));
  return (n>0. ? sqrt(n)/(2.*sqrt(s)) : 0.);
}

Complex VA_P_PP::VertexCoupling(Propagator_Base * prop,const int & L) const {
  // g from the resonance's own pole mass and on-shell width.
  if (prop==NULL) return 1.;
  const double M = prop->Mass(), G = prop->OnShellWidth();
  const double p0 = Breakup(M*M);
  if (M<=0. || G<=0. || p0<=0.) {
    msg_Info()<<"Warning in "<<METHOD<<": cannot derive a "<<(L==1?"V":"T")
              <<" -> PP coupling for "<<m_name<<" (M = "<<M<<", Gamma = "<<G
              <<", p = "<<p0<<"); using 1. The absolute normalisation of "
              <<"this wave is then arbitrary.\n";
    return Complex(1.,0.);
  }
  double g = 1.;
  //     S -> PP :  Gamma = g^2 p / (8 pi M^2)  ->  g^2 = 8 pi M^2 Gamma / p
  if      (L==0) g = sqrt(8.*M_PI*M*M*G/p0);
  else if (L==1) g = sqrt(6.*M_PI*M*M*G/pow(p0,3));
  else if (L==2) g = sqrt(15.*M_PI*M*M*G/(4.*pow(p0,5)));
  // These couplings are defined against the BARE pole 1/(M^2-s-i...),
  // but the propagator objects carry a constant numerator (M^2 for a
  // Breit-Wigner, s_pole for a complex pole) so that they equal 1 at
  // s=0. Divide it out, or the rate is off by |numerator|^2 - which
  // is exactly the M^4 that left D -> K pi at 0.60 of PDG and
  // D0 -> pi- pi0 at 0.349 = m_rho^4.
  const Complex num = prop->PoleNumerator();
  if (std::abs(num)<1.e-12) return Complex(g,0.);
  return Complex(g,0.)/num;
}

//////////////////////////////////////////////////////////////////////////////
// Isospin Clebsch-Gordan
//
// This is group theory, not a tunable, so it is computed from the
// flavours rather than read from the decay table. The quasi-two-body
// table got the 2:1 charge split for free from the K* decay block; a
// current producing the K pi pair directly has to supply it, and when
// it was left to a YAML key it was simply forgotten - both charge
// modes came out with equal rate, ratio 0.98 where isospin and PDG
// both demand 2.0.
//
// For the I=1/2 K pi system reached by c -> s,
//   |1/2,-1/2> = sqrt(2/3)|K0bar pi-> - sqrt(1/3)|K- pi0>     (D^0)
//   |1/2,+1/2> = sqrt(2/3)|K- pi+>   - sqrt(1/3)|K0bar pi0>   (D^+)
// The charged-pion modes get sqrt(2/3), the neutral-pion modes
// sqrt(1/3), with a relative minus sign that is unobservable between
// distinct final states but is kept for correctness.
//
// IMPORTANT for the pi pi and K K currents to come: a SINGLE overall
// factor is only right when every partial wave carries the same
// isospin, which is the case for K pi (all waves are I=1/2) and is
// NOT the case for pi pi, where the S wave (sigma, f_0(980)) is I=0
// and the P wave (rho) is I=1. Those need PER-WAVE factors, which is
// why the overrides below are per wave rather than global.
//////////////////////////////////////////////////////////////////////////////

double VA_P_PP::IsospinCG(const int & wave,bool & known) const
{
  known = true;
  const int par = int(m_flavs[p_i[0]].Kfcode());
  const int kf1 = int(m_flavs[p_i[1]].Kfcode());
  const int kf2 = int(m_flavs[p_i[2]].Kfcode());
  const bool k1 = (kf1==321||kf1==311), k2 = (kf2==321||kf2==311);
  const bool p1 = (kf1==211||kf1==111), p2 = (kf2==211||kf2==111);
  const bool e1 = (kf1==221), e2 = (kf2==221);          // eta

  // ---- K pi : every wave is I=1/2, one coefficient serves all ----
  if ((k1&&p2)||(k2&&p1)) {
    const int kfpi = (p1 ? kf1 : kf2);
    if (kfpi==211) return  sqrt(2./3.);
    return                -sqrt(1./3.);
  }

  // ---- K eta : eta is I=0, so the pair is I=1/2 with a single state
  //      per charge. Both D0 -> K- eta and D+ -> K0bar eta therefore
  //      carry amplitude 1, and their rates differ only by the D
  //      lifetimes. That is the relation used to infer the D+ mode,
  //      which PDG bounds but does not measure. Only even waves
  //      contribute: rho-like K eta P waves are G-parity forbidden and
  //      K*(892) sits below the K eta threshold (0.892 < 1.04 GeV).
  if ((k1&&e2)||(k2&&e1)) {
    if (wave==1) return 0.;
    return 1.;
  }

  // ---- eta pi : eta is I=0, so the pair is pure I=1 ----
  if ((e1&&p2)||(e2&&p1)) {
    if (wave==1) return 0.;   // rho -> eta pi is G-parity forbidden
    if (par==421) return  1.;                 // (d ubar), pure I=1
    if (par==411) return -SQRT_05;            // I=1 part of (d dbar)
    // D_s -> eta pi would need I=1 out of (s sbar): isospin violating,
    // which is why PDG has only an upper limit on it.
    known = false;
    return 1.;
  }

  // ---- pi pi : Bose symmetry ties L to I. The pair must be symmetric
  //      overall, spatial symmetry is (-1)^L, and I=0,2 are symmetric
  //      while I=1 is antisymmetric, so
  //          L even <-> I = 0 or 2,     L odd <-> I = 1.
  //      This is not a modelling choice, it is why a_0(980) - I=1,
  //      J^P=0^+ - cannot decay to pi pi at all and appears in eta pi.
  if (p1&&p2) {
    const bool mixed = (kf1==111)!=(kf2==111);         // pi^+- pi^0
    if (mixed) {
      // I_3 = +-1 excludes I=0 and no I=2 resonance exists, so the S
      // and D waves are FORBIDDEN, not merely unmeasured. PDG agrees
      // independently: B(D0 -> pi- pi0 e nu) = 1.45e-3 against
      // B(D0 -> rho- e nu) = 1.46e-3.
      if (wave==1) return 1.;
      return 0.;
    }
    if (par==431) {
      // c -> s leaves (s sbar), which is I=0. The pi pi pair is then
      // pure I=0: the P wave is forbidden outright (no rho from D_s,
      // and PDG has none), and the RATES split pi+pi- : pi0pi0 = 2 : 1,
      // as the measured f_0(980) products confirm, 1.64e-3/7.9e-4.
      //
      // The AMPLITUDES are EQUAL, not sqrt(2/3) : sqrt(1/3). The
      // remaining factor of two is the identical-particle 1/2! in the
      // pi0 pi0 phase space, which Sherpa already applies - measured
      // directly: with sqrt(2/3) : sqrt(1/3) the generated width ratio
      // came out 3.994 instead of 2, i.e. exactly one factor of two
      // too much, and that factor is the symmetry weight being counted
      // both here and in the phase space. Books distribute the sqrt2
      // between the Clebsch and the symmetry factor differently; the
      // only unambiguous check is the observed 2:1, and that fixes it
      // to equal amplitudes given Sherpa's convention.
      if (wave==1) return 0.;
      return sqrt(2./3.);
    }
    if (kf1==211 && par==411) { if (wave==1) return -SQRT_05; return SQRT_05; }
    if (kf1==111 && par==411) { if (wave==1) return 0.; return SQRT_05; }
  }

  // ---- K Kbar : from D_s the source is (s sbar), I=0, so
  //      |0,0> = (K+K- + K0 K0bar)/sqrt2 and the two charge modes have
  //      EQUAL amplitudes. The observed phi -> K+K- (49.1%) versus
  //      K0 K0bar (34%) split is pure phase space, m_K0 > m_K+, and
  //      comes out of the running width rather than the Clebsch.
  if (k1&&k2) {
    if (par==431) return SQRT_05;
    known = false;
    return 1.;
  }

  known = false;
  return 1.;
}

void VA_P_PP::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  const Vec4D p  = moms[p_i[0]];
  const Vec4D p1 = moms[p_i[1]], p2 = moms[p_i[2]];
  const Vec4D Q  = p1+p2;
  const double s = Q.Abs2();
  if (s<=sqr(m_m1+m_m2)) { Insert(Vec4C(0.,0.,0.,0.),0); return; }

  Vec4C J(0.,0.,0.,0.);

  // S wave: the VA_P_S structure times A_S(s).
  if (m_useS) {
    const Complex ph(cos(m_phase),sin(m_phase));
    J = J + m_isoS*m_cS*m_gS*p_S->Value(s)*ph*ScalarKernel(p,Q,p_fplus,p_fzero);
  }
  // P wave: the VA_P_V structure with eps* -> Lhat, times A_P(s).
  if (m_useP) {
    const Vec4C Lh = WaveP(p1,p2,Q,s);
    J = J + m_isoP*m_cP*m_gP*p_P->Value(s)
          *VectorKernel(p,Q,Lh,p_V,p_A0,p_A1,p_A2,m_epssign);
  }
  // D wave: the same kernel with the spin-2 projector, times A_D(s).
  if (m_useD) {
    const Vec4C Td = WaveD(p1,p2,Q,s,p);
    J = J + m_isoD*m_cD*m_gD*p_D->Value(s)
          *VectorKernel(p,Q,Td,p_VT,p_A0T,p_A1T,p_A2T,m_epssign);
  }

  J = m_norm*J;
  Insert(anti?conj(J):J,0);
}

FF_PP_Wave * VA_P_PP::MakeWave(const std::string & name,GeneralModel & model)
{
  std::vector<int> pp; pp.push_back(p_i[1]); pp.push_back(p_i[2]);
  std::map<std::string,double> pmap;
  FF_Parameters params(ff_model::none,m_flavs,pp,pmap,name,&model);
  FormFactor_Base * ff = FF_Getter::GetObject("FF_PP_Wave",params);
  if (ff==NULL)
    THROW(fatal_error,"Could not build wave lineshape '"+name+"' for current "
          +m_name+". FF_Getter returned NULL for tag 'FF_PP_Wave' - it takes "
          "exactly the two pseudoscalar indices, check FF_PP_Wave.C is built "
          "and registered.");
  FF_PP_Wave * w = dynamic_cast<FF_PP_Wave*>(ff);
  if (w==NULL) THROW(fatal_error,"'"+name+"' is not an FF_PP_Wave.");
  return w;
}

void VA_P_PP::SetModelParameters(struct GeneralModel model)
{
  CheckArity(2);
  ReadCommonParameters(model);
  m_phase = model("SCALAR_PHASE",0.5*M_PI);
  m_m1 = m_flavs[p_i[1]].HadMass();
  m_m2 = m_flavs[p_i[2]].HadMass();

  // With a genuinely running two-body mass there is no nominal
  // resonance mass to fall back on, so every form-factor shape must
  // follow sqrt(s) or t_+ becomes inconsistent with Delta = M^2 - s.
  const char * ffs[] = {"Fplus","Fzero","V","A0","A1","A2",
                        "VT","A0T","A1T","A2T"};
  for (size_t i=0;i<10;i++) model[std::string(ffs[i])+"_DYNAMIC_MASS"] = 1.;

  // Per-wave isospin factors: computed by default, overridable per
  // wave for systems whose waves do not share one isospin.
  // Per WAVE, because pi pi does not share one isospin across waves:
  // the S wave is I=0 (sigma, f_0), the P wave I=1 (rho).
  bool kS=false, kP=false, kD=false;
  m_isoS = model("ISO_S",IsospinCG(0,kS));
  m_isoP = model("ISO_P",IsospinCG(1,kP));
  m_isoD = model("ISO_D",IsospinCG(2,kD));
  const bool known = kS&&kP&&kD;
  // Never fall back on 1 silently: that is exactly how the K pi 2:1
  // charge split was lost the first time round.
  if (!known && model("ISO_S",-1.e30)<-1.e29)
    THROW(fatal_error,"Current "+m_name+" cannot derive isospin "
          "Clebsch-Gordan coefficients for the pair "
          +ToString(m_flavs[p_i[1]])+" "+ToString(m_flavs[p_i[2]])+
          " from a "+ToString(m_flavs[p_i[0]])+" parent, and none were "
          "supplied. Set ISO_S / ISO_P / ISO_D explicitly.");
  if (model("ISOSPIN",-1.)>0.)
    THROW(fatal_error,"ISOSPIN was set on current "+m_name+
          ". Use the per-wave keys ISO_S / ISO_P / ISO_D instead: a single "
          "global factor is wrong whenever the partial waves carry "
          "different isospin, and the K pi coefficients are derived "
          "automatically anyway.");

  m_useS = (model("USE_S",1.)>0.5);
  m_useP = (model("USE_P",1.)>0.5);
  m_useD = (model("USE_D",0.)>0.5);

  // Complex wave coefficients. A relative PHASE between waves is not
  // optional: it is precisely what the interference measures, and
  // leaving every phase at zero silently asserts that all waves add in
  // phase everywhere.
  const double aS = model("AMP_S",1.), phS = model("PHASE_S",0.);
  const double aP = model("AMP_P",1.), phP = model("PHASE_P",0.);
  const double aD = model("AMP_D",1.), phD = model("PHASE_D",0.);
  m_cS = Complex(aS*cos(phS),aS*sin(phS));
  m_cP = Complex(aP*cos(phP),aP*sin(phP));
  m_cD = Complex(aD*cos(phD),aD*sin(phD));

  if (m_useS) {
    p_fplus = MakeFF("Fplus",model);
    p_fzero = MakeFF("Fzero",model);
    const double mref = model("TIE_MASS",1.0);
    if (model("TIE_FZERO",1.)>0.5) {
      const double fp0 = p_fplus->Value(0.,mref);
      const double f00 = p_fzero->Value(0.,mref);
      if (std::abs(f00)>1.e-300)
        p_fzero->SetF0(p_fzero->GetF0()*fp0/f00);
    }
    p_S = MakeWave("Swave",model);
    if (p_S->L()!=0) THROW(fatal_error,"Swave_L must be 0 in "+m_name+".");
    // A LASS amplitude is already unitary-normalised, so it must NOT be
    // multiplied by a decay coupling; a bare resonance must.
    m_gS = (p_S->Model()==pp_wave_model::LASS
            ? Complex(1.,0.) : VertexCoupling(p_S->Prop(),0));
  }
  if (m_useP) {
    p_V  = MakeFF("V", model); p_A0 = MakeFF("A0",model);
    p_A1 = MakeFF("A1",model); p_A2 = MakeFF("A2",model);
    p_P  = MakeWave("Pwave",model);
    if (p_P->L()!=1) THROW(fatal_error,"Pwave_L must be 1 in "+m_name+".");
    m_gP = VertexCoupling(p_P->Prop(),1);
    // Same {A_1(0), r_V, r_2} input route as VA_P_V, evaluated at the
    // resonance pole mass, so published on-shell numbers transfer.
    if (model("USE_RATIOS",1.)>0.5) {
      const double A1_0 = model("A1_0",1.), rV = model("rV",1.),
                   r2   = model("r2",1.), mref = model("PWAVE_MASS",0.89555);
      const double v0=p_V->Value(0.,mref), a10=p_A1->Value(0.,mref),
                   a20=p_A2->Value(0.,mref);
      if (std::abs(a10)>1.e-300) p_A1->SetF0(p_A1->GetF0()*A1_0/a10);
      if (std::abs(v0) >1.e-300) p_V ->SetF0(p_V ->GetF0()*rV*A1_0/v0);
      if (std::abs(a20)>1.e-300) p_A2->SetF0(p_A2->GetF0()*r2*A1_0/a20);
    }
  }
  if (m_useD) {
    p_VT  = MakeFF("VT", model); p_A0T = MakeFF("A0T",model);
    p_A1T = MakeFF("A1T",model); p_A2T = MakeFF("A2T",model);
    p_D   = MakeWave("Dwave",model);
    if (p_D->L()!=2) THROW(fatal_error,"Dwave_L must be 2 in "+m_name+".");
    m_gD = VertexCoupling(p_D->Prop(),2);
    if (model("USE_RATIOS",1.)>0.5) {
      const double A1_0 = model("A1T_0",1.), rV = model("rVT",1.),
                   r2   = model("r2T",1.), mref = model("DWAVE_MASS",1.4324);
      const double v0=p_VT->Value(0.,mref), a10=p_A1T->Value(0.,mref),
                   a20=p_A2T->Value(0.,mref);
      if (std::abs(a10)>1.e-300) p_A1T->SetF0(p_A1T->GetF0()*A1_0/a10);
      if (std::abs(v0) >1.e-300) p_VT ->SetF0(p_VT ->GetF0()*rV*A1_0/v0);
      if (std::abs(a20)>1.e-300) p_A2T->SetF0(p_A2T->GetF0()*r2*A1_0/a20);
    }
  }

  PrintCommonParameters();
  msg_Tracking()<<"###   coherent two-body current, waves:"
                <<(m_useS?" S":"")<<(m_useP?" P":"")<<(m_useD?" D":"")
                <<", threshold sqrt(s) = "<<(m_m1+m_m2)<<" GeV\n"
                <<"###   vertex couplings |g_S| = "<<std::abs(m_gS)
                <<", |g_P| = "<<std::abs(m_gP)
                <<", |g_D| = "<<std::abs(m_gD)<<"\n"
                <<"###   |c_S| = "<<std::abs(m_cS)<<", |c_P| = "
                <<std::abs(m_cP)<<", |c_D| = "<<std::abs(m_cD)<<"\n"
                <<"###   NOTE: this covers the FULL two-body final state. "
                <<"Do not also enable quasi-two-body channels for any "
                <<"resonance included here.\n";
}

DEFINE_CURRENT_GETTER(METOOLS::VA_P_PP,"VA_P_PP")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_PP>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $D^+\\rightarrow K^-\\pi^+ e^+\\nu$, all waves coherent \n\n"
    <<"Order: 0 = decaying $D$, 1 and 2 = the two pseudoscalars \n\n"
    <<"$H^\\mu=\\sum_L c_L A_L(s) K_L^\\mu$, with $K_0$ the $D\\to S$ \n"
    <<"kernel, $K_1$ the $D\\to V$ kernel with $\\eps^*\\to\\hat L$ and \n"
    <<"$K_2$ the same kernel with the spin-2 projector built from \n"
    <<"$\\hat L$. Enable waves with {\\tt USE\\_S/P/D}; set relative \n"
    <<"strengths with {\\tt AMP\\_*} and {\\tt PHASE\\_*}. \n\n"
    <<"Lineshapes {\\tt Swave}, {\\tt Pwave}, {\\tt Dwave} are \n"
    <<"{\\tt FF\\_PP\\_Wave} objects; their $L$ must match the wave. \n\n"
    <<"Covers the WHOLE two-body final state - do not combine with \n"
    <<"quasi-two-body channels for the resonances included here. \n"
    <<std::endl;
}
