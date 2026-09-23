#include "METOOLS/HadronCurrents/FormFactors/FF_P_X.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/MathTools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

//////////////////////////////////////////////////////////////////////////////
//
// All formulae below are taken directly from charm_semileptonic.tex,
// Sec."Generic form-factor parametrizations" and Sec."z-expansion
// parametrizations". Equation labels in the comments refer to that
// document.
//
// Parameter keys are prefixed with the form factor's own name
// (m_name, set by the owning current: "Fplus", "Fzero", "V", "A0",
// "A1", "A2", "A", "V0", "V1", "V2", "VT", "A0T", "A1T", "A2T"), so a
// single decay channel can carry a completely independent shape and
// parameter set for each invariant form factor without any key
// collisions - the failure mode that had to be repaired by hand in
// the tau sector (see the "Parameter naming" block in Decaydata.yaml).
//
//////////////////////////////////////////////////////////////////////////////

FF_P_X::FF_P_X(const FF_Parameters & params) :
  FormFactor_Base(params),
  m_shape(ffq2_shape::constant),
  m_F0(1.),
  m_mpole(0.), m_mpole2(0.),
  m_alpha(0.), m_beta(1.), m_r(0.),
  m_c1(0.), m_c2(0.),
  m_a(0.), m_b(0.), m_Lambda2(1.),
  m_rISGW2(0.),
  m_superconv(false),
  m_mc(1.3), m_m0plus(0.), m_useBlaschke(false),
  m_t0user(0.), m_hast0user(false),
  m_M(0.), m_m(0.),
  m_dynamicmass(false),
  m_q2min(1.e-8)
{
  // m_pi[0] = parent (D), m_pi[1] = daughter hadron X.
  m_M = m_masses[m_pi[0]];
  m_m = m_masses[m_pi[1]];
  ReadParameters();
}

void FF_P_X::ReadParameters() {
  if (p_model==NULL) return;
  const GeneralModel & md = *p_model;
  const string & n = m_name;

  m_shape = ffq2_shape(int(md(n+"_SHAPE",double(int(ffq2_shape::constant)))));

  m_F0      = md(n+"_F0",      1.);
  m_mpole   = md(n+"_MPOLE",   0.);
  m_mpole2  = m_mpole*m_mpole;
  m_alpha   = md(n+"_ALPHA",   0.);
  m_beta    = md(n+"_BETA",    1.);
  m_r       = md(n+"_R",       0.);
  m_c1      = md(n+"_C1",      0.);
  m_c2      = md(n+"_C2",      0.);
  m_a       = md(n+"_A",       0.);
  m_b       = md(n+"_B",       0.);
  double L  = md(n+"_LAMBDA",  1.);
  m_Lambda2 = L*L;
  m_rISGW2  = md(n+"_RISGW2",  0.);
  m_mc      = md(n+"_MC",      1.3);
  m_m0plus  = md(n+"_M0PLUS",  0.);
  m_dynamicmass = (md(n+"_DYNAMIC_MASS",0.)>0.5);
  m_q2min   = md(n+"_Q2MIN",   1.e-8);

  // --- multipole / residue expansion ---
  int npole = int(md(n+"_NPOLE",0.));
  for (int i=0;i<npole;i++) {
    m_mn.push_back(md(n+"_m"+ToString(i),0.));
    m_Rn.push_back(md(n+"_R"+ToString(i),0.));
  }
  m_superconv = (md(n+"_SUPERCONVERGENT",0.)>0.5);
  if (m_superconv && m_Rn.size()>0) {
    // Impose R_0+...+R_{N-1}=0 by fixing the LAST residue, rather than
    // trusting the YAML to satisfy it numerically (note item 2).
    double sum = 0.;
    for (size_t i=0;i+1<m_Rn.size();i++) sum += m_Rn[i];
    m_Rn.back() = -sum;
  }

  // --- z-expansion coefficient block. Which object m_an holds is
  //     shape-dependent, by design (they are genuinely different
  //     quantities and must not be silently interchanged - see the
  //     warning under the experimental table in the note):
  //       BGL          : a_0 ... a_N            (key _a0.._aN)
  //       charm_series : r_1 ... r_N = a_n/a_0  (key _r1.._rN)
  //       BCL_plus/zero: b_0 ... b_{N-1}        (key _a0.._a(N-1))
  //       SSE          : a_1 ... a_N, a_0=F(0)  (key _a1.._aN)
  int nc = int(md(n+"_NCOEFF",0.));
  for (int i=0;i<nc;i++) {
    switch (m_shape) {
    case ffq2_shape::charm_series:
      m_an.push_back(md(n+"_r"+ToString(i+1),0.)); break;
    case ffq2_shape::SSE:
      m_an.push_back(md(n+"_a"+ToString(i+1),0.)); break;
    default:
      m_an.push_back(md(n+"_a"+ToString(i),0.));   break;
    }
  }
  int nb = int(md(n+"_NBLASCHKE",0.));
  for (int i=0;i<nb;i++) {
    const double mR = md(n+"_mR"+ToString(i),0.);
    if (mR<=0.)
      THROW(fatal_error,"Blaschke pole "+n+"_mR"+ToString(i)+" of form factor '"
            +n+"' is missing or non-positive. A zero pole mass would make "
            "P_F(0)=0 and the form factor identically zero.");
    m_blaschke.push_back(mR);
  }
  m_useBlaschke = (nb>0);

  m_hast0user = false;
  const double t0in = md(n+"_T0",-1.e30);
  if (t0in>-1.e29) { m_t0user = t0in; m_hast0user = true; }

  // A z-expansion with NO Blaschke block is a legitimate choice (it is
  // the correct one for D->pi, where the D* sits above t_+), but it is
  // also exactly what a mistyped or dropped key looks like, and the two
  // are indistinguishable at run time from the shape alone. Since the
  // consequence is a form factor that rises far too slowly - roughly a
  // 25% deficit in f_+(q^2_max) for D->K, i.e. >40% in the rate at the
  // top of the spectrum - say so out loud rather than proceeding
  // silently. Nothing here changes the physics; it only makes the
  // choice visible in the log.
  if ((m_shape==ffq2_shape::BGL || m_shape==ffq2_shape::charm_series ||
       m_shape==ffq2_shape::BCL_plus) && !m_useBlaschke)
    msg_Info()<<"Warning in "<<METHOD<<": form factor '"<<n<<"' uses a "
              <<"z-expansion with NO Blaschke factor ("<<n<<"_NBLASCHKE "
              <<"unset or 0). Correct for D->pi; a silent key typo "
              <<"otherwise. Check this is intended.\n";

  // p_model points at a by-value GeneralModel owned by the caller's
  // SetModelParameters frame and is dangling from here on. Null it so
  // that any future per-event use segfaults immediately instead of
  // quietly reading freed memory.
  p_model = NULL;
  PrintSetup();
}

void FF_P_X::PrintSetup() const {
  // One-shot summary so that the ACTUAL configuration - not the
  // intended one - is visible in the log, including which subthreshold
  // poles survived the m_R^2 < t_+ test and what the resulting form
  // factor looks like at the two ends of the physical range.
  const double tp = sqr(m_M+m_m), q2max = sqr(m_M-m_m);
  msg_Tracking()<<"### FF_P_X '"<<m_name<<"': shape "<<int(m_shape)
                <<", F(0) stored = "<<m_F0
                <<", M = "<<m_M<<", m = "<<m_m<<", t_+ = "<<tp<<"\n";
  if (m_an.size())
    msg_Tracking()<<"###   "<<m_an.size()<<" z-expansion coefficient(s), "
                  <<"first = "<<m_an[0]<<"\n";
  if (m_useBlaschke) {
    for (size_t j=0;j<m_blaschke.size();j++) {
      const bool sub = (sqr(m_blaschke[j])<tp);
      msg_Tracking()<<"###   Blaschke pole m_R = "<<m_blaschke[j]<<" GeV -> "
                    <<(sub?"KEPT (subthreshold)":"DROPPED (m_R^2 >= t_+)")<<"\n";
    }
  }
  else msg_Tracking()<<"###   no Blaschke factor\n";
  msg_Tracking()<<"###   F(0) = "<<Value(0.,m_m)
                <<", F(q^2max="<<q2max<<") = "<<Value(q2max,m_m)
                <<", ratio = "<<(Value(0.,m_m)!=0.?Value(q2max,m_m)/Value(0.,m_m):0.)
                <<"\n";
}

//////////////////////////////////////////////////////////////////////////////
// z-expansion machinery
//////////////////////////////////////////////////////////////////////////////

void FF_P_X::SetupZ(const double & m,double & tp,double & tm,
                    double & t0) const {
  tp = sqr(m_M+m);
  tm = sqr(m_M-m);
  // Eq.(t0opt): the standard |z|_max-minimising choice.
  t0 = tp*(1.-sqrt(1.-tm/tp));
  // The user override is CACHED at construction. It must not be read
  // from p_model here: Current_Base::SetModelParameters takes its
  // GeneralModel BY VALUE, so the map this pointer refers to is
  // destroyed the moment the owning current finishes its setup, and
  // any per-event dereference is undefined behaviour. p_model is
  // therefore nulled at the end of ReadParameters().
  if (m_hast0user) t0 = m_t0user;
}

double FF_P_X::Zvar(const double & q2,const double & t0,
                    const double & tp) const {
  // Eq.(zdef).
  double A = sqrt(Max(0.,tp-q2)), B = sqrt(Max(0.,tp-t0));
  if (A+B<1.e-12) return 0.;
  return (A-B)/(A+B);
}

double FF_P_X::ZOverDelta(const double & q2,const double & a,
                          const double & tp) const {
  // z(q^2,a)/(a-q^2). Both numerator and denominator vanish at q^2=a;
  // the analytic limit is 1/(4(t_+-a)), which is exactly what the note
  // means by "evaluated by their analytic limits when numerator and
  // denominator vanish simultaneously" (text under Eq.(charmphi)).
  double d = a-q2;
  if (std::abs(d)<1.e-9*Max(1.,std::abs(tp))) return 1./(4.*(tp-a));
  return Zvar(q2,a,tp)/d;
}

double FF_P_X::OuterPhiPlus(const double & q2,const double & t0,
                            const double & tp,const double & tm) const {
  // Eq.(charmphi). NOTE: in the F(0)-normalised charm series below,
  // every q^2-INDEPENDENT prefactor of phi (including sqrt(pi m_c^2/3)
  // and (t_+-t_0)^{-1/4}) cancels identically between numerator and
  // denominator. They are kept here anyway so that this function is
  // also usable for a genuine BGL fit, where the absolute
  // normalisation of the a_n does matter.
  double pref = sqrt(M_PI*m_mc*m_mc/3.);
  double f1   = pow(ZOverDelta(q2,0. ,tp), 2.5 );
  double f2   = pow(ZOverDelta(q2,t0 ,tp),-0.5 );
  double f3   = pow(ZOverDelta(q2,tm ,tp),-0.75);
  double f4   = (tp-q2)/pow(tp-t0,0.25);
  return pref*f1*f2*f3*f4;
}

double FF_P_X::BlaschkeP(const double & q2,const double & tp) const {
  // P_F(q^2) = prod_j z(q^2, m_Rj^2), over SUBTHRESHOLD poles only.
  //
  // The note's own criterion is that "a pole factor should only be
  // treated as a true Blaschke/subthreshold pole when its pole
  // location lies below the appropriate crossed-channel branch point",
  // i.e. iff m_R^2 < t_+. That test is applied here automatically
  // rather than being left to the YAML, because it is exactly what
  // distinguishes the two prescriptions the note quotes explicitly:
  //   D -> K   : m_{D_s^*}^2 = 4.46 < t_+ = 5.56  -> pole kept
  //   D -> pi  : m_{D^*}^2   = 4.04 > t_+ = 4.02  -> P_+ = 1
  // and it then extends the same rule, unprompted and correctly, to
  // the channels the note does not spell out. Numerically: with the
  // D^* pole switched OFF, the D_s -> K^0 charm series gives
  // f_+(q^2_max) = 1.19 against 2.03 from the same analysis' own
  // modified-pole fit; with it switched ON (m_{D^*}^2 = 4.04 <
  // t_+ = 6.08) the series gives 1.97, i.e. the two parametrizations
  // of that measurement agree. The same switch reproduces the
  // published D_s -> eta / eta' alternative fits to ~0.5%.
  //
  // A pole sitting ABOVE threshold is silently skipped (it is not a
  // Blaschke factor at all there, and z(q^2,m_R^2) would not even be
  // real), so one uniform "lowest c-qbar vector pole" entry can be
  // given for every channel.
  if (!m_useBlaschke) return 1.;
  double P = 1.;
  for (size_t j=0;j<m_blaschke.size();j++) {
    const double mR2 = sqr(m_blaschke[j]);
    if (mR2>=tp) continue;
    P *= Zvar(q2,mR2,tp);
  }
  return P;
}

//////////////////////////////////////////////////////////////////////////////
// The form factor itself
//////////////////////////////////////////////////////////////////////////////

double FF_P_X::Value(const double & q2,const double & mX) const {
  const double m = (m_dynamicmass ? mX : m_m);
  double tp,tm,t0;

  switch (m_shape) {

  case ffq2_shape::constant:
    return m_F0;

  case ffq2_shape::taylor: {
    double x = q2/m_Lambda2;
    return m_F0*(1.+m_c1*x+m_c2*x*x);
  }

  case ffq2_shape::simple_pole: {
    if (m_mpole2<=0.) return m_F0;
    return m_F0/(1.-q2/m_mpole2);
  }

  case ffq2_shape::modified_pole: {
    if (m_mpole2<=0.) return m_F0;
    double x = q2/m_mpole2;
    return m_F0/((1.-x)*(1.-m_alpha*x));
  }

  case ffq2_shape::BK_plus: {
    // f+ = F/[(1-x)(1-alpha x)], x=q^2/m_{1^-}^2. Identical in form to
    // modified_pole but kept as its own enum value because BK ties f+
    // and f0 to ONE normalisation F and one pole mass - see BK_zero.
    if (m_mpole2<=0.) return m_F0;
    double x = q2/m_mpole2;
    return m_F0/((1.-x)*(1.-m_alpha*x));
  }

  case ffq2_shape::BK_zero: {
    // f0 = F/(1-x/beta). The BK constraint f+(0)=f0(0)=F is automatic.
    if (m_mpole2<=0.) return m_F0;
    double x = q2/m_mpole2;
    if (std::abs(m_beta)<1.e-12) return m_F0;
    return m_F0/(1.-x/m_beta);
  }

  case ffq2_shape::BZ: {
    if (m_mpole2<=0.) return m_F0;
    double x = q2/m_mpole2;
    return m_F0*(1./(1.-x) + m_r*x/((1.-x)*(1.-m_alpha*x)));
  }

  case ffq2_shape::rational: {
    double x = q2/m_Lambda2;
    return m_F0/(1.-m_a*x+m_b*x*x);
  }

  case ffq2_shape::exponential: {
    double x = q2/m_Lambda2;
    return m_F0*exp(m_a*x+m_b*x*x);
  }

  case ffq2_shape::multipole: {
    double F = 0.;
    for (size_t i=0;i<m_Rn.size();i++) {
      double den = sqr(m_mn[i])-q2;
      if (std::abs(den)<1.e-12) continue;
      F += m_Rn[i]/den;
    }
    return m_F0*F;
  }

  case ffq2_shape::ISGW2: {
    // f+(q^2)=f+(qmax^2)[1+r^2/12 (qmax^2-q^2)]^{-2}, qmax^2=(M-m)^2.
    double qmax2 = sqr(m_M-m);
    double br    = 1.+m_rISGW2*m_rISGW2/12.*(qmax2-q2);
    if (std::abs(br)<1.e-12) return 0.;
    return m_F0/(br*br);
  }

  case ffq2_shape::BGL: {
    SetupZ(m,tp,tm,t0);
    double z = Zvar(q2,t0,tp), zn = 1., sum = 0.;
    for (size_t i=0;i<m_an.size();i++) { sum += m_an[i]*zn; zn *= z; }
    double PF = BlaschkeP(q2,tp), ph = OuterPhiPlus(q2,t0,tp,tm);
    if (std::abs(PF*ph)<1.e-30) return 0.;
    return m_F0*sum/(PF*ph);
  }

  case ffq2_shape::charm_series: {
    // F(q^2) = F(0) P(0)phi(0)[1+sum r_n z(q^2)^n]
    //          / { P(q^2)phi(q^2)[1+sum r_n z(0)^n] }.
    // This is the parametrization the BESIII "2-parameter series"
    // fits in the note's experimental tables actually use, with the
    // single shape parameter r_1. It is NOT a BCL b_1 and not the
    // a_1 of Eq.(SSE) - do not transplant values between them.
    SetupZ(m,tp,tm,t0);
    double z   = Zvar(q2,t0,tp);
    double z00 = Zvar(0.,t0,tp);
    double num = 1., den = 1., zn = z, zn0 = z00;
    for (size_t i=0;i<m_an.size();i++) {
      num += m_an[i]*zn;  zn  *= z;
      den += m_an[i]*zn0; zn0 *= z00;
    }
    double PFq = BlaschkeP(q2,tp), phq = OuterPhiPlus(q2,t0,tp,tm);
    double PF0 = BlaschkeP(0.,tp), ph0 = OuterPhiPlus(0.,t0,tp,tm);
    if (std::abs(PFq*phq*den)<1.e-30) return 0.;
    return m_F0*PF0*ph0*num/(PFq*phq*den);
  }

  case ffq2_shape::BCL_plus: {
    // Eq.(BCLplus). m_an = b_0..b_{N-1}; the n/N z^N subtraction
    // enforces the correct threshold behaviour.
    SetupZ(m,tp,tm,t0);
    int N = int(m_an.size());
    if (N==0) return m_F0;
    double z = Zvar(q2,t0,tp);
    double zN = pow(z,double(N)), zn = 1., sum = 0.;
    for (int n=0;n<N;n++) {
      double sign = ((n-N)%2==0 ? 1. : -1.);   // (-1)^{n-N}
      sum += m_an[n]*(zn - sign*double(n)/double(N)*zN);
      zn *= z;
    }
    double pole = (m_mpole2>0. ? 1.-q2/m_mpole2 : 1.);
    if (std::abs(pole)<1.e-12) return 0.;
    return m_F0*sum/pole;
  }

  case ffq2_shape::BCL_zero: {
    // Eq.(BCLzero). B_0=1, or 1-q^2/m_{0^+}^2 if a scalar pole is
    // factored out explicitly.
    SetupZ(m,tp,tm,t0);
    double z = Zvar(q2,t0,tp), zn = 1., sum = 0.;
    for (size_t n=0;n<m_an.size();n++) { sum += m_an[n]*zn; zn *= z; }
    double B0 = (m_m0plus>0. ? 1.-q2/sqr(m_m0plus) : 1.);
    if (std::abs(B0)<1.e-12) return 0.;
    return m_F0*sum/B0;
  }

  case ffq2_shape::SSE: {
    // Eq.(SSE): F=sum a_n [z(q^2,t0)-z(0,t0)]^n /(1-q^2/m_R^2), with
    // a_0=F(0) by construction of the shifted variable. This is the
    // note's own recommended "cleanest implementation-ready choice".
    SetupZ(m,tp,tm,t0);
    double dz = Zvar(q2,t0,tp)-Zvar(0.,t0,tp);
    double sum = m_F0, dzn = dz;
    for (size_t i=0;i<m_an.size();i++) { sum += m_an[i]*dzn; dzn *= dz; }
    double pole = (m_mpole2>0. ? 1.-q2/m_mpole2 : 1.);
    if (std::abs(pole)<1.e-12) return 0.;
    return sum/pole;
  }

  case ffq2_shape::unknown:
  default:
    msg_Error()<<"Error in "<<METHOD<<": unknown q^2 shape ("
               <<int(m_shape)<<") for form factor '"<<m_name
               <<"'. Falling back to a constant F(0)="<<m_F0<<".\n";
    break;
  }
  return m_F0;
}

Complex FF_P_X::operator()(const ATOOLS::Vec4D_Vector& moms) {
  // q = p_parent - p_daughter. NOTE the sign difference to
  // FF_0_PP/FF_0_PPP, which build Q from the SUM of the two final-state
  // momenta because there the current is a 0->hadrons production
  // current, not a P->X transition.
  Vec4D q  = moms[m_pi[0]]-moms[m_pi[1]];
  double q2 = q.Abs2();
  double mX = sqrt(Max(0.,moms[m_pi[1]].Abs2()));
  return Complex(Value(q2,mX),0.);
}

DECLARE_FF_GETTER(FF_P_X,"FF_P_X")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,FF_P_X>::
operator()(const METOOLS::FF_Parameters &params) const
{
  if (params.m_pi.size()!=2) return NULL;
  return new FF_P_X(params);
}
