
// LoopTools is an OPTIONAL Sherpa dependency (SHERPA_ENABLE_LOOPTOOLS, OFF by
// default). This file is compiled either way: without it the getter below
// still matches e+e- -> pi+pi-, but throws an explanatory error instead of
// silently declining and leaving the user with a generic "no virtual ME"
// message that never mentions LoopTools.
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#ifdef USING__LOOPTOOLS

#include "clooptools.h"
#include <complex>

namespace {
  inline double ReC(const ComplexType& z) { return std::real(z); }

  // Shorthands of arXiv:2409.03469 Eqs. (3.14a)-(3.14f).
  inline ComplexType C0e(double s, double me2, double x, double y)
  { return C0(me2, me2, s, x, me2, y); }
  inline ComplexType C0pi(double s, double mp2, double x, double y)
  { return C0(mp2, mp2, s, x, mp2, y); }
  inline ComplexType C0ee(double s, double me2, double x)
  { return C0(me2, me2, s, me2, x, me2); }
  inline ComplexType C0pipi(double s, double mp2, double x)
  { return C0(mp2, mp2, s, mp2, x, mp2); }
  inline ComplexType C0epi(double me2, double mp2, double z, double x)
  { return C0(me2, mp2, z, me2, x, mp2); }
  inline ComplexType D0epi(double s, double me2, double mp2,
                           double z, double x, double y)
  { return D0(me2, me2, mp2, mp2, s, z, x, me2, y, mp2); }
}

namespace {

  void VirtInit() { ltini(); setcmpbits(64); }

  //! arXiv:2409.03469 Eq. (3.16), with kappa from Eq. (3.17).
  double DeltaVISR(double s, double t, double u, double alpha,
                           double me2, double mp2, double lam2)
  {
    const double kappa =
        3.*s/(4.*me2 - s)
      - 4.*me2*(12.*mp2*s - 3.*s*s + 2.*(t-u)*(t-u))
        /((4.*me2 - s)*(4.*mp2*s - s*s + (t-u)*(t-u)));

    const ComplexType r =
        kappa*(B0(s, me2, me2) - B0(me2, lam2, me2))
      + 2.*(2.*me2 - s)*C0ee(s, me2, lam2)
      + 4.*me2*B0i(dbb0, me2, lam2, me2);

    return alpha/(2.*M_PI)*ReC(r);
  }

  //! arXiv:2409.03469 Eq. (3.18).
  double DeltaVFSR(double s, double alpha, double mp2, double lam2)
  {
    const ComplexType r =
        (2.*mp2 - s)*( (4.*mp2 - s)*C0pipi(s, mp2, lam2)
                       - 2.*B0(s, mp2, mp2)
                       + 2.*B0(mp2, lam2, mp2) )
      + 2.*mp2*(4.*mp2 - s)*B0i(dbb0, mp2, lam2, mp2);

    return alpha/M_PI/(4.*mp2 - s)*ReC(r);
  }

  //! arXiv:2409.03469 Eq. (3.19), with kappa_t^0 / kappa_u^0 from Eq. (3.20).
  double DeltaVIFI(double s, double t, double u, double alpha,
                           double me2, double mp2, double lam2)
  {
    const double kt = (2.*mp2 - t + u)*me2 + mp2*mp2 + (mp2 - t)*(mp2 - t) - t*u;
    const double ku = (2.*mp2 - u + t)*me2 + mp2*mp2 + (mp2 - u)*(mp2 - u) - t*u;

    const ComplexType r =
        4.*me2*(t-u)/(4.*me2 - s)*( B0(me2, lam2, me2) - B0(s, lam2, lam2) )
      + (8.*me2*me2 - 8.*me2*s + s*s)/(4.*me2 - s)*(t-u)*C0e(s, me2, lam2, lam2)
      + (2.*mp2 - s)*(t-u)*C0pi(s, mp2, lam2, lam2)
      - 2.*(me2 - t)*(me2 + mp2 - t)*C0epi(me2, mp2, t, lam2)
      + 2.*(me2 - u)*(me2 + mp2 - u)*C0epi(me2, mp2, u, lam2)
      + kt*(me2 + mp2 - t)*D0epi(s, me2, mp2, t, lam2, lam2)
      - ku*(me2 + mp2 - u)*D0epi(s, me2, mp2, u, lam2, lam2);

    return alpha/(2.*M_PI)*4.*s/(4.*mp2*s - s*s + (t-u)*(t-u))*ReC(r);
  }

  /*!
    Total delta_V for the YFS mode in use. The three subsets are separately
    gauge invariant and their IR divergences cancel within each subset
    (arXiv:2409.03469 Sect. 3.2), so taking only the ones YFS is resumming is
    consistent rather than a truncation:

      isr     ISR only
      fsr     FSR only
      isrfsr  ISR + FSR + IFI

    IFI belongs with isrfsr alone: it is the interference between radiation off
    the initial and final legs, so it has no meaning when only one side
    radiates.
  */
  double DeltaV(int yfsmode, double s, double t, double u, double alpha,
                double me2, double mp2, double lam2)
  {
    clearcache();
    Setlambda(lam2);

    double d(0.);
    const bool isr(yfsmode == 1 || yfsmode == 2);   // yfsmode::isr, ::isrfsr
    const bool fsr(yfsmode == 3 || yfsmode == 2);   // yfsmode::fsr, ::isrfsr
    const bool ifi(yfsmode == 2);                   // yfsmode::isrfsr

    if (isr) d += DeltaVISR(s, t, u, alpha, me2, mp2, lam2);
    if (fsr) d += DeltaVFSR(s, alpha, mp2, lam2);
    if (ifi) d += DeltaVIFI(s, t, u, alpha, me2, mp2, lam2);
    return d;
  }

}  // anonymous namespace -- nothing above this line may see a Sherpa header

// ... and nothing below it may see LoopTools' macros. clooptools.h defines
// Conjugate(x) as std::conj(x), which rewrites the declaration of
// ATOOLS::CMatrix::Conjugate() into a syntax error; Re and Im are just as
// likely to collide. Neither is used above, so drop all three here.
#undef Conjugate
#undef Re
#undef Im

#endif  // USING__LOOPTOOLS

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "YFS/Main/YFS_Base.H"


using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;

namespace EXTRAXS {

  /*!
    One-loop QED correction to e+e- -> pi+ pi- in scalar QED, interfered with
    the Born.

    delta_V is taken from arXiv:2409.03469 Sect. 3.2 (see the header comment),
    split into the ISR, FSR and IFI subsets so that only the ones the current
    YFS MODE resums are summed.

    Vacuum polarisation is NOT included, matching the paper. An OpenLoops
    comparison therefore needs photon_selfenergy: 0.
  */
  class PionPionVirtual: public PHASIC::Virtual_ME2_Base {
    double ME2, MP2;
    double m_s, m_t, m_u;
    double m_alpha, m_photonmass;
    int    m_yfsmode, m_check_poles;

  public:
    PionPionVirtual(const Process_Info& pi, const Flavour_Vector& flavs);
    ~PionPionVirtual() {}

    void Calc(const ATOOLS::Vec4D_Vector& momenta);
  };

}// namespace EXTRAXS

using namespace EXTRAXS;

PionPionVirtual::PionPionVirtual(const Process_Info& pi,
                                 const Flavour_Vector& flavs)
    : Virtual_ME2_Base(pi, flavs)
{
#ifdef USING__LOOPTOOLS
  VirtInit();
#endif

  ME2 = sqr(Flavour(kf_e).Mass());
  MP2 = sqr(Flavour(kf_pi).Mass());

  m_mode = 0;          // YFS reads ME_Finite() and multiplies by the Born
  m_IRscale = 1.;
  m_UVscale = 1.;

  Scoped_Settings s{Settings::GetMainSettings()["YFS"]};
  m_photonmass = s["PHOTON_MASS"].Get<double>();
  if (s["Dim_Reg"].Get<int>()) m_photonmass = 0.;
  m_check_poles = s["CHECK_POLES"].Get<int>();
  // Which gauge-invariant subsets to include follows YFS: MODE, so the virtual
  // covers exactly what is being resummed. Kept as an int so DeltaV and the
  // amplitude block above stay free of Sherpa headers; the values are
  // YFS::yfsmode::code (off=0, isr=1, isrfsr=2, fsr=3).
  m_yfsmode = (int)s["MODE"].Get<YFS::yfsmode::code>();
  // Validated rather than defaulted: an unset or unrecognised mode would
  // otherwise select no subset at all and return a silently wrong delta_V
  // instead of failing. (m_yfsmode was briefly left unassigned here, which did
  // exactly that.) yfsmode::off means YFS is not running, so reaching this
  // provider at all is a configuration error.
  if (m_yfsmode != YFS::yfsmode::isr && m_yfsmode != YFS::yfsmode::isrfsr &&
      m_yfsmode != YFS::yfsmode::fsr) {
    THROW(fatal_error, "The internal e+e- -> pi+pi- virtual needs YFS: MODE to "
                       "be one of ISR, FSR or ISRFSR; got '" +
                       ToString(s["MODE"].Get<YFS::yfsmode::code>()) + "'.");
  }
}

void PionPionVirtual::Calc(const Vec4D_Vector& momenta)
{
  m_s = (momenta[0] + momenta[1]).Abs2();
  m_t = (momenta[0] - momenta[2]).Abs2();
  m_u = (momenta[0] - momenta[3]).Abs2();
  m_alpha = MODEL::s_model->ScalarConstant("alpha_QED");

  // YFS::Virtual::Calc_V forms V = (alpha/2pi) * ME_Finite() * Born for
  // m_mode == 0, so a correction delta_V given relative to the Born has to be
  // handed over as delta_V * 2pi/alpha.
  const double tofinite = 2.*M_PI/m_alpha;
  const double lam2 = m_photonmass > 0. ? sqr(m_photonmass) : 0.;

#ifdef USING__LOOPTOOLS
  m_res.Finite() = tofinite * DeltaV(m_yfsmode, m_s, m_t, m_u, m_alpha,
                                     ME2, MP2, lam2);
#else
  m_res.Finite() = 0.;   // unreachable: the getter throws before constructing
#endif

  // The soft divergence here is regulated by the photon mass, not by eps, so
  // there is no pole to report in normal running -- hence IR() = 0 unless the
  // pole check asks for one.
  //
  // delta_V = A*ln(lam^2) + finite, exactly linear in ln(lam^2) up to O(lam^2),
  // and YFS's dim-reg subtraction maps ln(lam^2) -> -1/eps (YFS_Form_Factor.C,
  // `DivArrD massph(0,-1,0,0,0,0)`), so the 1/eps coefficient is -A.
  //
  // A is taken by finite difference rather than from a second, hand-written
  // eikonal formula: the dependence is exactly linear, so two points two
  // e-foldings apart give A to full precision, and -- unlike a transcribed
  // formula -- it cannot drift out of step with the expressions it describes.
  // One extra evaluation, only when CHECK_POLES is on.
  m_res.IR()  = 0.;
  m_res.IR2() = 0.;   // both e and pi are massive: soft pole only, no collinear
#ifdef USING__LOOPTOOLS
  if (m_check_poles == 1 && lam2 > 0.) {
    const double lam2b = lam2*std::exp(-2.);   // ln(lam2) - ln(lam2b) = 2
    const double db = DeltaV(m_yfsmode, m_s, m_t, m_u, m_alpha, ME2, MP2, lam2b);
    const double A  = 0.5*(m_res.Finite()/tofinite - db);
    m_res.IR() = -tofinite * A;
  }
#endif
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::PionPionVirtual, "PionPionVirtual")
Virtual_ME2_Base*
ATOOLS::Getter<PHASIC::Virtual_ME2_Base, PHASIC::Process_Info,
               EXTRAXS::PionPionVirtual>::operator()(const Process_Info& pi) const
{
  if (pi.m_loopgenerator.find("Internal") != 0) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size() != 4) return NULL;
  if (fl[0] == Flavour(kf_e) && fl[1] == fl[0].Bar() &&
      (fl[2].Kfcode() == kf_pi_plus || fl[2].Kfcode() == -kf_pi_plus) &&
      fl[3] == fl[2].Bar()) {
#ifdef USING__LOOPTOOLS
    return new PionPionVirtual(pi, fl);
#else
    THROW(fatal_error,
          "The internal one-loop virtual for e+e- -> pi+pi- needs LoopTools, "
          "and this Sherpa was built without it. Reconfigure with "
          "-DSHERPA_ENABLE_LOOPTOOLS=ON -DLOOPTOOLS_DIRS=<prefix> and rebuild, "
          "or set Loop_Generator to an external provider such as OpenLoops.");
#endif
  }
  return NULL;
}
