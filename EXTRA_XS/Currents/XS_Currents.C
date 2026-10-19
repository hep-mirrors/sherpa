#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "METOOLS/HadronCurrents/Exchange_ME.H"

#include <algorithm>

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

//
// A cross section from a Currents: block in the process declaration.
//
// The amplitude is NOT built here. It is built by METOOLS::Exchange_ME,
// which a hadron decay channel uses through the same block -- so a current
// written for a decay works in a scattering process unchanged, and there is
// one implementation of
//
//   M = sum_V  g_V  D_V(q^2) [ J1_V . J2_V - (J1_V.q)(J2_V.q)/M_V^2 ]
//
// rather than two that have to be kept in step.
//
// What is left here is what is genuinely EXTRA_XS and has no meaning for a
// decay: squaring the amplitude, the spin and colour averaging, the colour
// flow, the phase-space channel hint, and the getter.
//

namespace EXTRAXS {

  class XS_Currents : public ME2_Base {
  private:

    METOOLS::Exchange_ME * p_me;
    // Exchange_ME keeps a REFERENCE to the flavours, so they have to outlive
    // it. Same trap as in ee_3pi.
    Flavour_Vector m_ffflavs;
    double         m_colfac, m_alpha;

    // Validation scaffolds, both default off. XS_CURRENTS_CHECK builds the
    // hardcoded reference ME alongside and compares at every phase-space
    // point; XS_CURRENTS_ANALYTIC drops the Z and compares against the
    // closed-form massless QED result. An integrated cross section can hide
    // a wrong angular dependence; neither of these can.
    PHASIC::Tree_ME2_Base * p_ref;
    double                  m_maxdev;
    long                    m_ncheck;
    bool                    m_analytic;
    double                  m_adev;
    long                    m_ancheck;

    Flavour LineFlavour(const std::vector<int> & idx) const;
    void    SetColours();

  public:

    XS_Currents(const External_ME_Args& args);
    ~XS_Currents();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

// A line's two legs are the same fermion up to charge conjugation, and the
// particle is the one quoted. A flavour-changing current has no such thing,
// but nothing below needs one: this is only used to ask whether the line is
// coloured, which is false for every hadron.
Flavour XS_Currents::LineFlavour(const std::vector<int> & idx) const
{
  for (size_t j(0);j<idx.size();++j)
    if (!m_flavs[idx[j]].IsAnti()) return m_flavs[idx[j]];
  return m_flavs[idx[0]].Bar();
}

void XS_Currents::SetColours()
{
  // A boson that is itself a colour singlet leaves each fermion current
  // diagonal in colour, so summing that index gives N_c for every strongly
  // charged LINE -- once, whether the line is incoming or outgoing. For an
  // incoming one the framework averages it back out through m_symfac, which
  // already folds in StrongCharge. Without this, e+ e- -> q qbar comes out
  // at exactly 1/3 of the right answer, silently.
  m_colfac = 1.;
  const std::vector<METOOLS::Exchange_ME::Exchange> & x = p_me->Exchanges();
  for (size_t k(0);k<x.size();++k)
    if (x[k].p_prop && x[k].p_prop->Flav().Strong())
      THROW(fatal_error,"XS_Currents: exchanged boson "
            +ToString(x[k].p_prop->Flav())+" is coloured. Only colour "
            "singlet exchange is supported; an octet needs colour algebra "
            "the Currents: schema does not carry.");

  for (int line(0);line<2;++line) {
    const std::vector<int> & idx = line ? p_me->J2Indices() : p_me->J1Indices();
    if (idx.empty() || !LineFlavour(idx).Strong()) continue;
    m_colfac *= 3.;
    // One colour line threaded through the two legs: the quark leg carries
    // the colour index, the antiquark leg the anticolour one.
    if (idx.size()==2) {
      const int a(idx[0]), b(idx[1]);
      const size_t ca(m_flavs[a].IsAnti()?1:0);
      m_colours[a][ca] = m_colours[b][1-ca] = 500;
    }
  }
}

XS_Currents::XS_Currents(const External_ME_Args& args) : ME2_Base(args)
{
  m_oew   = 2;
  m_oqcd  = 0;
  for (size_t i(0);i<m_flavs.size();++i)
    m_colours[i][0] = m_colours[i][1] = 0;

  m_ffflavs = m_flavs;
  m_alpha   = MODEL::s_model->ScalarConstant("alpha_QED");

  Settings & set = Settings::GetMainSettings();
  m_analytic = set["XS_CURRENTS_ANALYTIC"].SetDefault(0).Get<int>();
  m_adev = 0.; m_ancheck = 0;

  std::vector<int> pi(m_flavs.size());
  for (size_t i(0);i<pi.size();++i) pi[i] = i;
  // 2: this is a scattering process, so legs 0 and 1 are the beams. The
  // decay path leaves this at its default of 1, where leg 0 alone is
  // incoming. Everything that depends on which legs are incoming -- the
  // sign of q, the Mandelstam classification, each current's topology
  // check -- follows from it.
  p_me = new METOOLS::Exchange_ME(m_ffflavs,pi,"XS_Currents",2);
  Scoped_Settings cur(args.m_currents);
  p_me->SetModelParameters(cur);

  // The closed form compared against below is pure QED, so under
  // XS_CURRENTS_ANALYTIC every exchange but the photon is switched off by
  // zeroing its coupling. Done from here rather than through a flag inside
  // Exchange_ME, which has no business knowing about this code's test
  // scaffolds.
  if (m_analytic) {
    std::vector<METOOLS::Exchange_ME::Exchange> & xs = p_me->Exchanges();
    for (size_t k(0);k<xs.size();++k)
      if (xs[k].p_prop==NULL || xs[k].p_prop->Flav().Kfcode()!=kf_photon)
        xs[k].m_g = 0.;
  }

  SetColours();

  // Which phase-space mapping matches this topology. Exchange_ME classifies
  // the exchange by which two legs J1 pairs; FillIntegrator reads SIntType
  // as a bitmask, and 0 means "none of S1/T1/U1", which for more than two
  // outgoing legs gets flat Rambo. That is required, not merely tolerated:
  // S1Channel throws for nout > 2.
  switch (p_me->Channel()) {
  case 1:  m_sintt = 1; break;     // q^2 = s
  case 2:  m_sintt = 2; break;     // q^2 = t
  case 3:  m_sintt = 4; break;     // q^2 = u
  default: m_sintt = 0; break;
  }

  p_ref = NULL; m_maxdev = 0.; m_ncheck = 0;
  if (set["XS_CURRENTS_CHECK"].SetDefault(0).Get<int>()) {
    // XS_ee_ffbar is a MASSLESS formula: it reads m_flavs[2].Mass() in its
    // constructor and never uses it again. Comparing a massive final state
    // against it reports a large deviation that says nothing about this
    // code -- it cost an hour the first time. Refuse rather than mislead.
    //
    // It is kept for the massless case because it is the only check on the
    // Z couplings and the gamma/Z interference; XS_CURRENTS_ANALYTIC drops
    // the Z and so cannot see them.
    bool massive(false);
    for (size_t i(0);i<m_flavs.size();++i)
      if (m_flavs[i].Mass()>0.) massive = true;
    if (massive) {
      msg_Info()<<"XS_Currents: XS_ee_ffbar is a massless matrix element, "
                <<"so it is not a valid reference for this final state; "
                <<"comparison skipped. Use XS_CURRENTS_ANALYTIC=1.\n";
    }
    else {
      p_ref = ATOOLS::Getter_Function<PHASIC::Tree_ME2_Base,
                                      PHASIC::External_ME_Args>::
              GetObject("XS_ee_ffbar",args);
      if (p_ref==NULL) msg_Info()<<"XS_Currents: no hardcoded reference for "
                                   "this process; comparison skipped.\n";
    }
  }
}

XS_Currents::~XS_Currents()
{
  if (m_analytic && m_ancheck)
    msg_Out()<<"### XS_Currents vs analytic QED over "<<m_ancheck
             <<" points: max relative deviation "<<m_adev<<"\n";
  if (p_ref) {
    msg_Out()<<"### XS_Currents vs XS_ee_ffbar over "<<m_ncheck
             <<" points: max relative deviation "<<m_maxdev<<"\n";
    delete p_ref;
  }
  if (p_me) delete p_me;
}

double XS_Currents::operator()(const ATOOLS::Vec4D_Vector& mom)
{
  const double q2 = p_me->Q(mom).Abs2();
  if (!(std::abs(q2)>0.)) return 0.;

  p_me->Calculate(mom,false);

  double sum(0.);
  for (size_t n(0);n<p_me->size();++n) sum += std::norm((*p_me)[n]);

  // m_symfac is ISSymmetryFactor(in) * FSSymmetryFactor(out). Despite the
  // name the first is not a symmetry factor at all: it is the product of
  // (2s+1) over the incoming legs, times their colour charge, i.e. exactly
  // the spin/colour averaging denominator. ME2_Base::Calc multiplies it
  // straight back out, so an ME hands back the AVERAGED value and lets the
  // framework undo it.
  const double res = m_colfac*sum/m_symfac;

  if (m_analytic && m_flavs.size()==4) {
    const double s = (mom[0]+mom[1]).Abs2();
    const double u = (mom[0]-mom[3]).Abs2();
    const double tt= (mom[0]-mom[2]).Abs2();
    const double e2 = 4.*M_PI*m_alpha, e4 = e2*e2;
    // One formula for all three channels: with X the invariant the exchange
    // carries and Y, Z the other two,
    //   <|M|^2> = 2 e^4 / X^2 [ (Y-M)^2 + (Z-M)^2 + 2 X M ],
    // M being the summed squared masses of the two fermion LINES, so it
    // does not care how the legs are ordered. Reduces to 2e^4(Y^2+Z^2)/X^2
    // at zero mass.
    const double M = sqr(LineFlavour(p_me->J1Indices()).Mass())
                   + sqr(LineFlavour(p_me->J2Indices()).Mass());
    double X(s), Y(tt), Z(u);
    if      (p_me->Channel()==2) { X = tt; Y = s; Z = u; }
    else if (p_me->Channel()==3) { X = u;  Y = s; Z = tt; }
    const double ana = 2.*e4*(sqr(Y-M)+sqr(Z-M)+2.*X*M)/(X*X);
    const double den = std::max(std::abs(ana),std::abs(res));
    if (den>0.) {
      const double dev = std::abs(ana-res)/den;
      if (dev>m_adev) m_adev = dev;
    }
    ++m_ancheck;
  }

  if (p_ref) {
    ME2_Base * rme = dynamic_cast<ME2_Base*>(p_ref);
    const double ref = rme ? (*rme)(mom) : p_ref->Calc(mom);
    const double den = std::max(std::abs(ref),std::abs(res));
    if (den>0.) {
      const double dev = std::abs(ref-res)/den;
      if (dev>m_maxdev) m_maxdev = dev;
    }
    ++m_ncheck;
  }
  return res;
}

// DECLARE_TREEME2_GETTER stringifies its tag for PrintInfo, so it has to be
// a literal here. This pins that literal to the constant Single_Process
// looks the ME up by: if the two ever drift, the build stops rather than
// the lookup quietly failing and a hardcoded ME answering instead.
namespace {
  constexpr bool SameTag(const char * a,const char * b)
  { return *a==*b && (*a=='\0' || SameTag(a+1,b+1)); }
}
static_assert(SameTag(EXTRAXS::s_currents_tag,"XS_Currents"),
              "getter tag must match EXTRAXS::s_currents_tag");

DECLARE_TREEME2_GETTER(EXTRAXS::XS_Currents,"XS_Currents")
Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,
                              EXTRAXS::XS_Currents>::
operator()(const External_ME_Args &args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  // A Currents: block is the only way to this matrix element. Without one
  // there is nothing to build from, so decline and let the usual getter
  // search run -- it is the hardcoded MEs' job from here.
  if (!args.m_hascurrents) return NULL;
  return new XS_Currents(args);
}
