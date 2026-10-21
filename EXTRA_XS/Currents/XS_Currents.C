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
    // XS_CURRENTS_ANALYTIC for e+e- -> f fbar gamma: the legs, found once.
    // -1 when the process is not of that shape, and the check is skipped.
    int                     m_ia, m_ib, m_ic, m_id, m_ig;
    // XS_CURRENTS_GAUGE: eps* -> k at every point. The amplitude must
    // vanish; the worst |M(k)| / (k^0 |M|) seen is reported.
    bool                    m_gauge;
    double                  m_gdev;
    long                    m_gcheck;
    // A maximum alone cannot tell rounding in a soft or collinear corner
    // from a small systematic error, so both radiative checks also count
    // deviations by decade (bin d holds [1e-(d+1), 1e-d); the last bin is
    // everything below 1e-12) and record where the worst point sits: photon
    // energy fraction 2k0/sqrt(s) and the smallest 1-cos(angle) to any
    // charged leg.
    std::vector<long>       m_ahist, m_ghist;
    double                  m_aworst[2], m_gworst[2];
    void Record(const double dev,double & max,std::vector<long> & hist,
                double * worst,const ATOOLS::Vec4D_Vector & mom) const;
    void Report(const std::string & what,const double max,const long n,
                const std::vector<long> & hist,const double * worst) const;
    void   FindRadiativeLegs();
    double AnalyticRadiative(const ATOOLS::Vec4D_Vector& mom) const;

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
    // The LINE, photons removed: a photon carries no colour, and every
    // diagram joins the same two lines, so this is done once.
    const std::vector<int> & idx = line ? p_me->Line2() : p_me->Line1();
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
  // Set from the card once the block has been read, below.
  m_oew   = 2;
  m_oqcd  = 0;
  for (size_t i(0);i<m_flavs.size();++i)
    m_colours[i][0] = m_colours[i][1] = 0;

  m_ffflavs = m_flavs;
  m_alpha   = MODEL::s_model->ScalarConstant("alpha_QED");

  Settings & set = Settings::GetMainSettings();
  m_analytic = set["XS_CURRENTS_ANALYTIC"].SetDefault(0).Get<int>();
  m_adev = 0.; m_ancheck = 0;
  m_gauge    = set["XS_CURRENTS_GAUGE"].SetDefault(0).Get<int>();
  m_gdev = 0.; m_gcheck = 0;
  m_ahist.assign(13,0); m_ghist.assign(13,0);
  m_aworst[0] = m_aworst[1] = m_gworst[0] = m_gworst[1] = 0.;
  m_ia = m_ib = m_ic = m_id = m_ig = -1;

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

  // COUPLING ORDERS come from the card's Order: block, which Sherpa treats
  // as exact; 99 means it was not stated. They are also DERIVABLE from what
  // was declared -- 2 for the exchange, one more per photon -- so the two
  // are compared rather than either being trusted. A card saying EW: 2 for
  // a process with a photon has been edited from its photon-less parent and
  // half-updated; it should not quietly integrate with the wrong coupling
  // power.
  //
  // QCD is 0 by construction: the exchange is a colour singlet (coloured
  // ones are refused in SetColours) and no current here has a strong
  // vertex.
  const int oew = p_me->OrderEW();
  if (args.m_orders.size()>1 && args.m_orders[1]<99. &&
      int(args.m_orders[1]+0.5)!=oew)
    THROW(fatal_error,"XS_Currents: the process' Order: block says EW: "
          +ToString(args.m_orders[1])+", but the Currents: block declares "
          +ToString(oew)+" -- 2 for the exchange plus one per photon ("
          +ToString(p_me->Photons().size())+" here). Make them agree.");
  if (args.m_orders.size()>0 && args.m_orders[0]<99. && args.m_orders[0]>0.)
    THROW(fatal_error,"XS_Currents: the process' Order: block says QCD: "
          +ToString(args.m_orders[0])+", but a matrix element built from "
          "currents has no strong vertex. Use QCD: 0.");
  m_oew  = oew;
  m_oqcd = 0;

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
  if (m_analytic) FindRadiativeLegs();

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
  if (m_analytic && m_ancheck) {
    if (m_ig>=0)
      Report("vs analytic e+e- -> f fbar gamma, relative deviation,",m_adev,
             m_ancheck,m_ahist,m_aworst);
    else
      msg_Out()<<"### XS_Currents vs analytic QED over "<<m_ancheck
               <<" points: max relative deviation "<<m_adev<<"\n";
  }
  if (m_gauge && m_gcheck)
    Report("gauge check (eps* -> k), |M(k)|/(k0 |M|),",m_gdev,m_gcheck,
           m_ghist,m_gworst);
  if (p_ref) {
    msg_Out()<<"### XS_Currents vs XS_ee_ffbar over "<<m_ncheck
             <<" points: max relative deviation "<<m_maxdev<<"\n";
    delete p_ref;
  }
  if (p_me) delete p_me;
}

double XS_Currents::operator()(const ATOOLS::Vec4D_Vector& mom)
{
  // A point on a massless propagator's pole is now refused inside
  // Exchange_ME, per exchange: with a photon in one current, different
  // diagrams have different q and there is no single q^2 to test here.
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
    const double M = sqr(LineFlavour(p_me->Line1()).Mass())
                   + sqr(LineFlavour(p_me->Line2()).Mass());
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

  if (m_analytic && m_ig>=0) {
    const double ana = AnalyticRadiative(mom);
    const double den = std::max(std::abs(ana),std::abs(res));
    if (den>0.) Record(std::abs(ana-res)/den,m_adev,m_ahist,m_aworst,mom);
    ++m_ancheck;
  }

  // Last, because it overwrites the amplitudes held by p_me.
  if (m_gauge && !p_me->Photons().empty() && sum>0.) {
    p_me->SetPhotonGaugeTest(true);
    p_me->Calculate(mom,false);
    p_me->SetPhotonGaugeTest(false);
    // Normalised by each photon's energy in the partonic rest frame, so the
    // figure does not depend on the frame the momenta arrive in.
    double sk(0.), k0(1.);
    const Vec4D P = mom[0]+mom[1];
    for (size_t n(0);n<p_me->size();++n) sk += std::norm((*p_me)[n]);
    for (size_t j(0);j<p_me->Photons().size();++j)
      k0 *= (mom[p_me->Photons()[j]]*P)/P.Abs();
    Record(sqrt(sk/sum)/k0,m_gdev,m_ghist,m_gworst,mom);
    ++m_gcheck;
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

void XS_Currents::Record(const double dev,double & max,
                         std::vector<long> & hist,double * worst,
                         const ATOOLS::Vec4D_Vector & mom) const
{
  const int d = dev>0. ? std::min(12,std::max(0,int(-std::log10(dev))))
                       : 12;
  ++hist[d];
  if (dev<=max) return;
  max = dev;
  const std::vector<int> & g = p_me->Photons();
  if (g.empty()) return;
  // In the partonic rest frame, written invariantly: with a lepton PDF the
  // momenta here are in the lab frame, and k0/sqrt(s_hat) is meaningless.
  const Vec4D & k = mom[g[0]];
  const Vec4D   P = mom[0]+mom[1];
  worst[0] = 2.*(k*P)/P.Abs2();
  worst[1] = 2.;
  for (size_t i(0);i<mom.size();++i) {
    if (m_flavs[i].IntCharge()==0) continue;
    worst[1] = std::min(worst[1],1.-mom[i].CosTheta(k));
  }
}

void XS_Currents::Report(const std::string & what,const double max,
                         const long n,const std::vector<long> & hist,
                         const double * worst) const
{
  msg_Out()<<"### XS_Currents "<<what<<" over "<<n<<" points: max "<<max
           <<"\n###   worst point: 2k0/sqrt(s) = "<<worst[0]
           <<", smallest 1-cos(angle) to a charged leg = "<<worst[1]
           <<"\n###   points per decade:";
  for (int d(0);d<13;++d) {
    if (!hist[d]) continue;
    if (d<12) msg_Out()<<"  [1e-"<<d+1<<",1e-"<<d<<"): "<<hist[d];
    else      msg_Out()<<"  <1e-12: "<<hist[d];
  }
  msg_Out()<<"\n";
}

// The reference for XS_CURRENTS_ANALYTIC with one photon needs an incoming
// antiparticle a and particle b of one massless flavour, an outgoing
// antiparticle c and particle d of another, and the photon g: e+e- ->
// mu+mu- gamma and its relatives (tau pairs, q qbar -> l+l- gamma, e+e- ->
// q qbar gamma). Anything else skips the check, saying so.
void XS_Currents::FindRadiativeLegs()
{
  if (m_flavs.size()!=5 || p_me->Photons().size()!=1) return;
  std::vector<int> in, out;
  for (int i(0);i<5;++i) {
    if (m_flavs[i].Kfcode()==kf_photon) { m_ig = i; continue; }
    (i<2 ? in : out).push_back(i);
  }
  const bool ok = m_ig>=2 && in.size()==2 && out.size()==2
    && m_flavs[in[0]].Kfcode()==m_flavs[in[1]].Kfcode()
    && m_flavs[in[0]].IsAnti()!=m_flavs[in[1]].IsAnti()
    && m_flavs[out[0]].Kfcode()==m_flavs[out[1]].Kfcode()
    && m_flavs[out[0]].IsAnti()!=m_flavs[out[1]].IsAnti()
    && m_flavs[in[0]].Mass()==0.;
  if (!ok) {
    m_ig = -1;
    msg_Info()<<"XS_Currents: no analytic reference for this radiative "
              <<"process (it covers f fbar -> f' fbar' gamma with massless "
              <<"incoming fermions); comparison skipped.\n";
    return;
  }
  m_ia = m_flavs[in[0]].IsAnti()  ? in[0]  : in[1];
  m_ib = m_flavs[in[0]].IsAnti()  ? in[1]  : in[0];
  m_ic = m_flavs[out[0]].IsAnti() ? out[0] : out[1];
  m_id = m_flavs[out[0]].IsAnti() ? out[1] : out[0];
}

// e+(p+) e-(p-) -> fbar(q+) f(q-) gamma(k) through photon exchange only,
// massless incoming fermions, any final-state mass; averaged over incoming
// and summed over outgoing spins.
//
// Berends, Kleiss, Jadach, Was, Acta Phys. Pol. B14 (1983) 413, Eqs.
// (3.2)-(3.7), in the dimensionful form of Silagadze, arXiv:2010.15553,
// Eqs. (48)-(50):
//
//   |M|^2 = e^6 Q^2 Q'^2 [ Q^2 R_ini + Q'^2 R_fin + Q Q' R_int ],
//
// Q and Q' the charges of the incoming and outgoing ANTIparticle. The
// arXiv paper prints "y^2" in R_int where the original has u^2; u^2 is
// right, and with it the two papers agree to machine precision for massive
// final states, and both reduce to the massless closed form of
// arXiv:2010.15553 Eq. (46) as m -> 0 (all checked numerically).
//
// Coupling: e^4 from the photon exchange (as built, so a Coupling: override
// is honoured) and e^2 at the emission vertex (Photon_Alpha). Colour: N_c
// summed over an outgoing coloured line, averaged out (1/N_c) over an
// incoming one.
double XS_Currents::AnalyticRadiative(const ATOOLS::Vec4D_Vector& mom) const
{
  const Vec4D & pp = mom[m_ia], & pm = mom[m_ib];
  const Vec4D & qp = mom[m_ic], & qm = mom[m_id], & k = mom[m_ig];
  const double s  = (pp+pm).Abs2(), sp = (qp+qm).Abs2();
  const double t  = -2.*(pp*qp), u  = -2.*(pp*qm);
  const double tp = -2.*(pm*qm), up = -2.*(pm*qp);
  const double x1 = k*pp, x2 = k*pm, y1 = k*qp, y2 = k*qm;
  const double m2 = sqr(m_flavs[m_ic].Mass());

  const double Ri = (t*t+u*u+tp*tp+up*up
                     +2.*m2*(sqr(t+u)+sqr(tp+up))/sp)/(sp*x1*x2);
  const double Rf = ((t*t+up*up+2.*m2*s)*(1.-2.*m2/s*(1.+y1/y2))
                    +(tp*tp+u*u+2.*m2*s)*(1.-2.*m2/s*(1.+y2/y1))
                    +8.*m2/s*(x1*x1+x2*x2)-8.*m2*(s-sp))/(s*y1*y2);
  const double Rn = ((u*x2*y1+up*x1*y2-t*x2*y2-tp*x1*y1)
                     *(t*t+tp*tp+u*u+up*up+2.*m2*(s+sp))
                    +2.*m2*x1*x2*((s-sp)*(u+up-t-tp)-4.*(x1-x2)*(y1-y2)))
                    /(s*sp*x1*x2*y1*y2);

  const double Q = m_flavs[m_ia].Charge(), Qp = m_flavs[m_ic].Charge();
  double gex(0.);
  const std::vector<METOOLS::Exchange_ME::Exchange> & xs = p_me->Exchanges();
  for (size_t j(0);j<xs.size();++j)
    if (xs[j].p_prop && xs[j].p_prop->Flav().Kfcode()==kf_photon &&
        xs[j].m_g!=0.) { gex = xs[j].m_g; break; }
  const double e2g = 4.*M_PI*p_me->PhotonAlpha();
  double col(1.);
  if (m_flavs[m_ia].Strong()) col /= 3.;
  if (m_flavs[m_ic].Strong()) col *= 3.;
  return col*gex*gex*e2g*Q*Q*Qp*Qp*(Q*Q*Ri+Qp*Qp*Rf+Q*Qp*Rn);
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
