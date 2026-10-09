#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Settings.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"

#include <algorithm>

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

//
// A matrix element built as a sum over exchanged bosons, each one joining
// two currents:
//
//    M = sum_V  g_V  D_V(q^2)  [ J1_V . J2_V - (J1_V.q)(J2_V.q)/M_V^2 ]
//
// declared in the process block of the run card:
//
//   PROCESSES:
//   - 11 -11 -> 13 -13:
//       Order: {QCD: 0, EW: 2}
//       Currents:
//         Exchange:
//         - Boson: 22
//           J1: {Type: VA_F_F, Indices: [1, 0]}
//           J2: {Type: VA_F_F, Indices: [2, 3]}
//         - Boson: 23
//           J1: {Type: VA_F_F, Indices: [1, 0]}
//           J2: {Type: VA_F_F, Indices: [2, 3]}
//
// Nothing here knows which particles it is looking at. q is DERIVED from
// the indices of J1 -- the sum of its legs' momenta with the sign that
// makes each outgoing -- so e+ e- -> mu+ mu- (both legs incoming, q^2 = s)
// and e- mu- -> e- mu- (one in, one out, q^2 = t) differ by their index
// assignment and by nothing else. There is no s-channel/t-channel switch
// anywhere in this file.
//
// Couplings default from the model for the photon and the Z, keyed on the
// boson and on the fermion line, and any of v, a or the exchange's overall
// Coupling may be given explicitly to override. Typing them is exact
// control; leaving them out cannot silently disagree with the model.
//

namespace EXTRAXS {

  class XS_Currents : public ME2_Base {
  private:

    // One exchanged boson. It owns its OWN pair of currents, because the
    // couplings belong to the (boson, current) pair and not to the current:
    // the photon wants (v,a) = (Q,0) on the same electron line where the Z
    // wants (v_e, a_e). One current object cannot be both.
    struct Exchange {
      METOOLS::Propagator_Base * p_prop;
      METOOLS::Current_Base    * p_j1;
      METOOLS::Current_Base    * p_j2;
      double                     m_g;    // overall coupling, both vertices
      double                     m_M2;   // 0 => skip the q^mu q^nu term
      Exchange() : p_prop(NULL), p_j1(NULL), p_j2(NULL), m_g(0.), m_M2(0.) {}
    };

    std::vector<Exchange> m_x;

    // Flavours have to outlive the currents: Current_Base keeps a REFERENCE,
    // not a copy. Same trap as in ee_3pi.
    Flavour_Vector   m_ffflavs;
    std::vector<int> m_i1, m_i2;      // index lists of J1 and J2
    std::vector<double> m_qsign;      // -1 for an incoming leg, +1 outgoing
    double           m_alpha;
    double           m_colfac;
    size_t           m_nhel;
    std::vector<size_t> m_nspin;  // spin states per leg, 1 for a scalar
    int              m_chan;   // 1 = s, 2 = t, 3 = u exchange

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

    METOOLS::Current_Base * MakeCurrent(const std::string & tag,
                                        const std::vector<int> & idx,
                                        const double & v,const double & a);
    Flavour LineFlavour(const std::vector<int> & idx) const;
    void    DefaultVA(const Flavour & bos,const Flavour & f,
                      double & v,double & a) const;
    double  DefaultCoupling(const Flavour & bos) const;
    void    AddExchange(const Flavour & bos,
                        const std::string & t1,const std::vector<int> & i1,
                        const std::string & t2,const std::vector<int> & i2,
                        const double & v1,const double & a1,
                        const double & v2,const double & a2,
                        const double & g);
    void    BuildFromSettings(ATOOLS::Scoped_Settings s);
    Vec4D   Q(const Vec4D_Vector & mom) const;

  public:

    XS_Currents(const External_ME_Args& args);
    ~XS_Currents();

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

static bool HasKey(Scoped_Settings s,const std::string & k)
{
  const std::vector<std::string> keys = s.GetKeys();
  return std::find(keys.begin(),keys.end(),k)!=keys.end();
}

METOOLS::Current_Base *
XS_Currents::MakeCurrent(const std::string & tag,const std::vector<int> & idx,
                         const double & v,const double & a)
{
  METOOLS::ME_Parameters pars(m_ffflavs,idx);
  METOOLS::Current_Base * c =
    METOOLS::Current_Getter_Function::GetObject(tag,pars);
  if (c==NULL) THROW(fatal_error,"Current '"+tag+"' not found.");
  METOOLS::GeneralModel model;
  model["v"] = v;
  model["a"] = a;
  c->SetModelParameters(model);
  return c;
}

// A fermion line has two legs of the same fermion, differing at most by
// particle/antiparticle. Couplings are quoted for the particle, so take
// whichever leg is not the antiparticle; if both are particles (a t-channel
// line) either will do.
Flavour XS_Currents::LineFlavour(const std::vector<int> & idx) const
{
  for (size_t j(0);j<idx.size();++j)
    if (!m_flavs[idx[j]].IsAnti()) return m_flavs[idx[j]];
  return m_flavs[idx[0]].Bar();
}

void XS_Currents::DefaultVA(const Flavour & bos,const Flavour & f,
                            double & v,double & a) const
{
  if (bos.Kfcode()==kf_photon) { v = f.Charge(); a = 0.; return; }
  if (bos.Kfcode()==kf_Z) {
    const double s2 =
      std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
    a = f.IsoWeak();
    v = a - 2.*f.Charge()*s2;
    return;
  }
  THROW(fatal_error,"No default couplings known for exchanged boson "
        +ToString(bos)+". Give v and a on both currents, and Coupling on "
        "the exchange, explicitly.");
}

double XS_Currents::DefaultCoupling(const Flavour & bos) const
{
  const double e2 = 4.*M_PI*m_alpha;
  if (bos.Kfcode()==kf_photon) return e2;
  if (bos.Kfcode()==kf_Z) {
    const double s2 =
      std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
    return e2/(4.*s2*(1.-s2));   // (g_Z/e)^2
  }
  THROW(fatal_error,"No default coupling known for exchanged boson "
        +ToString(bos)+". Give Coupling explicitly.");
}

void XS_Currents::AddExchange(const Flavour & bos,
                              const std::string & t1,const std::vector<int> & i1,
                              const std::string & t2,const std::vector<int> & i2,
                              const double & v1,const double & a1,
                              const double & v2,const double & a2,
                              const double & g)
{
  // Every exchange joins the same two lines, so they must all agree on the
  // momentum transfer. Silently using the first one's q for the rest would
  // be a hard bug to see in a cross section.
  if (m_x.empty()) {
    m_i1 = i1; m_i2 = i2;
    m_qsign.clear();
    for (size_t j(0);j<m_i1.size();++j)
      m_qsign.push_back(m_i1[j]<2 ? -1. : +1.);
  }
  else if (i1!=m_i1 || i2!=m_i2)
    THROW(fatal_error,"All Exchange entries must use the same J1/J2 indices: "
          "they share one momentum transfer.");

  Exchange x;
  x.p_prop = new METOOLS::Gauge_Propagator(bos.Mass(),bos.Width(),bos);
  x.p_j1   = MakeCurrent(t1,i1,v1,a1);
  x.p_j2   = MakeCurrent(t2,i2,v2,a2);
  x.m_g    = g;
  x.m_M2   = sqr(bos.Mass());
  m_x.push_back(x);
}

void XS_Currents::BuildFromSettings(Scoped_Settings s)
{
  Scoped_Settings ex = s["Exchange"];
  const size_t n = ex.GetItemsCount();
  if (n==0)
    THROW(fatal_error,"A Currents: block must contain an Exchange: list.");

  for (size_t k(0);k<n;++k) {
    Scoped_Settings e = ex[k];
    const int kf = e["Boson"].GetScalarWithOtherDefault<int>(0);
    if (kf==0)
      THROW(fatal_error,"Exchange entry "+ToString(k)+" has no Boson.");
    const Flavour bos(std::abs(kf),kf<0);

    Scoped_Settings j1 = e["J1"], j2 = e["J2"];
    const string t1 = j1["Type"].GetScalarWithOtherDefault<string>("VA_F_F");
    const string t2 = j2["Type"].GetScalarWithOtherDefault<string>("VA_F_F");
    const vector<int> i1 = j1["Indices"].SetDefault(vector<int>()).GetVector<int>();
    const vector<int> i2 = j2["Indices"].SetDefault(vector<int>()).GetVector<int>();
    if (i1.empty() || i2.empty())
      THROW(fatal_error,"Exchange entry "+ToString(k)+
            ": both J1 and J2 need Indices.");
    for (size_t j(0);j<i1.size();++j)
      if (i1[j]<0 || size_t(i1[j])>=m_flavs.size())
        THROW(fatal_error,"J1 index out of range in Exchange entry "+ToString(k));
    for (size_t j(0);j<i2.size();++j)
      if (i2[j]<0 || size_t(i2[j])>=m_flavs.size())
        THROW(fatal_error,"J2 index out of range in Exchange entry "+ToString(k));

    double v1,a1,v2,a2;
    DefaultVA(bos,LineFlavour(i1),v1,a1);
    DefaultVA(bos,LineFlavour(i2),v2,a2);
    // GetScalarWithOtherDefault, not SetDefault().Get(): the settings key
    // path does NOT include the list index, so every entry of Exchange:
    // shares the path ...Currents:Exchange:J1:v. Registering a per-entry
    // default there is a conflict as soon as two exchanges want different
    // couplings on the same line -- which is the normal case, since that is
    // the whole reason each boson carries its own currents.
    if (HasKey(j1,"v")) v1 = j1["v"].GetScalarWithOtherDefault<double>(v1);
    if (HasKey(j1,"a")) a1 = j1["a"].GetScalarWithOtherDefault<double>(a1);
    if (HasKey(j2,"v")) v2 = j2["v"].GetScalarWithOtherDefault<double>(v2);
    if (HasKey(j2,"a")) a2 = j2["a"].GetScalarWithOtherDefault<double>(a2);

    double g = DefaultCoupling(bos);
    if (HasKey(e,"Coupling"))
      g = e["Coupling"].GetScalarWithOtherDefault<double>(g);

    if (m_analytic && bos.Kfcode()!=kf_photon) continue;
    AddExchange(bos,t1,i1,t2,i2,v1,a1,v2,a2,g);
    msg_Info()<<"XS_Currents: exchange "<<bos<<", J1 = "<<t1
              <<" (v="<<v1<<", a="<<a1<<"), J2 = "<<t2
              <<" (v="<<v2<<", a="<<a2<<"), coupling "<<g<<"\n";
  }
}

// q^mu = sum over J1's legs, with the sign that makes every momentum
// outgoing. For two incoming legs this is -(p0+p1) and q^2 = s; for one in
// and one out it is p_out - p_in and q^2 = t.
Vec4D XS_Currents::Q(const Vec4D_Vector & mom) const
{
  Vec4D q(0.,0.,0.,0.);
  for (size_t j(0);j<m_i1.size();++j) q += m_qsign[j]*mom[m_i1[j]];
  return q;
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

  Scoped_Settings cur(args.m_currents);
  BuildFromSettings(cur);

  // Spin states per leg, counted the way Spin_Structure counts them: a
  // massless vector has 2, everything else 2s+1. A SCALAR THEREFORE HAS
  // ONE, not two. Enumerating 2^N instead would, for e+e- -> pi+ pi- pi0,
  // walk 32 combinations where only 4 exist and hand each pion a spin
  // index of 1 that its one-entry structure does not have. Every leg of
  // the two leptonic test processes is a fermion, where the two counts
  // coincide, so nothing here notices until a meson appears.
  m_nspin.clear();
  m_nhel = 1;
  for (size_t i(0);i<m_flavs.size();++i) {
    const size_t n = (m_flavs[i].IsVector() && m_flavs[i].IsMassive()==0)
                     ? 2 : size_t(m_flavs[i].IntSpin()+1);
    m_nspin.push_back(n);
    m_nhel *= n;
  }

  // Tell the phase-space integrator which mapping matches this topology.
  // Single_Process::FillIntegrator reads SIntType() as a bitmask: 1 picks
  // S1Channel, 2 T1Channel, 4 U1Channel, and 0 means "none of those", which
  // for more than two outgoing legs gets flat Rambo.
  //
  // This is derived from J1's indices, the same thing q is derived from, so
  // it cannot drift out of step with the topology. Hardcoding it is not a
  // correctness bug -- the cross section stays unbiased either way -- but it
  // samples hardest where the matrix element is smallest: an s-channel
  // mapping on e- mu- -> e- mu- integrated at 0.14% unweighting efficiency
  // against Comix's 2.9% on the same process.
  // Which Mandelstam the exchange carries follows from WHICH TWO legs J1
  // pairs, not from whether it happens to contain leg 2 or leg 3. Testing
  // membership is only right while J1 also holds leg 0, and that stops
  // being true as soon as Sherpa reorders the initial state -- which it
  // does the moment the two incoming legs differ in mass, since its
  // canonical ordering (Order_Flavour) sorts heavier first.
  //
  //   {0,1} or {2,3}  ->  s        q = -(p0+p1)
  //   {0,2} or {1,3}  ->  t        p3-p1 = p0-p2
  //   {0,3} or {1,2}  ->  u        p2-p1 = p0-p3
  //
  // Anything else -- three-leg currents, 2->3 -- leaves m_sintt at 0, which
  // FillIntegrator reads as "use Rambo". That is required, not merely
  // tolerated: S1Channel throws for more than two outgoing legs.
  m_chan  = 0;
  m_sintt = 0;
  if (m_flavs.size()==4 && m_i1.size()==2) {
    const int a = std::min(m_i1[0],m_i1[1]), b = std::max(m_i1[0],m_i1[1]);
    if      ((a==0&&b==1) || (a==2&&b==3)) { m_chan = 1; m_sintt = 1; }
    else if ((a==0&&b==2) || (a==1&&b==3)) { m_chan = 2; m_sintt = 2; }
    else if ((a==0&&b==3) || (a==1&&b==2)) { m_chan = 3; m_sintt = 4; }
    else THROW(fatal_error,"J1 does not pair two legs of one line.");
  }
  msg_Debugging()<<"XS_Currents: SIntType = "<<m_sintt<<"\n";

  // Colour. A boson that is itself a colour singlet leaves each fermion
  // current diagonal in colour, so summing that index gives N_c for every
  // strongly charged LINE -- once, whether the line is incoming or
  // outgoing. (For an incoming one the framework then averages it back out
  // through m_symfac, which already folds in StrongCharge.) Without this
  // e+ e- -> q qbar comes out at exactly 1/3 of the right answer.
  //
  // This is only valid while every exchanged boson is a colour singlet,
  // which is checked below. A gluon, or any octet exchange, needs real
  // colour algebra that the Currents: schema cannot express at all.
  m_colfac = 1.;
  for (size_t k(0);k<m_x.size();++k)
    if (m_x[k].p_prop->Flav().Strong())
      THROW(fatal_error,"XS_Currents: exchanged boson "
            +ToString(m_x[k].p_prop->Flav())+" is coloured. Only colour "
            "singlet exchange is supported; an octet needs colour algebra "
            "the Currents: schema does not carry.");
  for (int line(0);line<2;++line) {
    const std::vector<int> & idx = line ? m_i2 : m_i1;
    if (!LineFlavour(idx).Strong()) continue;
    m_colfac *= 3.;
    // One colour line threaded through the two legs: the quark leg carries
    // the colour index, the antiquark leg the anticolour one.
    if (idx.size()==2) {
      const int a(idx[0]), b(idx[1]);
      const size_t ca(m_flavs[a].IsAnti()?1:0);
      m_colours[a][ca] = m_colours[b][1-ca] = 500;
    }
  }

  p_ref = NULL; m_maxdev = 0.; m_ncheck = 0;
  if (set["XS_CURRENTS_CHECK"].SetDefault(0).Get<int>()) {
    // XS_ee_ffbar is a MASSLESS formula: it reads m_flavs[2].Mass() in its
    // constructor and never uses it again. Comparing a massive final state
    // against it reports a large deviation that says nothing about this
    // code -- it cost an hour the first time. Refuse the comparison rather
    // than hand back a number that looks like a failure.
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
  for (size_t k(0);k<m_x.size();++k) {
    if (m_x[k].p_prop) delete m_x[k].p_prop;
    if (m_x[k].p_j1)   delete m_x[k].p_j1;
    if (m_x[k].p_j2)   delete m_x[k].p_j2;
  }
}

double XS_Currents::operator()(const ATOOLS::Vec4D_Vector& mom)
{
  const Vec4D  q  = Q(mom);
  const double q2 = q.Abs2();
  if (!(std::abs(q2)>0.)) return 0.;

  for (size_t k(0);k<m_x.size();++k) {
    m_x[k].p_j1->Calc(mom,false);
    m_x[k].p_j2->Calc(mom,false);
  }

  // helicity sum. Each current spans its own legs, so the global spin
  // vector is sliced per current exactly as Current_ME does it.
  double sum(0.);
  std::vector<int> spins(m_flavs.size(),0), s1, s2;
  for (size_t n(0);n<m_nhel;++n) {
    // mixed radix over m_nspin, not binary; identical when every leg is a
    // fermion, which is why the leptonic cards are unaffected.
    size_t m(n);
    for (size_t i(0);i<spins.size();++i) {
      spins[i] = m % m_nspin[i];
      m /= m_nspin[i];
    }
    Complex amp(0.,0.);
    for (size_t k(0);k<m_x.size();++k) {
      const Exchange & x = m_x[k];
      s1.clear(); s2.clear();
      for (size_t j(0);j<x.p_j1->DecayIndices().size();++j)
        s1.push_back(spins[x.p_j1->DecayIndices()[j]]);
      for (size_t j(0);j<x.p_j2->DecayIndices().size();++j)
        s2.push_back(spins[x.p_j2->DecayIndices()[j]]);
      const Vec4C J1 = x.p_j1->Get(s1), J2 = x.p_j2->Get(s2);
      Complex lor = J1*J2;
      if (x.m_M2>0.) lor -= (J1*q)*(J2*q)/x.m_M2;
      amp += x.m_g*(*x.p_prop)(q2)*lor;
    }
    sum += std::norm(amp);
  }

  // m_symfac is ISSymmetryFactor(in) * FSSymmetryFactor(out). Despite the
  // name the first is not a symmetry factor at all: it is the product of
  // (2s+1) over the incoming legs, times their colour charge, i.e. exactly
  // the spin/colour averaging denominator (4 for e+ e-, 4*9 for a q qbar
  // initial state). ME2_Base::Calc multiplies it straight back out, so an
  // ME is expected to hand back the AVERAGED value and let the framework
  // undo it. Using m_symfac rather than a hardcoded 1/4 is what makes this
  // correct for coloured currents and for identical final-state particles,
  // neither of which these processes exercise.
  const double res = m_colfac*sum/m_symfac;

  if (m_analytic && m_flavs.size()==4) {
    // Massless QED, spin-averaged. The two topologies are related by
    // crossing s <-> t, which is exactly the statement that the only
    // difference upstream was which legs went into which current.
    const double s = (mom[0]+mom[1]).Abs2();
    const double u = (mom[0]-mom[3]).Abs2();
    const double tt= (mom[0]-mom[2]).Abs2();
    const double e2 = 4.*M_PI*m_alpha, e4 = e2*e2;
    // Annihilation with BOTH pairs massive:
    //   <|M|^2> = 2e^4/s^2 [ (t-M)^2 + (u-M)^2 + 2 s (m1^2 + m3^2) ],
    // with M = m1^2 + m3^2. Reduces to 2e^4(t^2+u^2)/s^2 when both are
    // massless, and to the massless-electron form when only m1 = 0.
    //
    // XS_ee_ffbar cannot be the reference for any of this: it reads
    // m_flavs[2].Mass() in its constructor and then never uses it, so it
    // is a massless matrix element whatever the particle data says.
    // One formula for all three channels: with X the invariant the
    // exchange carries and Y, Z the other two,
    //
    //   <|M|^2> = 2 e^4 / X^2 [ (Y-M)^2 + (Z-M)^2 + 2 X M ],
    //
    // M being the summed squared masses of the two fermion LINES -- taken
    // from LineFlavour, so it does not care how the legs are ordered.
    // Reduces to 2e^4(Y^2+Z^2)/X^2 at zero mass.
    const double M = sqr(LineFlavour(m_i1).Mass())
                   + sqr(LineFlavour(m_i2).Mass());
    double X(s), Y(tt), Z(u);
    if      (m_chan==2) { X = tt; Y = s; Z = u; }
    else if (m_chan==3) { X = u;  Y = s; Z = tt; }
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
