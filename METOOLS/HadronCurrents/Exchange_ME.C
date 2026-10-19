#include "METOOLS/HadronCurrents/Exchange_ME.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

namespace {

  bool HasKey(Scoped_Settings s,const std::string & k)
  {
    const std::vector<std::string> keys = s.GetKeys();
    return std::find(keys.begin(),keys.end(),k)!=keys.end();
  }

  // A fermion line has two legs of the same fermion, differing at most by
  // particle/antiparticle, and couplings are quoted for the particle. A
  // FLAVOUR-CHANGING current (D0 -> K-) has no such thing, which is why the
  // defaults below only serve the neutral bosons and everything else has to
  // declare its couplings.
  Flavour LineFlavour(const Flavour_Vector & flavs,
                      const std::vector<int> & idx)
  {
    for (size_t j(0);j<idx.size();++j)
      if (!flavs[idx[j]].IsAnti()) return flavs[idx[j]];
    return flavs[idx[0]].Bar();
  }

  void DefaultVA(const Flavour & bos,const Flavour & f,double & v,double & a)
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
          "the exchange, explicitly. A flavour-changing current has no "
          "model-derivable coupling in any case.");
  }

  double DefaultCoupling(const Flavour & bos)
  {
    const double e2 = 4.*M_PI*MODEL::s_model->ScalarConstant("alpha_QED");
    if (bos.Kfcode()==kf_photon) return e2;
    if (bos.Kfcode()==kf_Z) {
      const double s2 =
        std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
      return e2/(4.*s2*(1.-s2));            // (g_Z/e)^2
    }
    if (bos.Kfcode()==kf_Wplus) return e2/(8.*
      std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW"))));
    THROW(fatal_error,"No default coupling known for exchanged boson "
          +ToString(bos)+". Give Coupling explicitly.");
  }
}

Exchange_ME::Exchange_ME(const Flavour_Vector& flavs,
                         const vector<int>& indices,
                         const string& name,
                         const size_t nin):
  ME_Base(flavs,indices,name), m_nhel(0), m_nin(nin), m_chan(0)
{}

Exchange_ME::~Exchange_ME()
{
  for (size_t k(0);k<m_x.size();++k) {
    if (m_x[k].p_prop) delete m_x[k].p_prop;
    if (m_x[k].p_j1)   delete m_x[k].p_j1;
    if (m_x[k].p_j2)   delete m_x[k].p_j2;
  }
}

Current_Base * Exchange_ME::MakeCurrent(Scoped_Settings s,
                                        const vector<int> & idx)
{
  const string tag = s["Type"].GetScalarWithOtherDefault<string>("Unknown");
  ME_Parameters pars(m_flavs,idx);
  Current_Base * c = Current_Getter_Function::GetObject(tag,pars);
  if (c==NULL) THROW(fatal_error,"Current '"+tag+"' not found.");

  // Everything but Type and Indices is the current's own: form-factor model
  // selectors, fitted couplings, CKM elements. Passed through verbatim, as
  // Current_ME has always done, so a current needs no knowledge of who is
  // calling it.
  GeneralModel model;
  for (const auto & key: s.GetKeys()) {
    if (key=="Type" || key=="Indices") continue;
    if (s[key].GetItemsCount()>1) {
      model[key+string("_abs")]   = s[key][0].GetScalarWithOtherDefault(-1.0);
      model[key+string("_phase")] = s[key][1].GetScalarWithOtherDefault(0.0);
    }
    else model[key] = s[key].GetScalarWithOtherDefault(-1.0);
  }
  // Which legs are incoming is context the current cannot work out for
  // itself -- leg 0 is the decaying particle in a decay and the first beam
  // in a scattering process. Setting it also runs the current's topology
  // check, so this must happen before the current is used.
  c->SetNIn(m_nin);

  c->SetModelParameters(model);
  return c;
}

void Exchange_ME::AddExchange(Scoped_Settings e)
{
  // No Boson means a CONTACT term: D = 1, and the coupling defaults to
  // G_F/sqrt(2). That is exactly what Current_ME always did, so a decay
  // channel carrying no Boson reproduces its old amplitude bit for bit.
  const int kf = e["Boson"].GetScalarWithOtherDefault<int>(0);
  const bool contact(kf==0);
  const Flavour bos(std::abs(kf),kf<0);

  Scoped_Settings j1 = e["J1"], j2 = e["J2"];
  const vector<int> i1 = j1["Indices"].SetDefault(vector<int>()).GetVector<int>();
  const vector<int> i2 = j2["Indices"].SetDefault(vector<int>()).GetVector<int>();
  if (i1.empty() || i2.empty())
    THROW(fatal_error,"Both J1 and J2 need Indices.");
  for (size_t j(0);j<i1.size();++j)
    if (i1[j]<0 || size_t(i1[j])>=m_flavs.size())
      THROW(fatal_error,"J1 index out of range.");
  for (size_t j(0);j<i2.size();++j)
    if (i2[j]<0 || size_t(i2[j])>=m_flavs.size())
      THROW(fatal_error,"J2 index out of range.");

  // Every exchange joins the same two lines, so they must agree on the
  // momentum transfer. Silently using the first one's q for the rest would
  // be a hard bug to see.
  if (m_x.empty()) {
    m_i1 = i1; m_i2 = i2;
    m_qsign.clear();
    for (size_t j(0);j<m_i1.size();++j)
      m_qsign.push_back(IsIncomingLeg(m_i1[j]) ? -1. : +1.);
  }
  else if (i1!=m_i1 || i2!=m_i2)
    THROW(fatal_error,"All Exchange entries must use the same J1/J2 "
          "indices: they share one momentum transfer.");

  // Couplings: defaulted from the model where it knows them, which is the
  // neutral bosons only, and overridden by anything the card states.
  double v1(0.),a1(0.),v2(0.),a2(0.),g(1.);
  const bool wantdefaults = !contact &&
    !(HasKey(j1,"v") && HasKey(j1,"a") && HasKey(j2,"v") && HasKey(j2,"a"));
  if (wantdefaults) {
    DefaultVA(bos,LineFlavour(m_flavs,i1),v1,a1);
    DefaultVA(bos,LineFlavour(m_flavs,i2),v2,a2);
  }
  if (HasKey(j1,"v")) v1 = j1["v"].GetScalarWithOtherDefault<double>(v1);
  if (HasKey(j1,"a")) a1 = j1["a"].GetScalarWithOtherDefault<double>(a1);
  if (HasKey(j2,"v")) v2 = j2["v"].GetScalarWithOtherDefault<double>(v2);
  if (HasKey(j2,"a")) a2 = j2["a"].GetScalarWithOtherDefault<double>(a2);

  if (contact)                  g = Tools::GF/sqrt(2.);
  else if (!HasKey(e,"Coupling")) g = DefaultCoupling(bos);
  if (HasKey(e,"Coupling")) g = e["Coupling"].GetScalarWithOtherDefault<double>(g);

  Exchange x;
  x.p_prop = contact ? NULL
           : new Gauge_Propagator(bos.Mass(),bos.Width(),bos);
  x.p_j1   = MakeCurrent(j1,i1);
  x.p_j2   = MakeCurrent(j2,i2);
  // v and a reach the current through MakeCurrent's sweep when the card
  // states them; when they came from the model defaults they have to be
  // pushed in separately.
  if (wantdefaults) {
    GeneralModel m1,m2;
    m1["v"]=v1; m1["a"]=a1; x.p_j1->SetModelParameters(m1);
    m2["v"]=v2; m2["a"]=a2; x.p_j2->SetModelParameters(m2);
  }
  x.m_g  = g;
  x.m_M2 = contact ? 0. : sqr(bos.Mass());

  CheckFlavourFlow(contact,bos,i1,i2);

  m_x.push_back(x);

  msg_Debugging()<<"Exchange_ME: "<<(contact?Flavour(kf_none):bos)
                 <<" coupling "<<g<<"\n";
}

// A current joining two legs is one continuous line, and whether that line
// may change flavour is decided by the boson on it -- which is why this
// check lives here and not in the current, who never sees the boson.
//
//   NEUTRAL boson (photon, Z): the line is flavour DIAGONAL. ubar(e) Gamma
//     u(mu) is a flavour-changing neutral current, which does not exist at
//     tree level in the Standard Model.
//   CHARGED boson (W), or a contact term standing in for one: the line
//     changes flavour by design. tau- -> nu_tau and D0 -> K- are exactly
//     this, so neither may be rejected.
//
// Nothing downstream catches a violation. XYZFunc asks only whether a
// flavour is an antiparticle, to pick u or v; the spinor algebra is well
// defined for any two legs and returns a finite number for a process with
// no diagram. The Indices: are the ONLY statement of which legs form a
// line, so an index slip silently computes something else -- and the legs
// are not necessarily in the order the run card lists them, because Sherpa
// sorts them by mass unless Sort_Flavors: 0 is set.
void Exchange_ME::CheckFlavourFlow(const bool contact,
                                   const ATOOLS::Flavour & bos,
                                   const std::vector<int> & i1,
                                   const std::vector<int> & i2) const
{
  // A contact term is a W in disguise, so it carries no constraint.
  if (contact || bos.Charge()!=0.) return;

  for (int w(0);w<2;++w) {
    const std::vector<int> & idx = w ? i2 : i1;
    if (idx.size()!=2) continue;          // many-body currents: not a line
    const ATOOLS::Flavour & fa = m_flavs[idx[0]], & fb = m_flavs[idx[1]];
    if (fa.Kfcode()==fb.Kfcode()) continue;
    THROW(fatal_error,"J"+ToString(w+1)+" of a "+ToString(bos)+" exchange "
          "joins legs "+ToString(idx[0])+" and "+ToString(idx[1])+", which "
          "are "+ToString(fa)+" and "+ToString(fb)+". "+ToString(bos)+" is "
          "electrically neutral, so the line it sits on cannot change "
          "flavour -- this is a flavour-changing neutral current. Check the "
          "Indices: against the leg order of the process, remembering that "
          "Sherpa sorts the legs by mass unless Sort_Flavors: 0 is set, so "
          "they need not be in the order the run card lists them.");
  }
}

void Exchange_ME::SetModelParameters(Scoped_Settings& s)
{
  Scoped_Settings ex = s["Exchange"];
  const size_t n = ex.GetItemsCount();
  if (n==0) THROW(fatal_error,"A Currents: block needs an Exchange: list.");
  for (size_t k(0);k<n;++k) AddExchange(ex[k]);

  // Spin states per leg, counted as Spin_Structure counts them: a massless
  // vector has 2, everything else 2s+1, so a SCALAR HAS ONE. Enumerating
  // 2^N instead would walk 32 combinations for e+e- -> pi+ pi- pi0 where
  // only 4 exist, and hand each pion a spin index it does not have.
  m_nspin.clear();
  m_nhel = 1;
  for (size_t i(0);i<m_flavs.size();++i) {
    const size_t ns = (m_flavs[i].IsVector() && m_flavs[i].IsMassive()==0)
                      ? 2 : size_t(m_flavs[i].IntSpin()+1);
    m_nspin.push_back(ns);
    m_nhel *= ns;
  }

  // Which invariant the exchange carries, by WHICH TWO LEGS J1 PAIRS:
  //   {0,1} or {2,3} -> s      {0,2} or {1,3} -> t      {0,3} or {1,2} -> u
  // Only a 2 -> 2 scattering process has Mandelstam channels at all. A
  // four-leg DECAY (tau- -> nu_tau e- nubar_e) also has m_flavs.size()==4
  // and would otherwise be classified as though it were one.
  m_chan = 0;
  if (m_nin==2 && m_flavs.size()==4 && m_i1.size()==2) {
    const int a = std::min(m_i1[0],m_i1[1]), b = std::max(m_i1[0],m_i1[1]);
    if      ((a==0&&b==1) || (a==2&&b==3)) m_chan = 1;
    else if ((a==0&&b==2) || (a==1&&b==3)) m_chan = 2;
    else if ((a==0&&b==3) || (a==1&&b==2)) m_chan = 3;
    else THROW(fatal_error,"J1 does not pair two legs of one line.");
  }
}

// q^mu = sum over J1's legs with the sign that makes each one outgoing. Two
// incoming legs give q^2 = s; one in and one out gives t or u. For a decay
// leg 0 is the decaying particle, which carries the -1, and q^2 comes out as
// the invariant mass of whatever J2 spans -- which is what it should be.
Vec4D Exchange_ME::Q(const Vec4D_Vector & moms) const
{
  Vec4D q(0.,0.,0.,0.);
  for (size_t j(0);j<m_i1.size();++j) q += m_qsign[j]*moms[m_i1[j]];
  return q;
}

void Exchange_ME::Calculate(const Vec4D_Vector & moms, bool anti)
{
  const Vec4D  q  = Q(moms);
  const double q2 = q.Abs2();

  for (size_t k(0);k<m_x.size();++k) {
    m_x[k].p_j1->Calc(moms,anti);
    m_x[k].p_j2->Calc(moms,anti);
  }

  vector<int> spins(m_flavs.size(),0), s1, s2;
  for (size_t n(0);n<m_nhel;++n) {
    // mixed radix over m_nspin, not binary: identical when every leg is a
    // fermion, which is why the leptonic cases never noticed.
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
      amp += x.m_g*(x.p_prop ? (*x.p_prop)(q2) : Complex(1.,0.))*lor;
    }
    // The AMPLITUDE is kept, phase and all. A decay contracts it with the
    // production density matrix; a cross section squares it. Neither belongs
    // here.
    (*this)[n] = amp;
  }
}
