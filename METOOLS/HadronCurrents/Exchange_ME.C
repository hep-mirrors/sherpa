#include "METOOLS/HadronCurrents/Exchange_ME.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/STL_Tools.H"
#include <algorithm>

using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

namespace {

  // "boson 22 (P)": Sherpa's short name alone is cryptic in a message.
  std::string BosonName(const int kf)
  {
    if (kf==0) return "the contact term";
    return "boson "+ToString(kf)+" ("+ToString(Flavour(std::abs(kf),kf<0))+")";
  }

  bool HasKey(Scoped_Settings s,const std::string & k)
  {
    const std::vector<std::string> keys = s.GetKeys();
    return std::find(keys.begin(),keys.end(),k)!=keys.end();
  }

  // A fermion line has two legs of the same fermion, differing at most by
  // particle/antiparticle, and couplings are quoted for the particle. A
  // FLAVOUR-CHANGING current (D0 -> K-) has no such thing; the neutral-boson
  // defaults below read this flavour, the W default (V-A) does not need it.
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
      // VA_0_FF builds i gamma^mu (v + a gamma_5), as HADRONS++ has always done, 
      //  so the SM vertex gamma ^ mu(c_V - c_A gamma_5) with c_A = T3 
      //  needs a = -T3.
      v = f.IsoWeak() - 2. * f.Charge() * s2;
      a = -f.IsoWeak();
      return;
    }
    // W on a lepton line only: V-A. Quark and hadron lines need a CKM element
   // and/or form factors, so they must state their couplings.
   if (bos.Kfcode() == kf_Wplus && f.IsLepton())
    {
      v = 1;
      a = -1;
      return;
    }
    THROW(fatal_error,"No default couplings known for exchanged boson "
          +ToString(bos)+" on a "+ToString(f)+" line (defaults exist for the photon, the Z, "
          "and the W on lepton lines only). "
          "Give v and a on both currents, and Coupling on the exchange, "
          "explicitly.");
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
  ME_Base(flavs,indices,name), m_photonalpha(0.), m_nhel(0), m_nin(nin),
  m_chan(0)
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
  // Every current gets the emission coupling, photon or not: it costs
  // nothing, and a current that later radiates cannot then find it unset.
  c->Photon().SetAlpha(m_photonalpha);
  return c;
}

void Exchange_ME::AddExchange(const Diagram & d)
{
  // No Boson means a CONTACT term: D = 1, and the coupling defaults to
  // G_F/sqrt(2). That is exactly what Current_ME always did, so a decay
  // channel carrying no Boson reproduces its old amplitude bit for bit.
  Scoped_Settings e = d.m_e;
  const int kf = d.m_kf;
  const bool contact(kf==0);
  const Flavour bos(std::abs(kf),kf<0);

  Scoped_Settings j1 = e["J1"], j2 = e["J2"];
  const vector<int> & i1 = d.m_i1, & i2 = d.m_i2;

  // Couplings: defaulted from the model where it knows them (photon, Z, and
  // a W on a lepton line; see DefaultVA), and overridden by anything the
  // card states. Asked of the LINE, not the current's legs: LineFlavour
  // returns the first non-antiparticle, and a photon is not an antiparticle,
  // so a current carrying a photon would otherwise be asked for the
  // photon's couplings.
  double v1(0.),a1(0.),v2(0.),a2(0.),g(1.);
  const bool wantdefaults = !contact &&
    !(HasKey(j1,"v") && HasKey(j1,"a") && HasKey(j2,"v") && HasKey(j2,"a"));
  // Per line: a line stating both v and a needs no default, so a hadron
  // line with explicit couplings next to a defaulted lepton line is fine.
  if (wantdefaults) {
    if (!(HasKey(j1,"v") && HasKey(j1,"a")))
      DefaultVA(bos,LineFlavour(m_flavs,m_line1),v1,a1);
    if (!(HasKey(j2,"v") && HasKey(j2,"a")))
      DefaultVA(bos,LineFlavour(m_flavs,m_line2),v2,a2);
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
  x.m_i1 = i1;
  x.m_i2 = i2;
  for (size_t j(0);j<i1.size();++j)
    x.m_qsign.push_back(IsIncomingLeg(i1[j]) ? -1. : +1.);

  CheckFlavourFlow(contact,bos,m_line1,m_line2);

  m_x.push_back(x);

  msg_Debugging()<<"Exchange_ME: "<<(contact?Flavour(kf_none):bos)
                 <<" coupling "<<g<<"\n";
}

bool Exchange_ME::IsPhotonLeg(const int i) const
{
  return std::find(m_photons.begin(),m_photons.end(),i)!=m_photons.end();
}

vector<int> Exchange_ME::StripPhotons(const vector<int> & idx) const
{
  vector<int> line;
  for (size_t j(0);j<idx.size();++j)
    if (!IsPhotonLeg(idx[j])) line.push_back(idx[j]);
  return line;
}

bool Exchange_ME::IsChargedLine(const vector<int> & line) const
{
  for (size_t j(0);j<line.size();++j)
    if (m_flavs[line[j]].IntCharge()!=0) return true;
  return false;
}

// The charge a line delivers to its vertex, counting each leg as outgoing:
// +Q for an outgoing leg, -Q for an incoming one. An e+e- pair annihilating
// gives 0; tau- -> nu_tau gives +1, and the l- nubar_l pair it decays into
// gives -1. A photon attached to a set of lines whose flows sum to zero is
// a gauge-invariant subset of diagrams, by charge conservation alone.
double Exchange_ME::ChargeFlow(const vector<int> & line) const
{
  double c(0.);
  for (size_t j(0);j<line.size();++j)
    c += (IsIncomingLeg(line[j]) ? -1. : 1.)*m_flavs[line[j]].Charge();
  return c;
}

vector<Exchange_ME::Diagram> Exchange_ME::ReadDiagrams(Scoped_Settings & s)
{
  // Every external photon, whether or not the block mentions it: one that
  // nobody places is caught by CheckDiagrams as an uncovered leg.
  m_photons.clear();
  for (size_t i(0);i<m_flavs.size();++i)
    if (m_flavs[i].Kfcode()==kf_photon) m_photons.push_back(i);

  vector<int> expand;
  if (HasKey(s,"Photon")) {
    expand = s["Photon"].SetDefault(vector<int>()).GetVector<int>();
    if (expand.empty())
      THROW(fatal_error,"Photon: is given but names no leg.");
    for (size_t p(0);p<expand.size();++p) {
      if (expand[p]<0 || size_t(expand[p])>=m_flavs.size())
        THROW(fatal_error,"Photon: leg "+ToString(expand[p])+" is out of "
              "range.");
      if (!IsPhotonLeg(expand[p]))
        THROW(fatal_error,"Photon: names leg "+ToString(expand[p])+", which "
              "is "+ToString(m_flavs[expand[p]])+", not a photon.");
      if (std::count(expand.begin(),expand.end(),expand[p])>1)
        THROW(fatal_error,"Photon: names leg "+ToString(expand[p])+" twice.");
    }
  }

  Scoped_Settings ex = s["Exchange"];
  const size_t n = ex.GetItemsCount();
  if (n==0) THROW(fatal_error,"A Currents: block needs an Exchange: list.");

  vector<Diagram> ds;
  for (size_t k(0);k<n;++k) {
    Diagram d(ex[k]);
    d.m_kf = d.m_e["Boson"].GetScalarWithOtherDefault<int>(0);
    Scoped_Settings j1 = d.m_e["J1"], j2 = d.m_e["J2"];
    d.m_i1 = j1["Indices"].SetDefault(vector<int>()).GetVector<int>();
    d.m_i2 = j2["Indices"].SetDefault(vector<int>()).GetVector<int>();
    if (d.m_i1.empty() || d.m_i2.empty())
      THROW(fatal_error,"Both J1 and J2 need Indices.");
    for (size_t j(0);j<d.m_i1.size();++j)
      if (d.m_i1[j]<0 || size_t(d.m_i1[j])>=m_flavs.size())
        THROW(fatal_error,"J1 index out of range.");
    for (size_t j(0);j<d.m_i2.size();++j)
      if (d.m_i2[j]<0 || size_t(d.m_i2[j])>=m_flavs.size())
        THROW(fatal_error,"J2 index out of range.");

    if (expand.empty()) { ds.push_back(d); continue; }

    // Photon: expansion. The entry is the photon-less diagram; each photon
    // goes onto every CHARGED line in turn -- a neutral one cannot emit --
    // and several photons give every combination. The leg is appended, so
    // it is always the last index of the current it lands in.
    const vector<int> l1 = StripPhotons(d.m_i1), l2 = StripPhotons(d.m_i2);
    vector<Diagram> placed(1,d);
    for (size_t p(0);p<expand.size();++p) {
      const int g = expand[p];
      if (std::find(d.m_i1.begin(),d.m_i1.end(),g)!=d.m_i1.end() ||
          std::find(d.m_i2.begin(),d.m_i2.end(),g)!=d.m_i2.end())
        THROW(fatal_error,"Exchange entry "+ToString(k)+" places photon leg "
              +ToString(g)+" explicitly, and Photon: asks for it to be "
              "placed automatically. Use one or the other: with Photon:, "
              "write the entry without the photon.");
      if (!IsChargedLine(l1) && !IsChargedLine(l2))
        THROW(fatal_error,"Photon: asks for leg "+ToString(g)+" to be "
              "attached to exchange entry "+ToString(k)+", but neither of "
              "its lines carries charge.");
      vector<Diagram> next;
      for (size_t c(0);c<placed.size();++c) {
        if (IsChargedLine(l1)) {
          next.push_back(placed[c]); next.back().m_i1.push_back(g);
        }
        if (IsChargedLine(l2)) {
          next.push_back(placed[c]); next.back().m_i2.push_back(g);
        }
      }
      placed.swap(next);
    }
    msg_Info()<<"Exchange_ME: entry "<<k<<" ("<<BosonName(d.m_kf)
              <<") expanded by Photon: into "<<placed.size()
              <<" diagram"<<(placed.size()==1?"":"s")<<".\n";
    ds.insert(ds.end(),placed.begin(),placed.end());
  }
  return ds;
}

// What can be checked from the legs alone, before any current exists. Each
// of these was previously either unchecked, or covered by the blunt rule
// that every entry use identical J1/J2 Indices -- which a photon breaks.
//
//   1. Every diagram covers every external leg exactly once. A leg left out
//      still has its spin states enumerated by Calculate, so the rate would
//      come out wrong by that leg's multiplicity, silently.
//   2. Every diagram joins the SAME two lines once photons are set aside.
//      This is the physical content of the old rule, and what lets a caller
//      assign colour once for the whole amplitude.
//   3. A photon sits on a line that can emit it.
//   4. No diagram is listed twice: easy to do by mixing an explicit entry
//      with Photon:, and it would double that diagram's amplitude.
void Exchange_ME::CheckDiagrams(const vector<Diagram> & ds)
{
  for (size_t k(0);k<ds.size();++k) {
    const Diagram & d = ds[k];
    vector<int> count(m_flavs.size(),0);
    for (size_t j(0);j<d.m_i1.size();++j) ++count[d.m_i1[j]];
    for (size_t j(0);j<d.m_i2.size();++j) ++count[d.m_i2[j]];
    for (size_t i(0);i<count.size();++i) {
      if (count[i]==1) continue;
      if (count[i]==0 && IsPhotonLeg(i))
        THROW(fatal_error,"Photon leg "+ToString(i)+" is not in any current "
              "of exchange entry "+ToString(k)+". Either list it in the "
              "Indices of the current it is emitted from, or give "
              "Photon: "+ToString(i)+" at the top of the Currents: block "
              "to have it attached to every charged line.");
      THROW(fatal_error,"Exchange entry "+ToString(k)+" has leg "
            +ToString(i)+" ("+ToString(m_flavs[i])+") "
            +(count[i]==0 ? string("in neither current")
                          : string("more than once"))+". Each diagram must "
            "hold every external leg exactly once, between J1 and J2.");
    }

    const vector<int> l1 = StripPhotons(d.m_i1), l2 = StripPhotons(d.m_i2);
    if (l1.empty() || l2.empty())
      THROW(fatal_error,"Exchange entry "+ToString(k)+" has a current made "
            "of nothing but photons.");
    if (k==0) { m_line1 = l1; m_line2 = l2; }
    else if (l1!=m_line1 || l2!=m_line2)
      THROW(fatal_error,"All Exchange entries must join the same two lines: "
            "J1 and J2 may differ only by where a photon is attached. Entry "
            +ToString(k)+" pairs the legs differently from entry 0.");

    for (int w(0);w<2;++w) {
      const vector<int> & idx = w ? d.m_i2 : d.m_i1;
      const vector<int> & line = w ? l2 : l1;
      if (idx.size()!=line.size() && !IsChargedLine(line))
        THROW(fatal_error,"Exchange entry "+ToString(k)+" attaches a photon "
              "to J"+ToString(w+1)+", whose legs carry no charge.");
    }

    for (size_t j(0);j<k;++j)
      if (ds[j].m_kf==d.m_kf && ds[j].m_i1==d.m_i1 && ds[j].m_i2==d.m_i2)
        THROW(fatal_error,"Exchange entries "+ToString(j)+" and "+ToString(k)
              +" are the same diagram. If Photon: is given, the entries "
              "should not also place the photon by hand.");
  }
  CheckPhotonCoverage(ds);
}

// Is the photon on every line it should be? Asked separately for each
// exchanged boson, since each boson's diagrams are a set of their own.
//
// Missing placements are not automatically wrong. The deciding quantity is
// the charge flow into the vertex from the lines the photon IS attached to:
//
//   attached to every charged line     complete                 silent
//   a subset with zero net flow        gauge invariant but      warn
//                                      incomplete (ISR alone)
//   a subset with nonzero flow         k.M != 0: depends on     throw
//                                      the photon's gauge
//
// With a neutral boson each line is neutral on its own, so leaving out
// initial- or final-state radiation is a legitimate approximation and only
// warns. With a contact term standing in for a W, tau- -> nu_tau carries +1
// and its l- nubar_l pair -1; a photon on only one of them has no meaning.
// Photon: never reaches either of the last two rows.
void Exchange_ME::CheckPhotonCoverage(const vector<Diagram> & ds) const
{
  if (m_photons.empty()) return;

  // A CHARGED boson carries charge between the two lines, so gauge
  // invariance also needs the photon emitted from the boson itself -- the
  // WW gamma vertex, which nothing here provides. The coverage rule below
  // cannot see the gap: the two lines' charge flows (+1, -1) cancel, so it
  // would call the photon complete. Without the W emission the result
  // depends on the photon's gauge at O(k.q/M_W^2). A contact term has no
  // boson to radiate from and is exactly gauge invariant, so use that.
  for (size_t k(0);k<ds.size();++k) {
    if (ds[k].m_kf==0) continue;
    if (Flavour(std::abs(ds[k].m_kf),ds[k].m_kf<0).IntCharge()==0) continue;
    THROW(not_implemented,"Exchange entry "+ToString(k)+" exchanges "
          +BosonName(ds[k].m_kf)+", which is charged, in a process with a "
          "photon. Emission from the boson itself (the WW gamma vertex) is "
          "not implemented, and without it the amplitude is not gauge "
          "invariant. Drop Boson: to use the contact (Fermi) term, which is "
          "exactly gauge invariant with the photon on the fermion lines.");
  }
  vector<int> kfs;
  for (size_t k(0);k<ds.size();++k)
    if (std::find(kfs.begin(),kfs.end(),ds[k].m_kf)==kfs.end())
      kfs.push_back(ds[k].m_kf);

  const bool charged[2] = { IsChargedLine(m_line1), IsChargedLine(m_line2) };
  for (size_t b(0);b<kfs.size();++b) {
    const string who = BosonName(kfs[b]);
    for (size_t p(0);p<m_photons.size();++p) {
      const int g = m_photons[p];
      bool on[2] = { false, false };
      for (size_t k(0);k<ds.size();++k) {
        if (ds[k].m_kf!=kfs[b]) continue;
        const vector<int> & i1 = ds[k].m_i1;
        on[std::find(i1.begin(),i1.end(),g)!=i1.end() ? 0 : 1] = true;
      }
      if (on[0]==charged[0] && on[1]==charged[1]) continue;

      const int bare = on[0] ? 2 : 1;
      const vector<int> & bareline = bare==1 ? m_line1 : m_line2;
      double flow(0.);
      if (on[0]) flow += ChargeFlow(m_line1);
      if (on[1]) flow += ChargeFlow(m_line2);
      if (std::abs(flow)>1.e-9)
        THROW(fatal_error,"For "+who+", photon leg "+ToString(g)+" is never "
              "attached to J"+ToString(bare)+" (legs "+ToString(bareline)
              +"), which is charged. The lines it IS attached to carry net "
              "charge "+ToString(flow)+" into the vertex, so these diagrams "
              "are not gauge invariant: the result would depend on the "
              "photon's gauge. Add the placements on J"+ToString(bare)
              +", or give Photon: "+ToString(g)+" to have every one "
              "generated.");
      msg_Info()<<"Warning in "<<METHOD<<": for "<<who<<", photon leg "<<g
                <<" is never attached to J"<<bare<<" (legs "<<bareline
                <<"), which is charged. The diagrams given are gauge "
                <<"invariant on their own, so the result is physical, but "
                <<"it is INCOMPLETE -- emission from that line and its "
                <<"interference are missing. Fine if intended (radiation "
                <<"from one side only); otherwise add the placements on J"
                <<bare<<", or give Photon: "<<g<<".\n";
    }
  }
}

void Exchange_ME::SetPhotonGaugeTest(const bool on)
{
  for (size_t k(0);k<m_x.size();++k) {
    m_x[k].p_j1->Photon().SetGaugeTest(on);
    m_x[k].p_j2->Photon().SetGaugeTest(on);
  }
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
  // Resolve and check every diagram BEFORE building any current: a block
  // with a missing leg or a misplaced photon should be refused for that
  // reason, not for whatever a current makes of the legs it was handed.
  m_photonalpha = MODEL::s_model->ScalarConstant("alpha_QED");
  if (HasKey(s,"Photon_Alpha"))
    m_photonalpha = s["Photon_Alpha"].SetDefault(m_photonalpha).Get<double>();
  const vector<Diagram> ds = ReadDiagrams(s);
  CheckDiagrams(ds);
  for (size_t k(0);k<ds.size();++k) AddExchange(ds[k]);

  // You have to specify why the amplitude vanishes when using this flag, otherwise its not very helpful.
  bool allzero(true);
  for (size_t k(0);k<m_x.size();++k)
    if (!m_x[k].p_j1->Vanishes() && !m_x[k].p_j2->Vanishes()) allzero = false;
  if (allzero)
    THROW(fatal_error,"Every exchange in this Currents: block vanishes "
          "identically, so the amplitude is zero "
          "everywhere and the integrator never finds ME != 0.");

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
  // and would otherwise be classified as though it were one. A process
  // with a photon has five legs or more, so never gets here, and every
  // exchange of a 2 -> 2 has the same J1 -- CheckDiagrams saw to that.
  m_chan = 0;
  if (m_nin==2 && m_flavs.size()==4 && m_line1.size()==2) {
    const int a = std::min(m_line1[0],m_line1[1]);
    const int b = std::max(m_line1[0],m_line1[1]);
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
//
// A photon listed in J1 is outgoing, so it enters with +1 and q comes out as
// k - p0 - p1 for initial-state radiation: (p0+p1-k)^2, as it should be. Only
// q^2 and the product q^mu q^nu are used, so the overall sign is immaterial.
Vec4D Exchange_ME::Q(const Exchange & x,const Vec4D_Vector & moms) const
{
  Vec4D q(0.,0.,0.,0.);
  for (size_t j(0);j<x.m_i1.size();++j) q += x.m_qsign[j]*moms[x.m_i1[j]];
  return q;
}

void Exchange_ME::Calculate(const Vec4D_Vector & moms, bool anti)
{
  // Each exchange's q and propagator, once per point rather than once per
  // helicity. A propagator evaluated exactly on its pole -- in practice a
  // photon at q^2 = 0, which only an unlucky phase-space point reaches --
  // zeroes the whole amplitude, as the scattering side has always done.
  vector<Vec4D>   q(m_x.size());
  vector<Complex> prop(m_x.size(),Complex(1.,0.));
  for (size_t k(0);k<m_x.size();++k) {
    q[k] = Q(m_x[k],moms);
    if (m_x[k].p_prop==NULL) continue;
    const double q2 = q[k].Abs2();
    if (!(std::abs(q2)>0.)) {
      for (size_t n(0);n<size();++n) (*this)[n] = Complex(0.,0.);
      return;
    }
    prop[k] = (*m_x[k].p_prop)(q2);
  }

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
      if (x.m_M2>0.) lor -= (J1*q[k])*(J2*q[k])/x.m_M2;
      amp += x.m_g*prop[k]*lor;
    }
    // The AMPLITUDE is kept, phase and all. A decay contracts it with the
    // production density matrix; a cross section squares it. Neither belongs
    // here.
    (*this)[n] = amp;
  }
}
