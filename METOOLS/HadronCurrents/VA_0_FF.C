#include "METOOLS/HadronCurrents/VA_0_FF.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_0_FF::VA_0_FF(const ATOOLS::Flavour_Vector& flavs,
                 const std::vector<int>& indices,
                 const std::string& name) :
  Current_Base(flavs,indices,name),
  m_cR(0.,0.), m_cL(1.,0.)
{
  // No topology check here. A current cannot know which of its legs are
  // incoming until it is told how many external legs are incoming at all,
  // and that is context the constructor does not have: leg 0 is the
  // decaying particle in a decay, but the FIRST BEAM in a scattering
  // process. Exchange_ME calls SetNIn, which then calls CheckTopology.
}

void VA_0_FF::CheckTopology() const {
  if (p_i.size()!=2 && p_i.size()!=3)
    THROW(fatal_error,"Current "+m_name+" needs two legs, or three with a "
          "photon, got "+ToString(p_i.size())+".");
  if (p_i.size()==3) {
    const Flavour & g = m_flavs[p_i[2]];
    if (g.Kfcode()!=kf_photon)
      THROW(fatal_error,"Current "+m_name+": a third leg must be a photon, "
            "and must come LAST in Indices; leg "+ToString(p_i[2])+" is "
            +ToString(g)+".");
    if (IsIncoming(p_i[2]))
      THROW(not_implemented,"Current "+m_name+": an incoming photon is not "
            "implemented, only an emitted one.");
    if (m_flavs[p_i[0]].IntCharge()==0 && m_flavs[p_i[1]].IntCharge()==0)
      THROW(fatal_error,"Current "+m_name+" was given a photon, but neither "
            "of its fermion legs is charged.");
  }

  // IS INDEX 0 THE BARRED SPINOR? ubar for an outgoing particle, vbar for an
  // incoming antiparticle; index 1 is then u (incoming particle) or v
  // (outgoing antiparticle). XYZFunc picks u/v but never reorders, so a
  // swapped pair silently gives the charge-conjugate line.
  for (size_t k(0);k<2;++k) {
    const bool barred = (IsIncoming(p_i[k])==m_flavs[p_i[k]].IsAnti());
    if (barred!=(k==0))
      THROW(fatal_error,"Current "+m_name+": leg "+ToString(p_i[k])+" ("
            +ToString(m_flavs[p_i[k]])+") is in Indices slot "+ToString(k)
            +", but slot 0 must be the barred spinor (outgoing particle or "
            "incoming antiparticle). Swap the first two Indices or check your logic.");
      // If this is too restrictive in the tau decays for an example I haven't thought of, feel free to make it an error message not a "stop the whole run".
      // msg_Error()<< "Current " + m_name + ": leg " + ToString(p_i[k]) + " (" + ToString(m_flavs[p_i[k]]) + ") is in Indices slot " + ToString(k) + ", but slot 0 must be the barred spinor (outgoing particle or "
      //                                                                                                                                                 "incoming antiparticle). Swap the first two Indices or check your logic.\n";
  }

  // DOES THE TAG MATCH THE TOPOLOGY?
  //
  // Measured, not assumed. The old test asked whether the index list
  // contained external index 0 -- true for the decaying particle, but also
  // true for the first BEAM of a scattering process, where it means
  // nothing. Counting incoming legs works in both contexts.
  const size_t nin = NIncomingLegs();
  const std::string got = (nin==2) ? "VA_FF_0" : (nin==1 ? "VA_F_F"
                                                         : "VA_0_FF");
  const std::string what = (nin==2) ? "both legs incoming (an annihilating "
                                      "pair)"
                         : (nin==1 ? "one leg incoming and one outgoing (a "
                                     "fermion line passing through)"
                                   : "both legs outgoing (a produced pair)");
  if (m_name!=got)
    msg_Info()<<"Warning in "<<METHOD<<": current tagged "<<m_name
              <<" was given "<<what<<", which by this directory's "
              <<"<in>_<out> naming is "<<got<<". The result is unaffected "
              <<"-- one class serves all three, and XYZFunc resolves the "
              <<"spinor types from the flavours -- but the tag is wrong "
              <<"and should be corrected.\n";
}

void VA_0_FF::SetModelParameters(struct GeneralModel model)
{
  // Unchanged from the legacy implementation: chiral couplings i(v +- a)
  // with defaults v=1, a=-1, i.e. the pure left-handed Standard-Model
  // current. The global i is a convention shared with the hadronic
  // currents and cancels in every rate.
  m_cR = Complex(0.,model("v",1.)+model("a",-1.));
  m_cL = Complex(0.,model("v",1.)-model("a",-1.));
  
  if (Vanishes())
    msg_Error()<<"Warning in "<<METHOD<<": current "<<m_name<<" on "
               <<m_flavs[p_i[0]]<<" has v = a = 0 (e.g. a photon coupled to "
               <<"a neutral lepton), so it vanishes identically. Check process.\n";

  const int ff = int(model("V_A_FORM_FACTOR",1.)+0.5);
  if (ff!=1)
    THROW(fatal_error,"V_A_FORM_FACTOR = "+ToString(ff)+" requested for "
          "current "+m_name+", but the only implemented option is 1 "
          "(no form factor).");
  msg_Tracking()<<"### "<<m_name<<": V-A fermion current, v = "
                <<model("v",1.)<<", a = "<<model("a",-1.)
                <<", legs incoming = "<<NIncomingLegs()
                <<"\n";
}

void VA_0_FF::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  if (p_i.size()==3) { CalcWithPhoton(moms,anti); return; }
  XYZFunc F(moms,m_flavs,anti,p_i);
  for (int h0=0;h0<2;h0++) {
    for (int h1=0;h1<2;h1++) {
      // Index 0 must be the barred spinor (for the non-anti case),
      // index 1 the unbarred one.
      const Vec4C amp = F.L(0,h0, 1,h1, m_cR,m_cL);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(0,h0));
      spins.push_back(make_pair(1,h1));
      Insert(amp,spins);
    }
  }
}

void VA_0_FF::CalcWithPhoton(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // Under anti, XYZFunc flips every spinor type and the photon's
  // polarisation would have to follow; nothing here has been checked for
  // that, and only the decay path asks for it.
  if (anti)
    THROW(not_implemented,"Current "+m_name+": photon emission under the "
          "charge-conjugated evaluation (anti) is not implemented.");

  const Vec4D p0 = moms[p_i[0]], p1 = moms[p_i[1]];
  const Vec4D k  = m_photon.LightLike(moms[p_i[2]]);

  // The pseudo-leg for the spin sum: k appended as an extra momentum, with
  // a massless FERMION flavour so that XYZFunc builds u(k,s). Which fermion
  // is immaterial -- XYZFunc takes the mass from sqrt(k^2) = 0 and the
  // spinor type from IsAnti -- so the line's own particle is used, as
  // PHOTONS++ does for its pseudo-particles. The real legs keep their own
  // momenta and flavours, so their spinors are bit for bit the ones the
  // photon-less current builds.
  Vec4D_Vector   ext(moms);
  Flavour_Vector efl(m_flavs);
  Flavour pseudo = m_flavs[p_i[1]];
  if (pseudo.IsAnti()) pseudo = pseudo.Bar();
  ext.push_back(k);
  efl.push_back(pseudo);
  const vector<int> idx = { p_i[0], p_i[1], int(moms.size()) };
  XYZFunc F(ext,efl,anti,idx);

  const Complex one(1.,0.);
  for (int lam=0;lam<2;lam++) {
    const Vec4C eps = m_photon.EpsStar(k,lam);
    for (int h0=0;h0<2;h0++) {
      for (int h1=0;h1<2;h1++) {
        const Vec4C J0 = F.L(0,h0, 1,h1, m_cR,m_cL);
        // S_0 = ubar_0 eps* k Gamma u_1,  S_1 = ubar_0 Gamma k eps* u_1,
        // with k-slash = sum_s u(k,s) ubar(k,s) on pseudo-leg 2.
        Vec4C S0(0.,0.,0.,0.), S1(0.,0.,0.,0.);
        for (int s=0;s<2;s++) {
          S0 += F.X(0,h0, eps, 2,s, one,one)*F.L(2,s, 1,h1, m_cR,m_cL);
          S1 += F.L(0,h0, 2,s, m_cR,m_cL)*F.X(2,s, eps, 1,h1, one,one);
        }
        const Vec4C J =
          m_photon.Emit(m_flavs[p_i[0]],p0,IsIncoming(p_i[0]),k,eps,J0,S0)
        + m_photon.Emit(m_flavs[p_i[1]],p1,IsIncoming(p_i[1]),k,eps,J0,S1);
        vector<pair<int,int> > spins;
        spins.push_back(make_pair(0,h0));
        spins.push_back(make_pair(1,h1));
        spins.push_back(make_pair(2,lam));
        Insert(J,spins);
      }
    }
  }
}

// Preferred tag: both fermions outgoing, e.g. the lepton pair of a
// semileptonic decay.
DEFINE_CURRENT_GETTER(METOOLS::VA_0_FF,"VA_0_FF")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_0_FF>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: the $\\ell\\nu_\\ell$ pair of $D\\rightarrow K\\ell\\nu$ \n\n"
    <<"Order: 0 = bar'ed spinor, 1 = non-bar'ed spinor \n\n"
    <<"\\[\\bar{u}(p_0) \\gamma_\\mu [ v-a\\gamma_5 ] u(p_1) \\] \n\n"
    <<"Parameters: {\\tt v} (default 1), {\\tt a} (default -1). \n\n"
    <<"Both legs must be outgoing. For a fermion transition current \n"
    <<"($\\tau\\rightarrow\\nu_\\tau$, $b\\rightarrow c$) use \n"
    <<"{\\tt VA\\_F\\_F} instead - same class, checked topology. \n"
    <<std::endl;
}

// Transition tag: one leg incoming, one outgoing. A first-class name, not
// a deprecation - tau -> nu_tau and a t-channel lepton line both really
// are F -> F currents and read correctly as VA_F_F.
DEFINE_CURRENT_GETTER(METOOLS::VA_F_F,"VA_F_F")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_F_F>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\tau^-\\rightarrow\\nu_\\tau$, or the electron line \n"
    <<"of $e^-\\mu^-\\rightarrow e^-\\mu^-$ \n\n"
    <<"Order: 0 = bar'ed spinor, 1 = non-bar'ed spinor \n\n"
    <<"Identical to {\\tt VA\\_0\\_FF}; use this tag when exactly ONE \n"
    <<"leg is incoming. \n"
    <<std::endl;
}

// Annihilation tag: both legs incoming. Only reachable in a scattering
// process, where there is no decaying particle and leg 0 is a beam.
DEFINE_CURRENT_GETTER(METOOLS::VA_FF_0,"VA_FF_0")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_FF_0>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: the $e^+e^-$ side of $e^+e^-\\rightarrow\\mu^+\\mu^-$ \n\n"
    <<"Order: 0 = bar'ed spinor, 1 = non-bar'ed spinor \n\n"
    <<"Identical to {\\tt VA\\_0\\_FF}; use this tag when BOTH legs \n"
    <<"are incoming. \n"
    <<std::endl;
}
