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
  if (p_i.size()!=2)
    THROW(fatal_error,"Current "+m_name+" needs exactly two legs, got "
          +ToString(p_i.size())+".");

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
