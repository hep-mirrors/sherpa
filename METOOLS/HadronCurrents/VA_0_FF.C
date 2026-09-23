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
  m_cR(0.,0.), m_cL(1.,0.), m_hasincoming(false)
{
  // External index 0 is the decaying particle, so its presence in the
  // index list is exactly what distinguishes F -> F from 0 -> F Fbar.
  for (size_t k=0;k<p_i.size();k++) if (p_i[k]==0) m_hasincoming = true;
  CheckTopology();
}

void VA_0_FF::CheckTopology() const {
  if (m_name=="VA_0_FF" && m_hasincoming)
    THROW(fatal_error,"Current VA_0_FF was given the decaying particle "
          "(external index 0) as one of its two legs. VA_0_FF is the "
          "'nothing in, two fermions out' current; for a fermion "
          "transition such as tau -> nu_tau use the tag VA_F_F.");
  if (m_name=="VA_F_F" && !m_hasincoming)
    msg_Info()<<"Warning in "<<METHOD<<": current tagged VA_F_F, but "
              <<"neither leg is the decaying particle - both fermions are "
              <<"outgoing. By this directory's <in>_<out> naming that is "
              <<"a 0 -> F Fbar current and should be tagged VA_0_FF. The "
              <<"result is unaffected (same class, XYZFunc resolves the "
              <<"spinor types), but please relabel the channel.\n";
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
                <<", topology = "<<(m_hasincoming?"F -> F":"0 -> F Fbar")
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

// Legacy/transition tag: one leg is the decaying fermion. Kept as a
// first-class name, not a deprecation - tau -> nu_tau really is an
// F -> F current and reads correctly as VA_F_F.
DEFINE_CURRENT_GETTER(METOOLS::VA_F_F,"VA_F_F")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_F_F>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\tau^-\\rightarrow\\nu_\\tau$ transition current \n\n"
    <<"Order: 0 = bar'ed spinor, 1 = non-bar'ed spinor \n\n"
    <<"Identical to {\\tt VA\\_0\\_FF}; use this tag when ONE leg is \n"
    <<"the decaying particle. If both are outgoing, use \n"
    <<"{\\tt VA\\_0\\_FF} - a warning is issued otherwise. \n"
    <<std::endl;
}
