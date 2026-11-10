#include "METOOLS/HadronCurrents/VA_P_0.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

void VA_P_0::CheckTopology() const
{
  if (p_i.size()!=1)
    THROW(fatal_error,"Current "+m_name+" takes exactly one index, the "
          "pseudoscalar, but was given "+ToString(p_i.size())+".");
  // DOES THE TAG MATCH THE TOPOLOGY? VA_P_0 annihilates an incoming meson
  // (the decaying D in D -> l nu), VA_0_P creates an outgoing one (the
  // pi in tau -> pi nu).
  const bool incoming = IsIncoming(p_i[0]);
  if (m_name=="VA_0_P" && incoming)
    THROW(fatal_error,"Current VA_0_P creates its pseudoscalar, but leg "
          +ToString(p_i[0])+" ("+ToString(m_flavs[p_i[0]])+") is incoming. "
          "For a decaying meson use VA_P_0.");
  if (m_name!="VA_0_P" && !incoming)
    THROW(fatal_error,"Current "+m_name+" annihilates its pseudoscalar, but "
          "leg "+ToString(p_i[0])+" ("+ToString(m_flavs[p_i[0]])+") is "
          "outgoing. For a produced meson (tau -> P nu) use VA_0_P.");
}

void VA_P_0::SetModelParameters(struct GeneralModel model)
{
  double fP = 1., Vxx = 1.;
  switch (m_flavs[p_i[0]].Kfcode()) {
  case kf_pi_plus:  fP = 0.1304; Vxx = Tools::Vud; break;
  case kf_K_plus:   fP = 0.1561; Vxx = Tools::Vus; break;
  case kf_D_plus:   fP = 0.2067; Vxx = Tools::Vcd; break;
  case kf_D_s_plus: fP = 0.2600; Vxx = Tools::Vcs; break;
  case kf_B_plus:   fP = 0.1760; Vxx = Tools::Vub; break;
  case kf_B_c:      fP = 0.3600; Vxx = Tools::Vcb; break;
  default:
    msg_Info()<<"Warning in "<<METHOD<<": no decay constant known for "
              <<m_flavs[p_i[0]]<<" in "<<m_name<<"; using f_P = 1, "
              <<"V = 1. Set fP and Vxx explicitly.\n";
    break;
  }
  m_Vxx = model("Vxx",Vxx);
  m_fP  = model("fP", fP);
  msg_Tracking()<<"### "<<m_name<<": "<<m_flavs[p_i[0]]<<", "
                <<"f_P = "<<m_fP<<" GeV, V = "<<m_Vxx<<"\n";
}

void VA_P_0::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // J^mu = i f_P V p^mu, conjugated when the meson is created. The whole
  // current is proportional to the meson's momentum, which is why the
  // rate is helicity suppressed: contracted with the lepton current it
  // survives only through the charged-lepton mass. That is why
  // D_s -> tau nu is ten times D_s -> mu nu despite the far smaller phase
  // space, and why tau -> pi nu scales as m_tau^3 (1-m_pi^2/m_tau^2)^2.
  Vec4C J = Complex(0.,m_fP*m_Vxx)*Vec4C(moms[p_i[0]]);
  if (!IsIncoming(p_i[0])) J = conj(J);
  Insert(anti?conj(J):J,0);
}

DEFINE_CURRENT_GETTER(METOOLS::VA_P_0,"VA_P_0")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_0>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D_s^+ \\rightarrow \\tau^+\\nu_\\tau $ \n\n"
    <<"Order: 0 = the decaying pseudoscalar \n\n"
    <<"\\[ \\langle 0| \\bar q\\gamma^\\mu(1-\\gamma_5)c|P(p)\\rangle "
    <<"= i f_P V_{qq'} p^\\mu \\] \n\n"
    <<"Parameters: {\\tt fP} [GeV], {\\tt Vxx}. Defaults are tabulated \n"
    <<"for $\\pi$, $K$, $D$, $D_s$, $B$, $B_c$. \n"
    <<std::endl;
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_P,"VA_0_P")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_0_P>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ \\tau^-\\rightarrow \\pi^-\\nu_\\tau $ \n\n"
    <<"Order: 0 = the produced pseudoscalar \n\n"
    <<"\\[ \\langle P(p)| \\bar q\\gamma^\\mu(1-\\gamma_5)q'|0\\rangle "
    <<"= -i f_P V_{qq'} p^\\mu \\] \n\n"
    <<"Same class as {\\tt VA\\_P\\_0}, for a created rather than an \n"
    <<"annihilated meson. Parameters: {\\tt fP} [GeV], {\\tt Vxx}. \n"
    <<std::endl;
}
