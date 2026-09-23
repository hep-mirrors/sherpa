#include "METOOLS/HadronCurrents/VA_P_0.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

void VA_P_0::SetModelParameters(struct GeneralModel model)
{
  if (p_i.size()!=1)
    THROW(fatal_error,"Current "+m_name+" takes exactly one index, the "
          "decaying pseudoscalar, but was given "+ToString(p_i.size())+".");
  if (p_i[0]!=0)
    msg_Info()<<"Warning in "<<METHOD<<": "<<m_name<<" was given index "
              <<p_i[0]<<" rather than 0. This current annihilates the "
              <<"DECAYING meson, which is external index 0.\n";
  if (m_name=="VA_0_P")
    msg_Info()<<"Warning in "<<METHOD<<": the tag VA_0_P reads as "
              <<"'nothing in, one pseudoscalar out' under the <in>_<out> "
              <<"convention of this directory, but the pseudoscalar here "
              <<"is the DECAYING particle. Use VA_P_0; the result is "
              <<"identical.\n";

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
  msg_Tracking()<<"### "<<m_name<<": "<<m_flavs[p_i[0]]<<" -> leptons, "
                <<"f_P = "<<m_fP<<" GeV, V = "<<m_Vxx<<"\n";
}

void VA_P_0::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // J^mu = i f_P V p^mu. The whole current is proportional to the
  // decaying meson's momentum, which is why the rate is helicity
  // suppressed: contracted with the lepton current it survives only
  // through the charged-lepton mass. That is why D_s -> tau nu is
  // ten times D_s -> mu nu despite the far smaller phase space.
  const Vec4C J = Complex(0.,m_fP*m_Vxx)*Vec4C(moms[p_i[0]]);
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
  st<<"Deprecated alias for {\\tt VA\\_P\\_0} - the pseudoscalar is the \n"
    <<"decaying particle, so P$\\to$0 is the correct reading. \n"
    <<std::endl;
}
