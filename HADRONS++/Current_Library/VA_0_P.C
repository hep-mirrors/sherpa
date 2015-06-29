#include "HADRONS++/Current_Library/VA_0_P.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;

void VA_0_P::SetModelParameters( struct GeneralModel _md )
{
  double fP(1.0), Vxx(1.0);
  switch(m_flavs[p_i[0]].Kfcode()) {
  case kf_pi_plus:
    fP=0.1307;
    Vxx=Tools::Vud;
    break;
  case kf_K_plus:
    fP=0.1598;
    Vxx=Tools::Vus;
    break;
  case kf_D_plus:
    fP=0.2226;
    Vxx=Tools::Vcd;
    break;
  case kf_D_s_plus:
    fP=0.294;
    Vxx=Tools::Vcs;
    break;
  case kf_B_plus:
    fP=0.176;
    Vxx=Tools::Vub;
    break;
  case kf_B_c:
    fP=0.36;
    Vxx=Tools::Vcb;
    break;
  default:
    fP=1.0;
    Vxx=1.0;
  }
  m_Vxx = _md("Vxx", Vxx);
  m_fP  = _md("fP", fP);
}

void VA_0_P::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  // 0 is Pseudoscalar
  double factor=m_fP*m_Vxx;
  Insert(Vec4C((factor*Complex(0.0,1.0))*moms[p_i[0]]),0);
}

DEFINE_CURRENT_GETTER(VA_0_P,"VA_0_P")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_P>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi $ \n\n"
    <<"Order: 0 = Pseudoscalar \n\n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
