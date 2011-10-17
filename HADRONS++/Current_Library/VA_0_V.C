#include "HADRONS++/Current_Library/VA_0_V.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "METOOLS/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;

void VA_0_V::SetModelParameters( struct GeneralModel _md )
{
  double fV(1.0), Vxx(1.0);
  switch(m_flavs[p_i[0]].Kfcode()) {
  case kf_rho_770_plus:
    Vxx=Tools::Vud;
    break;
  case kf_K_star_892_plus:
    Vxx=Tools::Vus;
    break;
  case kf_D_star_2010_plus:
    Vxx=Tools::Vcd;
    break;
  case kf_D_s_star_plus:
    Vxx=Tools::Vcs;
    break;
  case kf_B_star_plus:
    Vxx=Tools::Vub;
    break;
  default:
    fV=1.0;
    Vxx=1.0;
  }
  m_Vxx = _md("Vxx", Vxx);
  m_fV  = _md("fV", fV);
}

void VA_0_V::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  // 0 is Vector
  double M = m_flavs[p_i[0]].HadMass();
  double factor = m_fV*M;
  Polarization_Vector eps(moms[p_i[0]], sqr(m_flavs[p_i[0]].HadMass()));
  for(int h=0;h<3;h++) {
    Insert(factor*eps[h], h);
  }
}

DEFINE_CURRENT_GETTER(VA_0_V,VA_0_V_Getter,"VA_0_V")

void VA_0_V_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ J/\\Psi \\rightarrow 0 $ \n\n"
    <<"Order: 0 = Vector \n\n"
    <<"Reference: hep-ph/9503201 (11) \n"
    <<std::endl;
}
