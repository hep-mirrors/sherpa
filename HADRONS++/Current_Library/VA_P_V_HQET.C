#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class HQET : public FormFactor_Base {
// hep-ph/
  double m_R1;
  double m_R2;
  double m_RStar;
  double m_rho2;
public:
  HQET(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

HQET::HQET(GeneralModel model, double* masses, const Flavour_Vector& flavs,
           const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_rho2  = model("HQET_rho2",1.179);
  m_R1    = model("HQET_R1",1.417);
  m_R2    = model("HQET_R2",0.836);
  m_RStar = ( 2.0*sqrt(m_m0*m_m1))/(m_m0+m_m1);
}

void HQET::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  double w  = (m_m0*m_m0 + m_m1*m_m1 - q2) / (2.0 * m_m0 * m_m1);
  double xi = 1.0 - m_rho2*(w-1.0);

  m_A0 = 0.0;
  m_A1 = (1.0 - q2/(sqr(m_m0+m_m1)))*xi/m_RStar;
  m_A2 = m_R2/m_RStar*xi;
  m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
  m_V  = m_R1/m_RStar*xi;
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
