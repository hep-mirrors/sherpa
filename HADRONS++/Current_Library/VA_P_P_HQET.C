#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class HQET : public FormFactor_Base {
// hep-ph/
  double m_rho2;
  double m_c;
public:
  HQET(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

HQET::HQET(GeneralModel model, double* masses, const Flavour_Vector& flavs,
           std::vector<int>& indices) :
  FormFactor_Base(model, masses, flavs, indices)
{
  m_rho2  = model("HQET_rho2",0.7);
  m_c     = model("HQET_R1",0.0);
}

void HQET::CalcFFs( Vec4D p0, Vec4D p1 )
{
  Vec4D v0 = p0/m_m0;
  Vec4D v1 = p1/m_m1;
  double w = v0*v1;
  double R = 2.*sqrt(m_m0*m_m1)/(m_m0+m_m1);
  m_fplus  = 1.0-m_rho2*(w-1.0)+m_c*(w-1.0)*(w-1.0)/R;
  m_f0     = 0.0;
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
