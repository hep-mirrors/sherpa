#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class HQET2 : public FormFactor_Base {
  // hep-ph/9712417
  double m_rho2;
  double m_V1_1;
public:
  HQET2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

HQET2::HQET2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices) :
  FormFactor_Base(model, masses, flavs, indices)
{
  m_rho2  = model("HQET2_rho2",1.19);
  m_V1_1  = model("HQET2_V1_1",0.98);
}

void HQET2::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double w = (p0/m_m0)*(p1/m_m1);
  const double z = (sqrt(w+1.0)-sqrt(2.0))/(sqrt(w+1.0)+sqrt(2.0));
  double V1 = m_V1_1*(1.0-8.0*m_rho2*z+(51.0*m_rho2-10.0)*z*z-(252.0*m_rho2-84.0)*z*z*z);
  m_fplus  = V1;
  m_f0     = 0.0;
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
