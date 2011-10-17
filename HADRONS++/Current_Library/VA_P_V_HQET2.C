#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class HQET2 : public FormFactor_Base {
  // hep-ph/9712417
  double m_R1_1;
  double m_R2_1;
  double m_hA1_1;
  double m_rho2;
public:
  HQET2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

HQET2::HQET2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_rho2  = model("HQET2_rho2",1.34);
  m_hA1_1 = model("HQET2_hA1_1",0.91);
  m_R1_1  = model("HQET2_R1_1",1.18);
  m_R2_1  = model("HQET2_R2_1",0.71);
}

void HQET2::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2 = (p0-p1).Abs2();
  double w  = (m_m0*m_m0 + m_m1*m_m1 - q2) / (2.0 * m_m0 * m_m1);
  const double z = (sqrt(w+1.0)-sqrt(2.))/(sqrt(w+1.0)+sqrt(2.));
  double hA1 = m_hA1_1*(1.- 8.*m_rho2*z + (53.*m_rho2-15.)*sqr(z) - (231.*m_rho2-91.)*sqr(z)*z);
  double R1 = m_R1_1-0.12*(w-1.0)+0.05*sqr(w-1.0);
  double R2 = m_R2_1+0.11*(w-1.0)-0.06*sqr(w-1.0);
  double Rstar = ( 2.0*sqrt(m_m0*m_m1))/(m_m0+m_m1);

  m_A1 = (1.0 - (q2/((m_m0+m_m1)*(m_m0+m_m1))))*hA1;
  m_A1 = m_A1/Rstar;
  m_A2 = (R2/Rstar)*hA1;
  m_V  = (R1/Rstar)*hA1;
  m_A0 = 0.0;
  m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
