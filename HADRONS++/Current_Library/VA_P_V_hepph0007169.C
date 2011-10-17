#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class hepph0007169 : public FormFactor_Base {
  double m_V_0, m_V_mfit2, m_V_delta;
  double m_A1_0, m_A1_mfit2, m_A1_delta;
  double m_A2_0, m_A2_mfit2, m_A2_delta;
  double m_A3_0, m_A3_mfit2, m_A3_delta;
  double Fit(double q2, double F0, double mfit2, double delta);

public:
  hepph0007169(GeneralModel model, double* masses, const Flavour_Vector& flavs,
               const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

  hepph0007169::hepph0007169(GeneralModel model, double* masses,
                             const Flavour_Vector& flavs,
                             const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_V_0=0.0; m_V_mfit2=1.0; m_V_delta=0.0;
  m_A1_0=0.0; m_A1_mfit2=1.0; m_A1_delta=0.0;
  m_A2_0=0.0; m_A2_mfit2=1.0; m_A2_delta=0.0;
  m_A3_0=0.0; m_A3_mfit2=1.0; m_A3_delta=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B_c) {
    if (kf1==kf_B_s_star) {
      m_A1_0 = -0.33;
      m_A1_mfit2 = 1.86*1.86;
      m_A1_delta = -0.13;
      m_A2_0 = 0.40;
      m_A2_mfit2 = 3.44*3.44;
      m_A2_delta = -107;
      m_A3_0 = 10.4;
      m_A3_mfit2 = 1.73*1.73;
      m_A3_delta = -0.09;
      m_V_0 = 3.25;
      m_V_mfit2 = 1.76*1.76;
      m_V_delta = -0.052;
    }
    else if (kf1==kf_J_psi_1S) {
      m_A1_0 = 0.68;
      m_A1_mfit2 = 8.2*8.2;
      m_A1_delta = 1.4;
      m_A2_0 = 0.66;
      m_A2_mfit2 = 5.91*5.91;
      m_A2_delta = 0.052;
      m_A3_0 = -1.13;
      m_A3_mfit2 = 5.67*5.67;
      m_A3_delta = -0.004;
      m_V_0 = 0.96;
      m_V_mfit2 = 5.65*5.65;
      m_V_delta = 0.0013;
    }
  }

  m_V_0     = model("V_0",m_V_0);
  m_V_mfit2 = model("V_mfit2",m_V_mfit2);
  m_V_delta = model("V_delta",m_V_delta);
  m_A1_0     = model("A1_0",m_A1_0);
  m_A1_mfit2 = model("A1_mfit2",m_A1_mfit2);
  m_A1_delta = model("A1_delta",m_A1_delta);
  m_A2_0     = model("A2_0",m_A2_0);
  m_A2_mfit2 = model("A2_mfit2",m_A2_mfit2);
  m_A2_delta = model("A2_delta",m_A2_delta);
  m_A3_0     = model("A3_0",m_A3_0);
  m_A3_mfit2 = model("A3_mfit2",m_A3_mfit2);
  m_A3_delta = model("A3_delta",m_A3_delta);
}

void hepph0007169::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  m_A0 = 0.0;
  m_A1 = Fit(q2, m_A1_0, m_A1_mfit2, m_A1_delta);
  m_A2 = Fit(q2, m_A2_0, m_A2_mfit2, m_A2_delta);
  double A3 = Fit(q2, m_A3_0, m_A3_mfit2, m_A3_delta);
  m_A3 = A3*q2/(2.0*m_m1*(m_m0+m_m1));
  m_V  = Fit(q2, m_V_0, m_V_mfit2, m_V_delta);
  m_calced = true;
}

double hepph0007169::Fit(double q2, double F0, double mfit2, double delta)
{
  double a=q2/mfit2;
  return F0/(1.0-a-delta*sqr(a));
}

} // namespace VA_P_V
} // namespace HADRONS
