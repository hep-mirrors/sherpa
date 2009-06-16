#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class SumRules3 : public FormFactor_Base {
  //hep-ph/0608264
  double* p_masses;
  double m_fV_0, m_fV_a, m_fV_b, m_fV_c, m_fV_d;
  double m_f0_0, m_f0_a, m_f0_b, m_f0_c, m_f0_d;
  double m_fplus_0, m_fplus_a, m_fplus_b, m_fplus_c, m_fplus_d;
  double m_fminus_0, m_fminus_a, m_fminus_b, m_fminus_c, m_fminus_d;

  double Fit(double q2, double f0, double a, double b, double c, double d);
public:
  SumRules3(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

SumRules3::SumRules3(GeneralModel model, double* masses, Flavour* flavs,
                     int* i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_fV_0=0.0; m_fV_a=0.0; m_fV_b=0.0; m_fV_c=0.0; m_fV_d=0.0;
  m_f0_0=0.0; m_f0_a=0.0; m_f0_b=0.0; m_f0_c=0.0; m_f0_d=0.0;
  m_fplus_0=0.0; m_fplus_a=0.0; m_fplus_b=0.0; m_fplus_c=0.0; m_fplus_d=0.0;
  m_fminus_0=0.0; m_fminus_a=0.0; m_fminus_b=0.0; m_fminus_c=0.0;m_fminus_d=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B_s) {
    if (kf1==kf_D_s1_H) {
      m_fV_0 = 1.18;
      m_fV_a = -1.87;
      m_fV_b = -1.88;
      m_fV_c = -2.41;
      m_fV_d = 3.34;

      m_f0_0 = 0.076;
      m_f0_a = 1.85;
      m_f0_b = 0.89;
      m_f0_c = 19.0;
      m_f0_d = -79.3;

      m_fplus_0 = 0.13;
      m_fplus_a = -7.14;
      m_fplus_b = 11.6;
      m_fplus_c = 21.3;
      m_fplus_d = -59.8;

      m_fminus_0 = -0.26;
      m_fminus_a = -4.11;
      m_fminus_b = -3.27;
      m_fminus_c = 15.2;
      m_fminus_d = 18.6;
    }
  }

  m_fV_0 = model("fV_0",m_fV_0);
  m_fV_a = model("fV_a",m_fV_a);
  m_fV_b = model("fV_b",m_fV_b);
  m_fV_c = model("fV_c",m_fV_c);
  m_fV_d = model("fV_d",m_fV_d);
  
  m_f0_0 = model("f0_0",m_f0_0);
  m_f0_a = model("f0_a",m_f0_a);
  m_f0_b = model("f0_b",m_f0_b);
  m_f0_c = model("f0_c",m_f0_c);
  m_f0_d = model("f0_d",m_f0_d);
  
  m_fplus_0 = model("fplus_0",m_fplus_0);
  m_fplus_a = model("fplus_a",m_fplus_a);
  m_fplus_b = model("fplus_b",m_fplus_b);
  m_fplus_c = model("fplus_c",m_fplus_c);
  m_fplus_d = model("fplus_d",m_fplus_d);
  
  m_fminus_0 = model("fminus_0",m_fminus_0);
  m_fminus_a = model("fminus_a",m_fminus_a);
  m_fminus_b = model("fminus_b",m_fminus_b);
  m_fminus_c = model("fminus_c",m_fminus_c);
  m_fminus_d = model("fminus_d",m_fminus_d);      
}

void SumRules3::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  double fV = Fit(q2, m_fV_0, m_fV_a, m_fV_b, m_fV_c, m_fV_d);
  double f0 = Fit(q2, m_f0_0, m_f0_a, m_f0_b, m_f0_c, m_f0_d);
  double fplus = Fit(q2, m_fplus_0, m_fplus_a, m_fplus_b, m_fplus_c, m_fplus_d);
  double fminus = Fit(q2, m_fminus_0, m_fminus_a, m_fminus_b, m_fminus_c, m_fminus_d);

  m_A1 = -f0;
  m_A2 = fplus;
  m_V  = -0.5*fV;
  m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 -
         (m_m0-m_m1)/(2.0*m_m1)*m_A2;
  m_A0 = m_A3+q2/(2.0*m_m1)*fminus;
  m_calced = true;
}

double SumRules3::Fit(double q2, double f0, double a, double b, double c, double d)
{
  double qhat = q2/sqr(m_m0);
  return f0/(1.0+a*qhat+b*sqr(qhat)+c*pow(qhat,3)+d*pow(qhat,4));
}

} // namespace VA_P_V
} // namespace HADRONS
