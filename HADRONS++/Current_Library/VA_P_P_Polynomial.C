#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class Polynomial : public FormFactor_Base {
  // http://www.slac.stanford.edu/spires/find/hep/www?j=PHLTA,B36,521
  double m_fplus_0, m_fplus_lambda, m_fplus_m2;
  double m_f0_0, m_f0_lambda, m_f0_m2;
  double Fit(double q2, double m2, double f_0, double lambda);
public:
  Polynomial(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

Polynomial::Polynomial(GeneralModel model, double* masses,
                       const Flavour_Vector& flavs, std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_fplus_0=0.0; m_fplus_lambda=0.0; m_fplus_m2=0.0;
  m_f0_0=0.0; m_f0_lambda=0.0; m_f0_m2=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_K || kf0==kf_K_plus || kf0==kf_K_S || kf0==kf_K_L) {
    if (kf1==kf_pi || kf1==kf_pi_plus) {
      m_fplus_0  = 0.69;
      m_fplus_lambda  = 2.96e-2;
      m_fplus_m2 = 0.13957*0.13957;

      m_f0_0  = 0.69;
      m_f0_lambda = 1.96e-2;
      m_f0_m2 = 0.13957*0.13957;
    }
  }

  m_fplus_0 = model("fplus_0",m_fplus_0);
  m_fplus_lambda = model("fplus_lambda",m_fplus_lambda);
  m_fplus_m2 = model("fplus_m2",m_fplus_m2);
  
  m_f0_0 = model("f0_0",m_f0_0);
  m_f0_lambda = model("f0_lambda",m_f0_lambda);
  m_f0_m2 = model("f0_m2",m_f0_m2);
}

void Polynomial::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  m_fplus = Fit(q2, m_fplus_m2, m_fplus_0, m_fplus_lambda);
  m_f0 = Fit(q2, m_f0_m2, m_f0_0, m_f0_lambda);
  m_calced = true;
}

double Polynomial::Fit(double q2, double m2, double f_0, double lambda)
{
  return f_0*(1.0+lambda*q2/m2);
}

} // namespace VA_P_V
} // namespace HADRONS
