#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class BallZwicky : public FormFactor_Base {
  double m_fplus_0, m_fplus_r, m_fplus_m2, m_fplus_alpha;
public:
  BallZwicky(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

BallZwicky::BallZwicky(GeneralModel model, double* masses, Flavour* flavs,
                       int* i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_fplus_0=0.0; m_fplus_r=0.0; m_fplus_m2=0.0; m_fplus_alpha=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B || kf0==kf_B_plus) {
    if (kf1==kf_eta_prime_958) {
      m_fplus_0 = 0.189;
      m_fplus_alpha = 0.851;
      m_fplus_r = 0.411;
      m_fplus_m2 = 5.33*5.33;
    }
    else if (kf1==kf_eta) {
      m_fplus_0 = 0.231;
      m_fplus_alpha = 0.851;
      m_fplus_r = 0.411;
      m_fplus_m2 = 5.33*5.33;
    }
  }


  m_fplus_0 = model("fplus_0",m_fplus_0);
  m_fplus_r = model("fplus_r",m_fplus_r);
  m_fplus_m2 = model("fplus_m2",m_fplus_m2);
  m_fplus_alpha = model("fplus_alpha",m_fplus_alpha);
}

void BallZwicky::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  m_fplus = 1.0/(1.0-q2/m_fplus_m2)+
            m_fplus_r*q2/m_fplus_m2/(1.0-q2/m_fplus_m2)/(1.0-m_fplus_alpha*q2/m_fplus_m2);
  m_fplus *= m_fplus_0;
  m_f0 = 0.0;
  m_calced = true;
}

}
} // namespace HADRONS
