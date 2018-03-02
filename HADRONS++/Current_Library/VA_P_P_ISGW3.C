#include "HADRONS++/Current_Library/VA_P_P.H"
#include <cmath>
namespace HADRONS {
namespace VA_P_P_FFs {

class ISGW3 : public FormFactor_Base {
 double f_plus0; double r; double fq2_max;

public:
  ISGW3(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

ISGW3::ISGW3(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices) :
  FormFactor_Base(model, masses, flavs, indices)
{

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  //D0-Decay:
  if (kf0==kf_D) {
    if (kf1==kf_pi_plus) {
      f_plus0 = 0.635;
      r = 2.0688;
      fq2_max = 2.7306;
    }
    else if (kf1==kf_K_plus) {
      f_plus0 = 0.7499255211;
      r = 1.56;
      fq2_max = 1.560461461;
    }
  }
  //D+-Decay
  else if (kf0==kf_D_plus) {
    if (kf1==kf_pi) {
      f_plus0 = 0.635;
      r = 2.01;
      fq2_max = 2.0731;
    }
    else if (kf1==kf_K) {
      f_plus0 = 0.7427344545;
      r = 1.48;
      fq2_max = 1.476223252;
    }
  }	
  
  f_plus0 = model("f_plus0",f_plus0);
  r = model("r",r); 
  fq2_max = model("fq2_max",fq2_max);
}

void ISGW3::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2_max = pow(m_m0 - m_m1,2);
  double q2 = (p0-p1).Abs2();

  m_fplus = fq2_max * pow(1+(pow(r,2)/12)*(q2_max - q2),-2);
  m_f0 = 0.0;
  m_calced = true;
}

} // namespace VA_P_V
} // namespace HADRONS
// taken from http://arxiv.org/pdf/0810.3878.pdf
