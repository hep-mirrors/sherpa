#include "HADRONS++/Current_Library/VA_P_P.H"

namespace HADRONS {
namespace VA_P_P_FFs {

class PoleFit : public FormFactor_Base {
  double m_fplus_F0, m_fplus_m0, m_fplus_a0, m_fplus_b0, m_fplus_c0, m_fplus_d0;
  double m_fplus_F1, m_fplus_m1, m_fplus_a1, m_fplus_b1, m_fplus_c1, m_fplus_d1;
  double m_f0_F0, m_f0_m0, m_f0_a0, m_f0_b0, m_f0_c0, m_f0_d0;
  double m_f0_F1, m_f0_m1, m_f0_a1, m_f0_b1, m_f0_c1, m_f0_d1;
  double m_fmin_F0, m_fmin_m0, m_fmin_a0, m_fmin_b0, m_fmin_c0, m_fmin_d0;
  double m_fmin_F1, m_fmin_m1, m_fmin_a1, m_fmin_b1, m_fmin_c1, m_fmin_d1;

  double Fit(double q2,
             double F0, double m0, double a0, double b0, double c0, double d0,
             double F1, double m1, double a1, double b1, double c1, double d1);
public:
  PoleFit(struct GeneralModel model, double* masses, Flavour* flavs, int* i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

PoleFit::PoleFit(GeneralModel model, double* masses, Flavour* flavs, int* i) :
  FormFactor_Base(model, masses, flavs, i)
{
  m_fplus_F0=0.0; m_fplus_m0=0.0; m_fplus_a0=0.0; m_fplus_b0=0.0;
  m_fplus_c0=0.0; m_fplus_d0=0.0;
  m_fplus_F1=0.0; m_fplus_m1=0.0; m_fplus_a1=0.0; m_fplus_b1=0.0;
  m_fplus_c1=0.0; m_fplus_d1=0.0;
  m_f0_F0=0.0; m_f0_m0=0.0; m_f0_a0=0.0; m_f0_b0=0.0; m_f0_c0=0.0; m_f0_d0;
  m_f0_F1=0.0; m_f0_m1=0.0; m_f0_a1=0.0; m_f0_b1=0.0; m_f0_c1=0.0; m_f0_d1;
  m_fmin_F0=0.0; m_fmin_m0=0.0; m_fmin_a0=0.0; m_fmin_b0=0.0;
  m_fmin_c0=0.0; m_fmin_d0=0.0;
  m_fmin_F1=0.0; m_fmin_m1=0.0; m_fmin_a1=0.0; m_fmin_b1=0.0;
  m_fmin_c1=0.0; m_fmin_d1=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B || kf0==kf_B_plus) {
    if (kf1==kf_pi || kf1==kf_pi_plus) {
      m_fplus_F0 = 0.744;
      m_fplus_m0 = 5.32;
      m_fplus_a0 = -1;
      m_fplus_F1 = -0.486;
      m_fplus_m1 = 6.382;
      m_fplus_a1 = -1;

      m_f0_F0 = 0.258;
      m_f0_m0 = 5.815;
      m_f0_a0 = -1;
    }
    else if (kf1==kf_eta) {
      m_fplus_F0 = 0.122;
      m_fplus_m0 = 5.32;
      m_fplus_a0 = -1;
      m_fplus_F1 = 0.155;
      m_fplus_m1 = 5.32;
      m_fplus_a1 = -2;
      m_fplus_b1 = 1;

      m_f0_F0 = 0.273;
      m_f0_m0 = 5.57;
      m_f0_a0 = -1;
    }
    else if (kf1==kf_f_0_980) {
      // hep-ph 0701108
      m_fplus_F0 = 1.7*0.25;
      m_fplus_m0 = m_m0;
      m_fplus_a0 = -0.48;
      m_fplus_b0 = -0.30;
      m_fplus_c0 = 0.47;
      m_fplus_a0 = -0.99;

      m_fmin_F0 = -1.7*0.24;
      m_fmin_m0 = m_m0;
      m_fmin_a0 = -0.41;
      m_fmin_b0 = -0.42;
      m_fmin_c0 = 0.95;
      m_fmin_a0 = -1.55;
    }
  }
  else if (kf0==kf_D_plus) {
    if (kf1==kf_f_0_980) {
      // hep-ph 0701108
      m_fplus_F0 = 1.7*0.32;
      m_fplus_m0 = m_m0;
      m_fplus_a0 = -0.89;
      m_fplus_b0 = -0.40;
      m_fplus_c0 = -0.18;
      m_fplus_a0 = -1.00;
    }
  }
  else if (kf0==kf_D_s_plus) {
    if (kf1==kf_f_0_980) {
      // hep-ph 0701108
      m_fplus_F0 = 1.7*0.27;
      m_fplus_m0 = m_m0;
      m_fplus_a0 = -0.87;
      m_fplus_b0 = -0.17;
      m_fplus_c0 = -0.37;
      m_fplus_a0 = 1.46;
    }
  }
  else if (kf0==kf_B_c) {
    // hep-ph/0007169
    if (kf1==kf_B_s) {
      m_fplus_F0=-0.61;
      m_fplus_m0=1.73;
      m_fplus_a0=-1.0;
      m_fplus_b0=0.09;

      m_fmin_F0=1.83;
      m_fmin_m0=2.21;
      m_fmin_a0=-1.0;
      m_fmin_b0=-0.07;
    }
    else if (kf1==kf_eta_c_1S) {
      m_fplus_F0=0.76;
      m_fplus_m0=6.37;
      m_fplus_a0=-1.0;
      m_fplus_b0=-0.087;

      m_fmin_F0=-0.38;
      m_fmin_m0=6.22;
      m_fmin_a0=-1.0;
      m_fmin_b0=-0.06;
    }
  }

  m_fplus_F0 = model("fplus_F0",m_fplus_F0);
  m_fplus_F1 = model("fplus_F1",m_fplus_F1);
  m_fplus_m0 = model("fplus_m0",m_fplus_m0);
  m_fplus_m1 = model("fplus_m1",m_fplus_m1);
  m_fplus_a0 = model("fplus_a0",m_fplus_a0);
  m_fplus_a1 = model("fplus_a1",m_fplus_a1);
  m_fplus_b0 = model("fplus_b0",m_fplus_b0);
  m_fplus_b1 = model("fplus_b1",m_fplus_b1);
  m_fplus_c0 = model("fplus_c0",m_fplus_c0);
  m_fplus_c1 = model("fplus_c1",m_fplus_c1);
  m_fplus_d0 = model("fplus_d0",m_fplus_d0);
  m_fplus_d1 = model("fplus_d1",m_fplus_d1);

  m_fmin_F0 = model("fplus_F0",m_fmin_F0);
  m_fmin_F1 = model("fplus_F1",m_fmin_F1);
  m_fmin_m0 = model("fplus_m0",m_fmin_m0);
  m_fmin_m1 = model("fplus_m1",m_fmin_m1);
  m_fmin_a0 = model("fplus_a0",m_fmin_a0);
  m_fmin_a1 = model("fplus_a1",m_fmin_a1);
  m_fmin_b0 = model("fplus_b0",m_fmin_b0);
  m_fmin_b1 = model("fplus_b1",m_fmin_b1);
  m_fmin_c0 = model("fplus_c0",m_fmin_c0);
  m_fmin_c1 = model("fplus_c1",m_fmin_c1);
  m_fmin_d0 = model("fplus_d0",m_fmin_d0);
  m_fmin_d1 = model("fplus_d1",m_fmin_d1);

  m_f0_F0 = model("f0_F0",m_f0_F0); m_f0_F1 = model("f0_F1",m_f0_F1);
  m_f0_m0 = model("f0_m0",m_f0_m0); m_f0_m1 = model("f0_m1",m_f0_m1);
  m_f0_a0 = model("f0_a0",m_f0_a0); m_f0_a1 = model("f0_a1",m_f0_a1);
  m_f0_b0 = model("f0_b0",m_f0_b0); m_f0_b1 = model("f0_b1",m_f0_b1);
  m_f0_c0 = model("f0_c0",m_f0_c0); m_f0_c1 = model("f0_c1",m_f0_c1);
  m_f0_d0 = model("f0_d0",m_f0_d0); m_f0_d1 = model("f0_d1",m_f0_d1);
}

void PoleFit::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  m_fplus = Fit(q2,m_fplus_F0,m_fplus_m0,m_fplus_a0,m_fplus_b0,m_fplus_c0,m_fplus_d0,
                m_fplus_F1,m_fplus_m1,m_fplus_a1,m_fplus_b1,m_fplus_c1,m_fplus_d1);
  if (m_f0_F0!=0.0) {
    m_f0 = Fit(q2,m_f0_F0,m_f0_m0,m_f0_a0,m_f0_b0,m_f0_c0,m_f0_d0,
               m_f0_F1,m_f0_m1,m_f0_a1,m_f0_b1,m_f0_c1,m_f0_d1);
  }
  else if (m_fmin_F0!=0.0) {
    double fmin = Fit(q2,m_fmin_F0,m_fmin_m0,m_fmin_a0,m_fmin_b0,m_fmin_c0,
                      m_fmin_d0,m_fmin_F1,m_fmin_m1,m_fmin_a1,m_fmin_b1,
                      m_fmin_c1,m_fmin_d1);
    m_f0 = q2/(sqr(m_m0)-sqr(m_m1))*(fmin+m_fplus);
  }
  else m_f0=0.0;
  m_calced = true;
}

double PoleFit::Fit(double q2,
             double F0, double m0, double a0, double b0, double c0, double d0,
             double F1, double m1, double a1, double b1, double c1, double d1)
{
  double fit=0.0;
  double m02=m0*m0;
  double m12=m1*m1;
  if(F0!=0.0 && m02!=0.0) {
    double q=q2/m02;
    fit+=F0/(1.0 + a0*q + b0*sqr(q) + c0*pow(q,3) + d0*pow(q,4));
  }
  if(F1!=0.0 && m12!=0.0) {
    double q=q2/m12;
    fit+=F1/(1.0 + a1*q + b1*sqr(q) + c1*pow(q,3) + d1*pow(q,4));
  }
  return fit;
}

} // namespace VA_P_P
} // namespace HADRONS
