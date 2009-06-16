#include "HADRONS++/Current_Library/VA_P_V.H"

namespace HADRONS {
namespace VA_P_V_FFs {

class PoleFit : public FormFactor_Base {
  double m_V_F0, m_V_m0, m_V_a0, m_V_b0, m_V_c0, m_V_d0;
  double m_V_F1, m_V_m1, m_V_a1, m_V_b1, m_V_c1, m_V_d1;
  double m_A0_F0, m_A0_m0, m_A0_a0, m_A0_b0, m_A0_c0, m_A0_d0;
  double m_A0_F1, m_A0_m1, m_A0_a1, m_A0_b1, m_A0_c1, m_A0_d1;
  double m_A1_F0, m_A1_m0, m_A1_a0, m_A1_b0, m_A1_c0, m_A1_d0;
  double m_A1_F1, m_A1_m1, m_A1_a1, m_A1_b1, m_A1_c1, m_A1_d1;
  double m_A2_F0, m_A2_m0, m_A2_a0, m_A2_b0, m_A2_c0, m_A2_d0;
  double m_A2_F1, m_A2_m1, m_A2_a1, m_A2_b1, m_A2_c1, m_A2_d1;

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
  m_V_F0=0.0; m_V_m0=0.0; m_V_a0=0.0; m_V_b0=0.0; m_V_c0=0.0; m_V_d0=0.0;
  m_V_F1=0.0; m_V_m1=0.0; m_V_a1=0.0; m_V_b1=0.0; m_V_c1=0.0; m_V_d1=0.0;
  m_A0_F0=0.0; m_A0_m0=0.0; m_A0_a0=0.0; m_A0_b0=0.0; m_A0_c0=0.0; m_A0_d0=0.0;
  m_A0_F1=0.0; m_A0_m1=0.0; m_A0_a1=0.0; m_A0_b1=0.0; m_A0_c1=0.0; m_A0_d1=0.0;
  m_A1_F0=0.0; m_A1_m0=0.0; m_A1_a0=0.0; m_A1_b0=0.0; m_A1_c0=0.0; m_A1_d0=0.0;
  m_A1_F1=0.0; m_A1_m1=0.0; m_A1_a1=0.0; m_A1_b1=0.0; m_A1_c1=0.0; m_A1_d1=0.0;
  m_A2_F0=0.0; m_A2_m0=0.0; m_A2_a0=0.0; m_A2_b0=0.0; m_A2_c0=0.0; m_A2_d0=0.0;
  m_A2_F1=0.0; m_A2_m1=0.0; m_A2_a1=0.0; m_A2_b1=0.0; m_A2_c1=0.0; m_A2_d1=0.0;

  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B || kf0==kf_B_plus) {
    if (kf1==kf_rho_770 || kf1==kf_rho_770_plus) {
      m_V_F0 = 1.045;
      m_V_m0 = 5.32;
      m_V_F1 = -0.721;
      m_V_m1 = 6.192;
      m_V_a0 = -1.0;
      m_V_a1 = -1.0;

      m_A0_F0 = 1.527;
      m_A0_m0 = 5.28;
      m_A0_F1 = -1.220;
      m_A0_m1 = 5.776;
      m_A0_a0 = -1.0;
      m_A0_a1 = -1.0;

      m_A1_F0 = 0.240;
      m_A1_m0 = 6.125;
      m_A1_a0 = -1.0;

      m_A2_F0 = 0.009;
      m_A2_m0 = 6.389;
      m_A2_F1 = 0.212;
      m_A2_m1 = 6.389;
      m_A2_a0 = -1.0;
      m_A2_a1 = -2.0;
      m_A2_b1 = 1.0;
    }
    else if(kf1==kf_omega_782) {
      m_V_F0 = 1.006;
      m_V_m0 = 5.32;
      m_V_F1 = -0.713;
      m_V_m1 = 6.12;
      m_V_a0 = -1.0;
      m_V_a1 = -1.0;

      m_A0_F0 = 1.321;
      m_A0_m0 = 5.28;
      m_A0_F1 = -1.040;
      m_A0_m1 = 5.871;
      m_A0_a0 = -1.0;
      m_A0_a1 = -1.0;

      m_A1_F0 = -0.217;
      m_A1_m0 = 6.084;
      m_A1_a0 = -1.0;

      m_A2_F0 = 0.006;
      m_A2_m0 = 6.422;
      m_A2_F1 = 0.192;
      m_A2_m1 = 6.422;
      m_A2_a0 = -1.0;
      m_A2_a1 = -2.0;
      m_A2_b1 = 1.0;
    }
    else if(kf1==kf_a_1_1260 || kf1==kf_a_1_1260_plus) {
      m_V_F0=-0.67;
      m_V_m0=m_m0;
      m_V_a0=-0.72;
      m_V_b0=-0.2;
      
      m_A0_F0=-0.23;
      m_A0_m0=m_m0;
      m_A0_a0=-0.86;
      m_A0_b0=-0.38;

      m_A1_F0=-0.42;
      m_A1_m0=m_m0;
      m_A1_a0=0.44;
      m_A1_b0=0.45;

      m_A2_F0=-0.53;
      m_A2_m0=m_m0;
      m_A2_a0=-0.45;
      m_A2_b0=0.13;
    }
  }
  else if (kf0==kf_D_plus || kf0==kf_D) {
    if (kf1==kf_K_star_892 || kf1==kf_K_star_892_plus) {
      m_V_F0 = 1.62;
      m_V_a0 = -1.0;
      m_V_m0 = 2.1;

      m_A1_F0 = 1.0;
      m_A1_a0 = -1.0;
      m_A1_m0 = 2.5;

      m_A2_F0 = 0.83;
      m_A2_a0 = -1.0;
      m_A2_m0 = 2.5;
    }
  }
  else if (kf0==kf_D_s_plus) {
    if (kf1==kf_phi_1020) {
      m_V_F0 = 0.9;
      m_V_m0 = 1.9685;
      m_V_a0 = -2.82;
      m_V_b0 = 1.51;

      m_A0_F0 = 0.56;
      m_A0_m0 = 1.9685;
      m_A0_a0 = -0.13;
      m_A0_b0 = 0.46;

      m_A1_F0 = 0.65;
      m_A1_m0 = 1.9685;
      m_A1_a0 = -1.36;
      m_A1_b0 = -0.31;

      m_A2_F0 = 0.85;
      m_A2_m0 = 1.9685;
      m_A2_a0 = -4.5;
      m_A2_b0 = 5.55;
    }
  }

  m_V_F0 = model("V_F0",m_V_F0); m_V_F1 = model("V_F1",m_V_F1);
  m_V_m0 = model("V_m0",m_V_m0); m_V_m1 = model("V_m1",m_V_m1);
  m_V_a0 = model("V_a0",m_V_a0); m_V_a1 = model("V_a1",m_V_a1);
  m_V_b0 = model("V_b0",m_V_b0); m_V_b1 = model("V_b1",m_V_b1);
  m_V_c0 = model("V_c0",m_V_c0); m_V_c1 = model("V_c1",m_V_c1);
  m_V_d0 = model("V_d0",m_V_d0); m_V_d1 = model("V_d1",m_V_d1);

  m_A0_F0 = model("A0_F0",m_A0_F0); m_A0_F1 = model("A0_F1",m_A0_F1);
  m_A0_m0 = model("A0_m0",m_A0_m0); m_A0_m1 = model("A0_m1",m_A0_m1);
  m_A0_a0 = model("A0_a0",m_A0_a0); m_A0_a1 = model("A0_a1",m_A0_a1);
  m_A0_b0 = model("A0_b0",m_A0_b0); m_A0_b1 = model("A0_b1",m_A0_b1);
  m_A0_c0 = model("A0_c0",m_A0_c0); m_A0_c1 = model("A0_c1",m_A0_c1);
  m_A0_d0 = model("A0_d0",m_A0_d0); m_A0_d1 = model("A0_d1",m_A0_d1);

  m_A1_F0 = model("A1_F0",m_A1_F0); m_A1_F1 = model("A1_F1",m_A1_F1);
  m_A1_m0 = model("A1_m0",m_A1_m0); m_A1_m1 = model("A1_m1",m_A1_m1);
  m_A1_a0 = model("A1_a0",m_A1_a0); m_A1_a1 = model("A1_a1",m_A1_a1);
  m_A1_b0 = model("A1_b0",m_A1_b0); m_A1_b1 = model("A1_b1",m_A1_b1);
  m_A1_c0 = model("A1_c0",m_A1_c0); m_A1_c1 = model("A1_c1",m_A1_c1);
  m_A1_d0 = model("A1_d0",m_A1_d0); m_A1_d1 = model("A1_d1",m_A1_d1);

  m_A2_F0 = model("A2_F0",m_A2_F0); m_A2_F1 = model("A2_F1",m_A2_F1);
  m_A2_m0 = model("A2_m0",m_A2_m0); m_A2_m1 = model("A2_m1",m_A2_m1);
  m_A2_a0 = model("A2_a0",m_A2_a0); m_A2_a1 = model("A2_a1",m_A2_a1);
  m_A2_b0 = model("A2_b0",m_A2_b0); m_A2_b1 = model("A2_b1",m_A2_b1);
  m_A2_c0 = model("A2_c0",m_A2_c0); m_A2_c1 = model("A2_c1",m_A2_c1);
  m_A2_d0 = model("A2_d0",m_A2_d0); m_A2_d1 = model("A2_d1",m_A2_d1);
}

void PoleFit::CalcFFs( Vec4D p0, Vec4D p1 )
{
  double q2=(p0-p1).Abs2();
  m_V = Fit(q2,m_V_F0,m_V_m0,m_V_a0,m_V_b0,m_V_c0,m_V_d0,
               m_V_F1,m_V_m1,m_V_a1,m_V_b1,m_V_c1,m_V_d1);
  m_A0 = Fit(q2,m_A0_F0,m_A0_m0,m_A0_a0,m_A0_b0,m_A0_c0,m_A0_d0,
               m_A0_F1,m_A0_m1,m_A0_a1,m_A0_b1,m_A0_c1,m_A0_d1);
  m_A1 = Fit(q2,m_A1_F0,m_A1_m0,m_A1_a0,m_A1_b0,m_A1_c0,m_A1_d0,
               m_A1_F1,m_A1_m1,m_A1_a1,m_A1_b1,m_A1_c1,m_A1_d1);
  m_A2 = Fit(q2,m_A2_F0,m_A2_m0,m_A2_a0,m_A2_b0,m_A2_c0,m_A2_d0,
               m_A2_F1,m_A2_m1,m_A2_a1,m_A2_b1,m_A2_c1,m_A2_d1);
  m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 -
         (m_m0-m_m1)/(2.0*m_m1)*m_A2;
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

} // namespace VA_P_V
} // namespace HADRONS
