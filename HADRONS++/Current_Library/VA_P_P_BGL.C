#include "HADRONS++/Current_Library/VA_P_P.H"
#include <cmath>
// BGL parametrization -> arXiv:hep-ph/9508211

namespace HADRONS {
namespace VA_P_P_FFs {

class BGL : public FormFactor_Base {
 int case_number; // case_number selects the used phi- and f_plus function
 //Variables used by case 0 arXiv:1208.1253 [hep-ex] and case 2, arXiv:1510.03657 [hep-ex]
 // a0 is for m_fplus, a0_0 is for m_f0
 double mstar ; double a0; double a1; double a2; double a3; double a0_0; double a1_0; double a2_0; double a3_0; double q02;
 double fplus_0; double r1; double r2; // Variables used by case 1
 double zCalc(double q2, double t_0, double t_plus, double t_minus);

public:
  BGL(GeneralModel model, double* masses, const Flavour_Vector& flavs,
        std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

BGL::BGL(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices) :
  FormFactor_Base(model, masses, flavs, indices)
{
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  // B-Decay:
  if (kf0==kf_B) {
    if (kf1==kf_pi_plus) {
    case_number = 0;
    const double mminus = m_m0-m_m1;
    mstar = 5.32483;
    a0 = 0.026;
    a1 = -1.15*a0;
    a2 = -4.52*a0;
    a3 = 0;
    q02 = 0.65*mminus*mminus;  
    }
    else if (kf1==kf_D_plus) {
      case_number = 2;
      a0 = 0.0126;
      a1 = -0.094;
      a2 = 0.34;
      a3 = -0.1;

      a0_0  = 0.0115;
      a1_0  = -0.057;
      a2_0  = 0.12;
      a3_0  = 0.4;
     }
  }
  // B+-Decay:
   else if (kf0==kf_B_plus) {
     if (kf1==kf_pi) {
       const double mminus = m_m0-m_m1;
       case_number = 0;
       mstar = 5.32483;
       a0 = 0.026;
       a1 = -0.63*a0;
       a2 = -5.80*a0;
       a3 = 0;
       q02 = 0.65*mminus*mminus;
    }
    else if (kf1==kf_eta) {
      case_number = 0;
      mstar=5.32483;
      a0 = 0.0031;
      a1 = -0.005301;
      q02 = 14.14;
    }
    else if (kf1==kf_D) {
      case_number = 2;
      a0 = 0.0126;
      a1 = -0.094;
      a2 = 0.34;
      a3 = -0.1;

      a0_0 = 0.0115;
      a1_0 = -0.057;
      a2_0 = 0.12;
      a3_0 = 0.4;
     }
  }
  // D0-Decay:
  else if (kf0==kf_D) {
    if (kf1==kf_pi_plus) {
      case_number = 1;
      fplus_0 = 0.6748956576;
      r1 = -2.0;
      r2 = -1.6;
    }
    else if (kf1==kf_K_plus) {
      case_number = 1;
      fplus_0 = 0.7653349496;  
      r1 = -2.4;
      r2 = 15.6;
    }
  }
  // D+-Decay
  else if (kf0==kf_D_plus) {
    if (kf1==kf_pi) {
      case_number = 1;
      fplus_0 = 0.6340467099; 
      r1 = -1.95;
      r2 = -0.11;
    }
    else if (kf1==kf_K) {
      case_number = 1;
      fplus_0 = 0.738;
      r1 = -1.9;
      r2 = 16.6;
    }
  }
  
  mstar = model("mstar",mstar);
  a0 = model("a0",a0);
  a1 = model("a1",a1);
  a2 = model("a2",a2);
  a3 = model("a3",a3); 
  a0_0 = model("a0_0",a0_0);
  a1_0 = model("a1_0",a1_0);
  a2_0 = model("a2_0",a2_0);
  a3_0 = model("a3_0",a3_0); 
  fplus_0 = model("fplus_0",fplus_0);
  r1 = model("r1",r1);
  r2 = model("r2",r2);
  case_number = model("case_number",case_number);
}

void BGL::CalcFFs( Vec4D p0, Vec4D p1 )
{
  if (case_number==0) {
    const double mplus = m_m0+m_m1;
    const double mminus = m_m0-m_m1;
    const double chi = pow(6.889,-4);
    double q2 = (p0-p1).Abs2();
    double z = (sqrt(mplus*mplus-q2)-sqrt(mplus*mplus-q02))/(sqrt(mplus*mplus-q2)+sqrt(mplus*mplus-q02));
    // Blasche factor P
    double P = (sqrt(mplus*mplus-q2)-sqrt(mplus*mplus-mstar*mstar ))/(sqrt(mplus*mplus-q2)+sqrt(mplus*mplus-q02));
    double phi = sqrt(1/(32*M_PI*chi))*(sqrt(mplus*mplus-q2)+sqrt(mplus*mplus-q02))*pow(sqrt(mplus*mplus-q2)+sqrt(mplus*mplus-(mminus*mminus)),1.5)*pow(sqrt(mplus*mplus-q2)+mplus,-5)*(mplus*mplus-q2)/pow(mplus*mplus-q02,0.25);  
    m_fplus  = 1/(P*phi) * (a0 + a1*z + a2*z*z + a3*z*z*z);
    m_f0     = 0.0;
    m_calced = true;
  }

  else if (case_number==1) {
    const double m_Dstar_plus = 2.1121;
    const double zero = 0.0;
    const double m_c = 1.29;

    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();

    double q2 = (p0 - p1).Abs2();
    double m2_c = pow(m_c,2);
    double t_plus = pow(m_m0 + m_m1,2);
    double t_minus = pow(m_m0 - m_m1,2);
    double t_0 = t_plus * (1 - sqrt(1 - t_minus/t_plus));
    double m2_Dstar_plus = pow(m_Dstar_plus,2);
    double c = sqrt(M_PI * m2_c/3);
    // c fuer Parameter von Paper: http://arxiv.org/pdf/1412.5502v1.pdf
    //   double c = sqrt(4 * M_PI * m2_c/9);
    double P_q2;
    double P_0;

    if (kf1==kf_pi || kf1==kf_pi_plus) { P_q2 = 1; }
    else if (kf1==kf_K || kf1==kf_K_plus) { P_q2 = zCalc(q2, m2_Dstar_plus, t_plus, t_minus); }

    if (kf1==kf_pi || kf1==kf_pi_plus) { P_0 = 1; }
    else if (kf1==kf_K || kf1==kf_K_plus) { P_0 = zCalc(zero, m2_Dstar_plus, t_plus, t_minus); }
   
    // Berechnung von Phi:
    double phi = c * pow(zCalc(q2, zero, t_plus, t_minus)/(-1 * q2),2.5) *
      pow(zCalc(q2, t_0, t_plus, t_minus)/(t_0 - q2),-0.5) *
      pow(zCalc(q2, t_minus, t_plus, t_minus)/(t_minus - q2),-0.75) *
      (t_plus - q2)/(pow(t_plus - t_0,0.25));
    // Berechnung von Phi(0):
    double phi_0 = c * pow(zCalc(zero, t_0, t_plus, t_minus)/(t_0 - zero),-0.5) *
      pow(zCalc(zero, t_minus, t_plus, t_minus)/(t_minus - zero),-0.75) *
      (t_plus - zero)/(pow(t_plus - t_0,0.25));
    // Berechnung von z:
    double z = zCalc(q2, t_0, t_plus, t_minus);
    // Berechnung von z(0)
    double z_0 = zCalc(zero, t_0, t_plus, t_minus);

    m_fplus = (fplus_0 * P_0 * phi_0 * (1 + r1*z + r2*z*z)) / (P_q2 * phi * (1 + r1*z_0 + r2*z_0*z_0));
    m_f0 = 0.0;
    m_calced = true;
  }

  else if (case_number==2) {
    double r = m_m1/m_m0;
    double q2 = (p0-p1).Abs2();
    double w = (m_m0*m_m0+m_m1*m_m1-q2)/(2*m_m0*m_m1);
    double z = (sqrt(w+1)-sqrt(2))/(sqrt(w+1)+sqrt(2));
    double phi_plus = 1.1213*pow(1+z,2)*pow(1-z,0.5)*pow((1+r)*(1-z)+2*sqrt(r)*(1+z),-5);  
    double phi_0 = 0.5299*(1+z)*pow(1-z,1.5)*pow((1+r)*(1-z)+2*sqrt(r)*(1+z),-4); 
    m_fplus  = 1/phi_plus*(a0+a1*z+a2*z*z+a3*z*z*z);
    m_f0     = 1/phi_0*(a0_0+a1_0*z+a2_0*z*z+a3_0*z*z*z);
    m_calced = true;
  }
}

double BGL::zCalc(double q2, double t_0, double t_plus, double t_minus)
{
  double z = (sqrt(t_plus - q2) - sqrt(t_plus - t_0))/(sqrt(t_plus - q2) + sqrt(t_plus - t_0));
  return z;
}

} // namespace VA_P_V
} // namespace HADRONS
