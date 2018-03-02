#include "HADRONS++/Current_Library/VA_P_V.H"

// hep-ph  1607.00622v1
namespace HADRONS {
namespace VA_P_V_FFs {

class PoleFit2 : public FormFactor_Base {
double     F_0_V,m_fit_V,delta_V;
double     F_0_A_0,m_fit_A_0,delta_A_0;
double     F_0_A_1,m_fit_A_1,delta_A_1;
double     F_0_A_2,m_fit_A_2,delta_A_2;


public:
  PoleFit2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
          const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

PoleFit2::PoleFit2(GeneralModel model, double* masses,
                 const Flavour_Vector& flavs, const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)

{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();
    
    if (kf0==kf_B_c) {
        if (kf1==kf_B_s_star) {
        
            F_0_V=3.70,m_fit_V=1.57,delta_V=0.48;
            F_0_A_0=0.55,m_fit_A_0=1.49,delta_A_0=0.61;
            F_0_A_1=0.52,m_fit_A_1=1.90,delta_A_1=0.56;
            F_0_A_2=0.07,m_fit_A_2=-1.04,delta_A_2=0.37;
            }
        else if (kf1==10533){
            F_0_V=-18.60,m_fit_V=1.50,delta_V=0.48;
            F_0_A_0=-2.94,m_fit_A_0=1.47,delta_A_0=0.54;
            F_0_A_1=-2.89,m_fit_A_1=1.75,delta_A_1=0.48;
            F_0_A_2=-1.32,m_fit_A_2=-3.24,delta_A_2=9.56;
            }
        else if (kf1==kf_B_star){
            F_0_V=3.44,m_fit_V=1.50,delta_V=0.51;
            F_0_A_0=0.47,m_fit_A_0=1.42,delta_A_0=0.68;
            F_0_A_1=0.44,m_fit_A_1=1.84,delta_A_1=0.63;
            F_0_A_2=0.07,m_fit_A_2=-1.03,delta_A_2=0.37;
            }
    }
}

void PoleFit2::CalcFFs( Vec4D p0, Vec4D p1)
{
   
    double q2 = (p0-p1).Abs2();
    double m_fit_V2=m_fit_V*m_fit_V;
    double m_fit_A_02=m_fit_A_0*m_fit_A_0;
    double m_fit_A_12=m_fit_A_1*m_fit_A_1;
    double m_fit_A_22=m_fit_A_2*m_fit_A_2;
    if(m_fit_V<0){m_fit_V2=-m_fit_V2;}
    if(m_fit_A_0<0){m_fit_A_02=-m_fit_A_02;}
    if(m_fit_A_1<0){m_fit_A_12=-m_fit_A_12;}
    if(m_fit_A_2<0){m_fit_A_22=-m_fit_A_22;}

  
  m_V = F_0_V/(1-q2/(m_fit_V2)+delta_V*(q2/(m_fit_V2))*(q2/(m_fit_V2)));
  m_A0 = F_0_A_0/(1-q2/(m_fit_A_02)+delta_A_0*(q2/(m_fit_A_02))*(q2/(m_fit_A_02)));
  m_A1 = F_0_A_1/(1-q2/(m_fit_A_12)+delta_A_1*(q2/(m_fit_A_12))*(q2/(m_fit_A_12)));
  m_A2 = F_0_A_2/(1-q2/(m_fit_A_22)+delta_A_2*(q2/(m_fit_A_22))*(q2/(m_fit_A_22)));
  m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
  
  m_calced = true;

}
}
}
