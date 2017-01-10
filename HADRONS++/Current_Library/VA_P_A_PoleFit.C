#include "HADRONS++/Current_Library/VA_P_A.H"

namespace HADRONS {
namespace VA_P_A_FFs {

class PoleFit : public FormFactor_Base {
double     F_0_A,m_fit_A,delta_A;
double     F_0_V_0,m_fit_V_0,delta_V_0;
double     F_0_V_1,m_fit_V_1,delta_V_1;
double     F_0_V_2,m_fit_V_2,delta_V_2;
public:
  PoleFit(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

PoleFit::PoleFit(GeneralModel model, double* masses, const Flavour_Vector& flavs,
           const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode(); 
  if (kf0==kf_B_c)
    if (kf1==10533){
        F_0_A=0.19,m_fit_A=1.71,delta_A=0.45;
        F_0_V_0=0.10,m_fit_V_0=-0.75,delta_V_0=0.95;
        F_0_V_1=5.28,m_fit_V_1=-2.28,delta_V_1=2.08;
        F_0_V_2=0.07,m_fit_V_2=1.73,delta_V_2=0.32;
    }
    else if (kf1==10513){
        F_0_A=0.21,m_fit_A=1.64,delta_A=0.49;
        F_0_V_0=0.13,m_fit_V_0=-2.48,delta_V_0=51.50;
        F_0_V_1=4.97,m_fit_V_1=-3.14,delta_V_1=6.49;
        F_0_V_2=0.09,m_fit_V_2=1.64,delta_V_2=0.38;   
        
    }
    


}

void PoleFit::CalcFFs( Vec4D p0, Vec4D p1 )
{
    
  double q2 = (p0-p1).Abs2();
  double m_fit_A2=m_fit_A*m_fit_A;
  double m_fit_V_02=m_fit_V_0*m_fit_V_0;
  double m_fit_V_12=m_fit_V_1*m_fit_V_1;
  double m_fit_V_22=m_fit_V_2*m_fit_V_2;
  if(m_fit_A<0){m_fit_A2=-m_fit_A2;}
  if(m_fit_V_0<0){m_fit_V_02=-m_fit_V_02;}
  if(m_fit_V_1<0){m_fit_V_12=-m_fit_V_12;}
  if(m_fit_V_2<0){m_fit_V_22=-m_fit_V_22;}

   
  m_A = F_0_A/(1.0-q2/(m_fit_A2)+delta_A*(q2/(m_fit_A2))*(q2/(m_fit_A2)));
  m_V0 = F_0_V_0/(1.0-q2/(m_fit_V_02)+delta_V_0*(q2/(m_fit_V_02))*(q2/(m_fit_V_02)));
  m_V1 = F_0_V_1/(1.0-q2/(m_fit_V_12)+delta_V_1*(q2/(m_fit_V_12))*(q2/(m_fit_V_12)));
  m_V2 = F_0_V_2/(1.0-q2/(m_fit_V_22)+delta_V_2*(q2/(m_fit_V_22))*(q2/(m_fit_V_22)));
    
    
  m_calced = true;
}

} // namespace VA_P_A
} // namespace HADRONS
