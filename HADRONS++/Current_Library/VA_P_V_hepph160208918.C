#include "HADRONS++/Current_Library/VA_P_V.H"


namespace HADRONS {
namespace VA_P_V_FFs {

class hepph160208918 : public FormFactor_Base {
double     F_0_V,a_V,b_V;
double     F_0_A_0,a_A_0,b_A_0;
double     F_0_A_1,a_A_1,b_A_1;
double     F_0_A_2,a_A_2,b_A_2;
double Fit(double q2, double F_0, double a, double b);

public:
 hepph160208918(GeneralModel model, double* masses, const Flavour_Vector& flavs,
          const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

hepph160208918::hepph160208918(GeneralModel model, double* masses,
                 const Flavour_Vector& flavs, const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)

{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();
    
    if (kf0==kf_B_c) {
        if (kf1==kf_J_psi_1S) {        
            F_0_V=1.59,a_V=5.04,b_V=5.88;
            F_0_A_0=0.78,a_A_0=5.41,b_A_0=10.86;
            F_0_A_1=0.96,a_A_1=5.24,b_A_1=-15.18;
            F_0_A_2=1.36,a_A_2=7.60,b_A_2=-5.94;
            }
        else if (kf1==100443){
            F_0_V=1.71,a_V=3.43,b_V=9.79;
            F_0_A_0=0.80,a_A_0=5.14,b_A_0=-32.16;
            F_0_A_1=0.87,a_A_1=5.45,b_A_1=-100.23;
            F_0_A_2=1.22,a_A_2=12.74,b_A_2=-214.39;      
        }
    

    }
}

void hepph160208918::CalcFFs( Vec4D p0, Vec4D p1)
{
   
    double q2 = (p0-p1).Abs2();
    m_V=Fit(q2, F_0_V, a_V, b_V);
    m_A0=Fit(q2, F_0_A_0, a_A_0, b_A_0);
    m_A1=Fit(q2, F_0_A_1, a_A_1, b_A_1);
    m_A2=Fit(q2, F_0_A_2, a_A_2, b_A_2);
    m_A3 = (m_m0+m_m1)/(2.0*m_m1)*m_A1 - (m_m0-m_m1)/(2.0*m_m1)*m_A2;
 
  m_calced = true;
    
}
double hepph160208918::Fit(double q2, double F_0, double a, double b)
{

return F_0*exp(a*q2/sqr(m_m0)+b*sqr(q2/sqr(m_m0)));

}





}
}
