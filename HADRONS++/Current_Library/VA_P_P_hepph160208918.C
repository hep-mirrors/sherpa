#include "HADRONS++/Current_Library/VA_P_P.H"

  
namespace HADRONS {
namespace VA_P_P_FFs {

class hepph160208918 : public FormFactor_Base {
  double     F_0_plus,a_plus,b_plus,F_0_0,a_0,b_0;
  double     m_0,m_1;
public:
  hepph160208918(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

hepph160208918::hepph160208918(GeneralModel model, double* masses,
                       const Flavour_Vector& flavs, std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();
    
    if (kf0==kf_B_c) {
        
        if (kf1==kf_eta_c_1S) {
        F_0_plus=1.06,a_plus=4.18,b_plus=10.46;
        F_0_0=1.06,a_0=3.36,b_0=10.21;
        }
        else if (kf1==100441){
        F_0_plus=1.04,a_plus=7.62,b_plus=-96.66;
        F_0_0=1.04,a_0=6.40,b_0=-98.88;
        }
    
    }
}

void hepph160208918::CalcFFs( Vec4D p0, Vec4D p1)
{

    double q2 = (p0-p1).Abs2();
    
    
    m_fplus = F_0_plus*exp(a_plus*q2/sqr(m_m0)+b_plus*sqr(q2/sqr(m_m0)));
    m_f0 = F_0_0*exp(a_0*q2/sqr(m_m0)+b_0*sqr(q2/sqr(m_m0)));
    
    m_calced = true;

}
}
}
