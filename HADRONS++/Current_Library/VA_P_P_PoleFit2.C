#include "HADRONS++/Current_Library/VA_P_P.H"

//hep-ph 1607.00622v1   
namespace HADRONS {
namespace VA_P_P_FFs {

class PoleFit2 : public FormFactor_Base {
  double     F_0_plus,m_fit_plus,delta_plus,F_0_0,m_fit_0,delta_0;
public:
  PoleFit2(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

PoleFit2::PoleFit2(GeneralModel model, double* masses,
                       const Flavour_Vector& flavs, std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();
    
    if (kf0==kf_B_c) {
        if (kf1==kf_B_s) {
        F_0_plus=0.73,m_fit_plus=1.57,delta_plus=0.49;
        F_0_0=0.73,m_fit_0=2.07,delta_0=0.82;
        }
        else if (kf1==kf_B){
        
        F_0_plus=0.64,m_fit_plus=1.50,delta_plus=0.52;
        F_0_0=0.64,m_fit_0=1.94,delta_0=0.83;
        }
    
    }
}

void PoleFit2::CalcFFs( Vec4D p0, Vec4D p1)
{

    double q2 = (p0-p1).Abs2();
    double m_fit_plus2=m_fit_plus*m_fit_plus;
    double m_fit_02=m_fit_0*m_fit_0;
    
    if(m_fit_plus<0){m_fit_plus2=-m_fit_plus2;}
    if(m_fit_0<0){m_fit_02=-m_fit_02;}
    

    m_fplus = F_0_plus/(1.0-q2/m_fit_plus2+delta_plus*(q2/m_fit_plus2)*(q2/m_fit_plus2));
    m_f0 = F_0_0/(1.0-q2/m_fit_02+delta_0*(q2/m_fit_02)*(q2/m_fit_02));
    
    m_calced = true;

}
}
}
