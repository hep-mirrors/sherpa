#include "HADRONS++/Current_Library/VA_P_S.H"

//hep-ph 1607.00622v1   
namespace HADRONS {
namespace VA_P_S_FFs {

class PoleFit : public FormFactor_Base {
  double     F_0_plus,m_fit_plus,delta_plus,F_0_0,m_fit_0,delta_0;
public:
  PoleFit(GeneralModel model, double* masses, const Flavour_Vector& flavs,
             std::vector<int>& indices);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

PoleFit::PoleFit(GeneralModel model, double* masses,
                       const Flavour_Vector& flavs, std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode();
    
    if (kf0==kf_B_c) {
        if (kf1==10531) {
        F_0_plus=0.71,m_fit_plus=1.69,delta_plus=0.48;
        F_0_0=0.72,m_fit_0=-1.98,delta_0=1.43;
        
        }
        else if (kf1==10511){
        F_0_plus=0.69,m_fit_plus=1.61,delta_plus=0.51;
        F_0_0=0.69,m_fit_0=-2.83,delta_0=4.84;

        }
        
    
    }
}

void PoleFit::CalcFFs( Vec4D p0, Vec4D p1)
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
