#include "HADRONS++/Current_Library/VA_P_T.H"

namespace HADRONS {
namespace VA_P_T_FFs {

class vector : public FormFactor_Base {
    double F_0_A0,m_fit_A0,delta_A0;
    double F_0_A1,m_fit_A1,delta_A1;
    double F_0_A2,m_fit_A2,delta_A2;
    double F_0_V,m_fit_V,delta_V;

public:
  vector(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       const std::vector<int>& i);
  void CalcFFs(ATOOLS::Vec4D p0, ATOOLS::Vec4D p1);
};

vector::vector(GeneralModel model, double* masses, const Flavour_Vector& flavs,
           const std::vector<int>& i) :
  FormFactor_Base(model, masses, flavs, i)
{
    kf_code kf0=m_flavs[p_i[0]].Kfcode();
    kf_code kf1=m_flavs[p_i[1]].Kfcode(); 

if (kf0==kf_B_c){
    if (kf1==535){          //das ist B_s_2*
        F_0_V=-18.60,m_fit_V=1.50,delta_V=0.48;
        F_0_A0=-2.94,m_fit_A0=1.47,delta_A0=0.54;
        F_0_A1=-2.89,m_fit_A1=1.75,delta_A1=0.48;
        F_0_A2=-1.32,m_fit_A2=-3.24,delta_A2=9.56; 
    } 
    else if(kf1==515){      //das ist B_2*
        F_0_V=17.60,m_fit_V=1.43,delta_V=0.52;
        F_0_A0=2.64,m_fit_A0=1.40,delta_A0=0.59;
        F_0_A1=2.59,m_fit_A1=1.68,delta_A1=0.52;
        F_0_A2=1.31,m_fit_A2=-3.13,delta_A2=9.72;  
        
     } 

        
    }
}

void vector::CalcFFs( Vec4D p0, Vec4D p1 )
{
  
      
    double q2 = (p0-p1).Abs2();
    double m_fit_V2=m_fit_V*m_fit_V;
    double m_fit_A02=m_fit_A0*m_fit_A0;
    double m_fit_A12=m_fit_A1*m_fit_A1;
    double m_fit_A22=m_fit_A2*m_fit_A2;
    if(m_fit_V<0){m_fit_V2=-m_fit_V2;}
    if(m_fit_A0<0){m_fit_A02=-m_fit_A02;}
    if(m_fit_A1<0){m_fit_A12=-m_fit_A12;}
    if(m_fit_A2<0){m_fit_A22=-m_fit_A22;}

  
  double V = F_0_V/(1-q2/(m_fit_V2)+delta_V*(q2/(m_fit_V2))*(q2/(m_fit_V2)));
  double A0 = F_0_A0/(1-q2/(m_fit_A02)+delta_A0*(q2/(m_fit_A02))*(q2/(m_fit_A02)));
  double A1 = F_0_A1/(1-q2/(m_fit_A12)+delta_A1*(q2/(m_fit_A12))*(q2/(m_fit_A12)));
  double A2 = F_0_A2/(1-q2/(m_fit_A22)+delta_A2*(q2/(m_fit_A22))*(q2/(m_fit_A22)));
  
  

   m_h=-V/(m_m0*(m_m0+m_m1));
   m_k=A1*(m_m0+m_m1)/m_m0;
   m_bplus=-A2/(m_m0*(m_m0+m_m1));
   m_bminus=-2.0*m_m1/m_m0*A0+(sqr(m_m0)-sqr(m_m1))/(q2*m_m0*(m_m0+m_m1))*A2;
    
  m_calced = true;


  
}

} // namespace VA_P_T
} // namespace HADRONS
