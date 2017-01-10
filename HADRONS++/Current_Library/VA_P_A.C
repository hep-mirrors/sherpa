#include "HADRONS++/Current_Library/VA_P_A.H"
#include "METOOLS/Main/Polarization_Tools.H"


using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;

#include "HADRONS++/Current_Library/VA_P_A_PoleFit.C"

//Results are not correct.  # TODOFS

namespace HADRONS { namespace VA_P_A_FFs {
  FormFactor_Base::~FormFactor_Base()
  {
  }
class NoFF : public FormFactor_Base {
public:
  NoFF(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       std::vector<int>& i) : 
    FormFactor_Base(model, masses, flavs, i) {}
  void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) {
    m_V0=m_V1=m_V2=m_A=0.0;
    m_calced=true;
  }
};
} }

void VA_P_A::SetModelParameters( struct GeneralModel model )
{
  double Vxx(1.0);
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  
  if (kf0==kf_B_c) {
    if (kf1==10533)
      Vxx=Tools::Vcs;
    else if (kf1==10513)
      Vxx=Tools::Vcd;
  }
  
  m_Vxx = model("Vxx", Vxx);
  m_cV  = model("cV",1.0);
  switch( int(model("FORM_FACTOR", 1)+0.5) ) {
  case 0:
    p_ff = new VA_P_A_FFs::PoleFit(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using Pole Fit hep-ph 1607.00622v1 "<<m_name<<std::endl;
    break;
  
  default:
    msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}

void VA_P_A::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p0 = moms[p_i[0]], p1 = moms[p_i[1]];
  double m0 = p_masses[0], m1 = p_masses[1];
  
  Vec4D q = p0 - p1;
  double q2 = q.Abs2();
  p_ff->CalcFFs(p0,p1);
  Complex i = Complex(0.0,1.0);
  

  Polarization_Vector pol(p1, sqr(m1));
 
  for(int h_had=0;h_had<3;h_had++) {
    Vec4C eps = pol[h_had]; 
    
    Vec4C current(0.0,0.0,0.0,0.0);
    
    
    current +=-2.0*m1*(conj(eps)*q)/q2*q*p_ff->V0();
    current += -(m0+m1)*p_ff->V1()*(conj(eps)-(conj(eps)*q)/q2*q);
    current += +(conj(eps)*(p0+p1))/(m0+m1)*p_ff->V2()*((p0+p1)-(sqr(m0)-sqr(m1))/q2*q);
    current += +i/(m0-m1)*cross(conj(eps),(p0+p1),q)*p_ff->A();
    
    /*
    double V3=(m0+m1)/(2.0*m1)*p_ff->V1() - (m0-m1)/(2.0*m1)*p_ff->V2();
    current += p_ff->V1()*(m0+m1) * conj(eps);
    current += -p_ff->V2()*(conj(eps)*q)*((p0+p1)/(m0+m1));
    current += -(V3-p_ff->V0())*2.0*m1*(conj(eps)*q)*(q/q2);
    current += i*p_ff->A()*2.0/(m0+m1) * cross(conj(eps),p0,p1);
    */
   
    Insert(m_Vxx*m_cV*(m_anti?conj(current):current), h_had);
    
    
    
    }

}

DEFINE_CURRENT_GETTER(VA_P_A,"VA_P_A")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_P_A>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Results wrong"  <<std::endl;
}
