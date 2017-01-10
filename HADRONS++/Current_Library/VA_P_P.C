#include "HADRONS++/Current_Library/VA_P_P.H"

using namespace HADRONS;
using namespace ATOOLS;

#include "ATOOLS/Org/My_MPI.H"
#include "HADRONS++/Current_Library/VA_P_P_ISGW.C"
#include "HADRONS++/Current_Library/VA_P_P_ISGW2.C"
#include "HADRONS++/Current_Library/VA_P_P_HQET.C"
#include "HADRONS++/Current_Library/VA_P_P_HQET2.C"
#include "HADRONS++/Current_Library/VA_P_P_PoleFit.C"
#include "HADRONS++/Current_Library/VA_P_P_Polynomial.C"
#include "HADRONS++/Current_Library/VA_P_P_BallZwicky.C"
#include "HADRONS++/Current_Library/VA_P_P_BGL.C"
#include "HADRONS++/Current_Library/VA_P_P_ISGW3.C"
#include "HADRONS++/Current_Library/VA_P_P_PoleFit2.C"
#include "HADRONS++/Current_Library/VA_P_P_hepph160208918.C"

namespace HADRONS { namespace VA_P_P_FFs {
  FormFactor_Base::~FormFactor_Base()
  {
  }
class NoFF : public FormFactor_Base {
public:
  NoFF(GeneralModel model, double* masses, const Flavour_Vector& flavs,
       std::vector<int>& i) : 
    FormFactor_Base(model, masses, flavs, i) {}
  void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) {
    m_fplus=1.0; m_f0=0.0;
    m_calced=true;
  }
};
} }

void VA_P_P::SetModelParameters( struct GeneralModel model )
{
  double Vxx(1.0);
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  
  if (kf0==kf_B || kf0==kf_B_plus || kf0==kf_B_s) {
    if (kf1==kf_D || kf1==kf_D_plus || kf1==kf_D_s_plus || 
        kf1==kf_D_0_star || kf1==kf_D_0_star_plus)
      Vxx=Tools::Vcb;
    else if (kf1==kf_pi || kf1==kf_pi_plus || kf1==kf_eta || 
             kf1==kf_eta_prime_958 || kf1==kf_a_0_980 ||
             kf1==kf_a_0_980_plus ||  kf1==kf_f_0_980 || kf1==kf_f_0_1370)
      Vxx=Tools::Vub;
  }
  else if (kf0==kf_D || kf0==kf_D_plus) {
    if (kf1==kf_K_plus || kf1==kf_K || kf1==kf_K_S || kf1==kf_K_L ||
        kf1==kf_K_0_star_1430 || kf1==kf_K_0_star_1430_plus)
      Vxx=Tools::Vcs;
    else if (kf1==kf_pi || kf1==kf_pi_plus || kf1==kf_eta || 
             kf1==kf_eta_prime_958 || kf1==kf_a_0_980 ||
             kf1==kf_a_0_980_plus ||  kf1==kf_f_0_980 || kf1==kf_f_0_1370)
      Vxx=Tools::Vcd;
  }
  else if (kf0==kf_K_plus || kf0==kf_K || kf0==kf_K_S || kf0==kf_K_L) {
    if (kf1==kf_pi || kf1==kf_pi_plus)
      Vxx=Tools::Vus;
  }
  else if (kf0==kf_B_c) {
    
    if (kf1==kf_eta_c_1S||kf1==100441)
      Vxx=Tools::Vcb;
    else if (kf1==kf_B_s)
      Vxx=Tools::Vcs;
    else if (kf1=kf_B)
      Vxx=Tools::Vcd;
    
    
  }
  else if (kf0==kf_D_s_plus) {
    
    if (kf1==kf_eta_prime_958 || kf1==kf_eta)
      Vxx=Tools::Vcs;
    else if (kf1==kf_K || kf1==kf_K_S || kf1==kf_K_L)
      Vxx=Tools::Vcd;
  }
  m_Vxx = model("Vxx", Vxx);
  
  switch( int(model("FORM_FACTOR", 1)+0.5) ) {
  case 0:
    p_ff = new VA_P_P_FFs::NoFF(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using no form factor model for "<<m_name<<std::endl;
    break;
  case 1:
    p_ff = new VA_P_P_FFs::ISGW(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new VA_P_P_FFs::ISGW2(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using ISGW2 form factor model for "<<m_name<<std::endl;
    break;
  case 3:
    p_ff = new VA_P_P_FFs::HQET(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using HQET form factor model for "<<m_name<<std::endl;
    break;
  case 4:
    p_ff = new VA_P_P_FFs::HQET2(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using HQET2 form factor model for "<<m_name<<std::endl;
    break;
  case 5:
    p_ff = new VA_P_P_FFs::PoleFit(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using PoleFit form factor model for "<<m_name<<std::endl;
    break;
  case 7:
    p_ff = new VA_P_P_FFs::Polynomial(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using Polynomial form factor model for "<<m_name<<std::endl;
    break;
  case 8:
    p_ff = new VA_P_P_FFs::BallZwicky(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using BallZwicky form factor model for "<<m_name<<std::endl;
    break;
  case 9:
    p_ff = new VA_P_P_FFs::BGL(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using BGL form factor model for "<<m_name<<std::endl;
    break;   
  case 10:
    p_ff = new VA_P_P_FFs::ISGW3(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using ISGW3 form factor model for "<<m_name<<std::endl;
    break;
  case 11: 
    p_ff= new VA_P_P_FFs::PoleFit2(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<" Using PoleFit2"<<m_name<<std::endl;
    break;
  case 12:  
    p_ff= new VA_P_P_FFs::hepph160208918(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"hep-ph1602.08918"<<m_name<<std::endl;
    break;
  default:
    msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    Abort();
  }
}

void VA_P_P::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p0 = moms[p_i[0]], p1 = moms[p_i[1]];
  double m0 = p_masses[0], m1 = p_masses[1];
  
  Vec4D q = p0 - p1;
  double q2 = q.Abs2();
  p_ff->CalcFFs(p0,p1);

  Vec4C current(0.0,0.0,0.0,0.0);
  if(p_ff->fplus()!=0.0)
    current += p_ff->fplus()*(p0+p1 - (sqr(m0)-sqr(m1))/q2*q);
  if(p_ff->f0()!=0.0)
    current += p_ff->f0()*((sqr(m0)-sqr(m1))/q2*q);

  Insert(m_Vxx*current, 0);

}

DEFINE_CURRENT_GETTER(VA_P_P,"VA_P_P")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_P_P>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow D (l \\nu_l) $ \n\n"
    <<"Order: 0 = (Pseudo)Scalar, 1 = (Pseudo)Scalar \n\n"
    <<"\\[ \\langle P(p_1) | (V-A)_\\mu | P(p_0) \\rangle = "
    <<"  f_+(q^2)\\left\\{ (p_0+p_1)_\\mu - \\frac{m_0^2-m_1^2}{q^2} q_\\mu \\right\\} + "
    <<"  f_0(q^2)\\left\\{ \\frac{m_0^2-m_1^2}{q^2} q_\\mu \\right\\} \\] \n \n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 0 :} no form factor \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} ISGW http://www.slac.stanford.edu/spires/find/hep/www?j=PHRVA,D39,799 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} ISGW2 arXiv:hep-ph/9503486 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 3 :} HQET hep-ph/9306320 and arXiv:hep-ph/9508250 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 4 :} HQET2  \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 5 :} PoleFit arXiv:hep-ph/0406232 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 6 :} PoleFit arXiv:hep-ph/0007169 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 7 :} Polynomial, e.g. linear $K \\rightarrow \\pi$ in PDG \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 8 :} Ball/Zwicky arXiv:0706.3628 \n"
    <<"  \\end{itemize} \n"
    <<std::endl;
}
