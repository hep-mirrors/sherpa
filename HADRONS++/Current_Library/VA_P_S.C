#include "HADRONS++/Current_Library/VA_P_S.H"

using namespace HADRONS;
using namespace ATOOLS;

#include "HADRONS++/Current_Library/VA_P_S_PoleFit.C"


namespace HADRONS { namespace VA_P_S_FFs {
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

void VA_P_S::SetModelParameters( struct GeneralModel model )
{
  double Vxx(1.0);
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  

  if (kf0==kf_B_c) {
    

    if (kf1==10531)
      Vxx=Tools::Vcs;
    else if (kf1==10511)
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

  case 1:
    p_ff = new VA_P_S_FFs::PoleFit(model,p_masses,m_flavs,p_i);
    msg_Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  
  default:
    msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}

void VA_P_S::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p0 = moms[p_i[0]], p1 = moms[p_i[1]];
  double m0 = p_masses[0], m1 = p_masses[1];
  
  Vec4D q = p0 - p1;
  double q2 = q.Abs2();
  p_ff->CalcFFs(p0,p1);
  Complex i = Complex(0.0,1.0);

  Vec4C current(0.0,0.0,0.0,0.0);
  if(p_ff->fplus()!=0.0)
    current += -i*p_ff->fplus()*(p0+p1 - (sqr(m0)-sqr(m1))/q2*q);
  if(p_ff->f0()!=0.0) 
    current += -i*p_ff->f0()*((sqr(m0)-sqr(m1))/q2*q);
  
  
  Insert(m_Vxx*current, 0);

}

DEFINE_CURRENT_GETTER(VA_P_S,"VA_P_S")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_P_S>::
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
