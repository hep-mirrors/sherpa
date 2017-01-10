#include "HADRONS++/Current_Library/VA_P_V.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;

#include "HADRONS++/Current_Library/VA_P_V_ISGW.C"
#include "HADRONS++/Current_Library/VA_P_V_ISGW2.C"
#include "HADRONS++/Current_Library/VA_P_V_HQET.C"
#include "HADRONS++/Current_Library/VA_P_V_HQET2.C"
#include "HADRONS++/Current_Library/VA_P_V_SumRules3.C"
#include "HADRONS++/Current_Library/VA_P_V_PoleFit.C"
#include "HADRONS++/Current_Library/VA_P_V_hepph0007169.C"
#include "HADRONS++/Current_Library/VA_P_V_PoleFit2.C"
#include "HADRONS++/Current_Library/VA_P_V_hepph160208918.C"

namespace HADRONS { namespace VA_P_V_FFs {
    FormFactor_Base::~FormFactor_Base()
    {
    }
  class NoFF : public FormFactor_Base {
  public:
    NoFF(GeneralModel model, double* masses, const Flavour_Vector& flavs,
         const std::vector<int>& i) : 
      FormFactor_Base(model,masses, flavs, i) {}
    void CalcFFs( ATOOLS::Vec4D pB, ATOOLS::Vec4D pS ) {
      m_A0=m_A1=m_A2=m_A3=m_V=1.0;
      m_calced=true;
    }
  };
} }

void VA_P_V::SetModelParameters( struct GeneralModel model )
{
  double Vxx(1.0);
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B || kf0==kf_B_plus || kf0==kf_B_s) {
    if (kf1==kf_D_star_2007 || kf1==kf_D_star_2010_plus ||
        kf1==kf_D_s_star_plus || kf1==kf_D_1_2420_plus || kf1==kf_D_1_2420 ||
        kf1==kf_D_s1_2536_plus || kf1==kf_D_1_H_plus || kf1==kf_D_1_H ||
        kf1==kf_D_s1_H)
      Vxx=Tools::Vcb;
    else if (kf1==kf_rho_770 || kf1==kf_rho_770_plus || kf1==kf_omega_782 || 
             kf1==kf_b_1_1235 || kf1==kf_b_1_1235_plus || kf1==kf_h_1_1170 ||
             kf1==kf_a_1_1260 || kf1==kf_a_1_1260_plus || kf1==kf_f_1_1285 ||
             kf1==kf_h_1_1380)
      Vxx=Tools::Vub;
  }
  else if (kf0==kf_D || kf0==kf_D_plus) {
    if (kf1==kf_K_star_892 || kf1==kf_K_star_892_plus ||
        kf1==kf_K_1_1270 || kf1==kf_K_1_1270_plus ||
        kf1==kf_K_1_1400 || kf1==kf_K_1_1400_plus)
      Vxx=Tools::Vcs;
    else if (kf1==kf_rho_770 || kf1==kf_rho_770_plus || kf1==kf_omega_782 || 
             kf1==kf_b_1_1235 || kf1==kf_b_1_1235_plus || kf1==kf_h_1_1170 ||
             kf1==kf_a_1_1260 || kf1==kf_a_1_1260_plus || kf1==kf_f_1_1285 ||
             kf1==kf_h_1_1380)
      Vxx=Tools::Vcd;
  }
  else if (kf0==kf_B_c) {
    
    if (kf1==kf_J_psi_1S || kf1==100443)
      Vxx=Tools::Vcb;
    else if (kf1==533)       
      Vxx=Tools::Vcs;
    else if (kf1==kf_B_star)
      Vxx=Tools::Vcd;
  }
  else if (kf0==kf_D_s_plus) {
    if (kf1==kf_phi_1020)   
      Vxx=Tools::Vcs;
    else if (kf1==kf_K_star_892)
      Vxx=Tools::Vcd;
  }

  m_Vxx = model("Vxx", Vxx);
  m_cV  = model("cV",1.0);
  switch( int(model("FORM_FACTOR", 1)+0.5) ) {
  case 0:
    p_ff = new VA_P_V_FFs::NoFF(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using no form factor model for "<<m_name<<std::endl;
    break;
  case 1:
    p_ff = new VA_P_V_FFs::ISGW(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new VA_P_V_FFs::ISGW2(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using ISGW2 form factor model for "<<m_name<<std::endl;
    break;
  case 3:
    p_ff = new VA_P_V_FFs::HQET(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using HQET form factor model for "<<m_name<<std::endl;
    break;
  case 4:
    p_ff = new VA_P_V_FFs::HQET2(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using HQET2 form factor model for "<<m_name<<std::endl;
    break;
  case 5:
    p_ff = new VA_P_V_FFs::PoleFit(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using PoleFit form factor model for "<<m_name<<std::endl;
    break;
  case 6:
    p_ff = new VA_P_V_FFs::hepph0007169(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using hepph0007169 form factor model for "<<m_name<<std::endl;
    break;
  case 8:
    p_ff = new VA_P_V_FFs::SumRules3(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using SumRules3 form factor model for "<<m_name<<std::endl;
    break;
  case 9:
    p_ff = new VA_P_V_FFs::PoleFit2(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using PoleFit 2 hep-ph 1607.00622v1 for "<<m_name<<std::endl;
    break;
case 10:
    p_ff = new VA_P_V_FFs::hepph160208918(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    hepph160208918 "<<m_name<<std::endl;
    break;
  default:
    msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    Abort();
  }
}

void VA_P_V::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p0 = moms[p_i[0]], p1 = moms[p_i[1]];
  double m0 = p_masses[0], m1 = p_masses[1];
  
  Vec4D q = p0 - p1;
  double q2 = q.Abs2();
  p_ff->CalcFFs(p0,p1);
  Complex i = Complex(0.0,1.0);

  double fac_costhl(1.0);
  if(m_flavs[p_i[0]].Kfcode()==kf_D_plus ||
     m_flavs[p_i[0]].Kfcode()==kf_D ||
     m_flavs[p_i[0]].Kfcode()==kf_D_s_plus)
    fac_costhl=-1.0;

  Polarization_Vector pol(p1, sqr(m1));
  for(int h_had=0;h_had<3;h_had++) {
    Vec4C eps = pol[h_had];

    Vec4C current(0.0,0.0,0.0,0.0);

    if(p_ff->A1()!=0.0)
      current += -i*p_ff->A1()*(m0+m1) * conj(eps);
    if(p_ff->A2()!=0.0)
      current += i*p_ff->A2()*(conj(eps)*q)*((p0+p1)/(m0+m1));
    if(p_ff->A3()-p_ff->A0()!=0.0)
      current += i*(p_ff->A3()-p_ff->A0())*2.0*m1*(conj(eps)*q)*(q/q2);
    if(p_ff->V()!=0.0)
      current += fac_costhl*p_ff->V()*2.0/(m0+m1) * cross(conj(eps),p0,p1);
    
    Insert(m_Vxx*m_cV*(m_anti?conj(current):current), h_had);
  }
}

DEFINE_CURRENT_GETTER(VA_P_V,"VA_P_V")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_P_V>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow D^* (l \\nu_l) $ \n\n"
    <<"Order: 0 = (Pseudo)Scalar, 1 = Vector \n\n"
    <<"\\begin{eqnarray*} \\langle V(p_1,\\epsilon) | (V-A)_\\mu | P(p_0) \\rangle & = \n"
    <<"  & - i A_1(q^2) (m_0+m_1) \\epsilon^*_\\mu \\\\ \n"
    <<"& & + i A_2(q^2) \\frac{\\epsilon^* \\cdot q}{m_0+m_1} (p_0+p_1)_\\mu \\\\ \n"
    <<"& & + i (A_3(q^2)-A_0(q^2)) \\frac{2 m_1 \\epsilon^* \\cdot q}{q^2} q_\\mu \\\\ \n"
    <<"& & + V(q^2) \\frac{2}{m_0+m_1} \\varepsilon_{\\mu\\nu\\rho\\sigma} \\epsilon^{*\\nu} p_0^\\rho p_1^\\sigma \\end{eqnarray*} \n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 0 :} no form factor \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} ISGW http://www.slac.stanford.edu/spires/find/hep/www?j=PHRVA,D39,799 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} ISGW2 hep-ph/9503486 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 3 :} HQET hep-ph/9306320 and hep-ph/9508250 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 4 :} HQET2  \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 5 :} Sum Rules Ball/Zwicky hep-ph/0412079 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 6 :} Sum Rules hep-ph/9901395, hep-ph/9811259; with A vs. V \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 7 :} Sum Rules hep-ph/9811259 with A vs. V for a1 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 8 :} Sum Rules parametrization as in hep-ph/0608264 \n"
    <<"  \\end{itemize} \n"
    <<std::endl;
}
