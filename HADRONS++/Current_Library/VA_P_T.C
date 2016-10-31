#include "HADRONS++/Current_Library/VA_P_T.H"
#include "METOOLS/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;

#include "HADRONS++/Current_Library/VA_P_T_ISGW.C"
#include "HADRONS++/Current_Library/VA_P_T_ISGW2.C"
#include "ATOOLS/Org/My_MPI.H"

namespace HADRONS { namespace VA_P_T_FFs {
    FormFactor_Base::~FormFactor_Base()
    {
    }
  class NoFF : public FormFactor_Base {
  public:
    NoFF(GeneralModel model, double* masses, const Flavour_Vector& flavs,
         const std::vector<int>& i) :
      FormFactor_Base(model,masses, flavs, i) {}
    void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) {
      m_h=m_k=m_bplus=m_bminus=1.0;
      m_calced=true;
    }
  };
} }

void VA_P_T::SetModelParameters(struct GeneralModel model)
{
  double Vxx(1.0);
  kf_code kf0=m_flavs[p_i[0]].Kfcode();
  kf_code kf1=m_flavs[p_i[1]].Kfcode();
  if (kf0==kf_B || kf0==kf_B_plus || kf0==kf_B_s || kf0==kf_B_c) {
    if (kf1==kf_D_2_star_2460 || kf1==kf_D_2_star_2460_plus ||
        kf1==kf_D_s2_star_2573)
      Vxx=Tools::Vcb;
    else if (kf1==kf_f_2_1270 || kf1==kf_f_2_prime_1525 ||
             kf1==kf_a_2_1320 || kf1==kf_a_2_1320_plus)
      Vxx=Tools::Vub;
  }
  else if (kf0==kf_D || kf0==kf_D_plus || kf0==kf_D_s_plus) {
    if (kf1==kf_K_2_star_1430 || kf1==kf_K_2_star_1430_plus)
      Vxx=Tools::Vcs;
    else if (kf1==kf_f_2_1270 || kf1==kf_f_2_prime_1525 ||
             kf1==kf_a_2_1320 || kf1==kf_a_2_1320_plus)
      Vxx=Tools::Vcd;
  }
  m_Vxx = model("Vxx", Vxx);
  switch( int(model("FORM_FACTOR", 2)+0.5) ) {
  case 0:
    p_ff = new VA_P_T_FFs::NoFF(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using no form factor model for "<<m_name<<std::endl;
    break;
  case 1:
    p_ff = new VA_P_T_FFs::ISGW(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new VA_P_T_FFs::ISGW2(model, p_masses, m_flavs, p_i);
    msg_Tracking()<<"    Using ISGW2 form factor model for "<<m_name<<std::endl;
    break;
  default:
    msg_Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    Abort();
  }
}

void VA_P_T::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p0 = moms[p_i[0]];
  Vec4D p1 = moms[p_i[1]];
  Vec4D q = p0 - p1;
  double m1 = p_masses[1];

  // Get formfactors
  p_ff->CalcFFs(p0,p1);
  double h = p_ff->h();
  double k = p_ff->k();
  double bplus = p_ff->bplus();
  double bminus = p_ff->bminus();

  double fac_costhl(1.0);
  if(m_flavs[p_i[0]].Kfcode()==kf_D_plus ||
     m_flavs[p_i[0]].Kfcode()==kf_D ||
     m_flavs[p_i[0]].Kfcode()==kf_D_s_plus)
    fac_costhl=-1.0;

  Vec4C amp;
  Polarization_Tensor pol(p1, sqr(m1));
  for( int h_had=0; h_had<5; h_had++) {
    CMatrix eps = pol[h_had];
    Complex i = Complex(0.0,1.0);
    amp = m_Vxx*( fac_costhl*h*i*cross((eps.Conjugate()*p0),p0+p1,p0-p1)
                               - k*(eps.Conjugate()*p0) - bplus*(eps.Conjugate()*p0*p0)*(p0+p1)
                               - bminus*(eps.Conjugate()*p0*p0)*(p0-p1));
    Insert(m_anti?conj(amp):amp, h_had );
  }
}

DEFINE_CURRENT_GETTER(VA_P_T,"VA_P_T")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_P_T>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow D^* (l \\nu_l) $ \n\n"
    <<"Order: 0 = (Pseudo)Scalar, 1 = Tensor \n\n"
    <<"\\begin{eqnarray*} \\langle T(p_1,\\epsilon) | (V-A)_\\mu | P(p_0) \\rangle & = \n"
    <<"  & ih(q^2) \\varepsilon_{\\mu\\nu\\lambda\\rho} \\epsilon^{*\\nu\\alpha} p_{0\\alpha} (p_0+p_1)^{\\lambda} (p_0-p_1)^{\\rho} \\\\ \n"
    <<"& & -k(q^2) \\epsilon^*_{\\mu\\nu} p_0^{\\nu} \\\\ \n"
    <<"& & -b_+(q^2) \\varepsilon^*_{\\alpha\\beta} p_0^\\alpha p_0^\\beta (p_0+p_1)_{\\mu} \\\\ \n"
    <<"& & -b_-(q^2) \\varepsilon^*_{\\alpha\\beta} p_0^\\alpha p_0^\\beta (p_0-p_1)_{\\mu} \n"
    <<"\\end{eqnarray*} \n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 0 :} no form factor \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} ISGW http://www.slac.stanford.edu/spires/find/hep/www?j=PHRVA,D39,799 \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} ISGW2 hep-ph/9503486 \n"
    <<"  \\end{itemize} \n"
    <<std::endl;
}
