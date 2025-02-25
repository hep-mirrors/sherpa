#include "HADRONS++/Current_Library/VA_0_PiPiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// Form factors for pi pi pi, K K, K pi final-state currents from:
// - KS, 1 (Kuehn-Santamaria model):
//   * pi pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
// - none, 0: no form factor
//
///////////////////////////////////////////////////////////////////////////

VA_0_PiPiPi::VA_0_PiPiPi(const ATOOLS::Flavour_Vector& flavs,
		       const std::vector<int>& indices,
		       const std::string& name) :
  Current_Base(flavs, indices, name),
  m_restype(resonance_type::running),
  m_PSmode(PSmode::unknown),
  m_ffmodel(ffmodel::KS),
  m_global(1.), m_deltaM2(0.),
  m_fpi(0.13041), m_m2_pi(sqr(Flavour(kf_pi_plus).Mass(true))),
  m_mu2(sqr(Flavour(kf_rho_770_plus).Mass(true)))
{
  if (m_flavs[p_i[0]].Kfcode()==kf_pi &&
      m_flavs[p_i[1]].Kfcode()==kf_pi &&
      m_flavs[p_i[2]].Kfcode()==kf_pi_plus)             m_PSmode = PSmode::pipipi_p00;
  msg_Out()<<METHOD<<"("<<m_flavs[p_i[0]].Kfcode()<<", "
	   <<m_flavs[p_i[1]].Kfcode()<<", "<<m_flavs[p_i[2]].Kfcode()<<")\n";
  if (m_PSmode==PSmode::unknown)
    THROW(fatal_error,"Current called for illegal flavour combination.");
}

VA_0_PiPiPi::~VA_0_PiPiPi() {
}

void VA_0_PiPiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D  q012 = moms[p_i[2]]+moms[p_i[1]]+moms[p_i[0]];
  double s012 = q012.Abs2();
}

void VA_0_PiPiPi::SetModelParameters(struct GeneralModel model)
{
  m_fpi     = model("fpi", 0.13041 );
  m_ffmodel = ffmodel(model("FORM_FACTOR",1));
  if (m_ffmodel==ffmodel::KS) {
  }
  // global pre-factor: 1/sqrt(2) for pi_0 wave-function, V_ud for the
  // quark-level coupling producing a rho (or rho-resonance), 1/sqrt(2)
  // for the overall normalisation.
  double iso = ( (m_flavs[p_i[0]].Kfcode()==kf_pi_plus ||
		  m_flavs[p_i[1]].Kfcode()==kf_pi_plus) ? sqrt(0.5) : 1.);
  double CKM = ( (m_PSmode==PSmode::pipipi_pmp || m_PSmode==PSmode::pipipi_p00 ) ?
		 model("Vud", Tools::Vud) : 0. );
  switch (int(model("RUNNING_WIDTH",1))) {
  case 10: m_restype = resonance_type::bespoke; break; 
  case 1:  m_restype = resonance_type::running; break;
  case 0:  m_restype = resonance_type::fixed;   break;
  default: m_restype = resonance_type::running; break;
  }
  if (m_ffmodel!=ffmodel::none) { 
  }
}



// need to update the references once we have the tau's under control
DEFINE_CURRENT_GETTER(HADRONS::VA_0_PiPiPi,"VA_0_PiPiPi")

void ATOOLS::Getter<HADRONS::Current_Base,
		    HADRONS::ME_Parameters,HADRONS::VA_0_PiPiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
