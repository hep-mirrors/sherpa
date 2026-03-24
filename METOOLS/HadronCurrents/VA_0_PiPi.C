#include "METOOLS/HadronCurrents/VA_0_PiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// Comment to be deleted once we have it all checked ....
// Please take a look at a summary of what is implemented in Tauola
// - https://arxiv.org/pdf/1509.09140
// - https://arxiv.org/abs/1609.04617
// and "harvest" the references within.
// Please test/debug the different form factor models (1 and 2).
//
// RUNNING_WIDTH in Decay.yaml translates into different width schemes,
// like fixed, running, Gounaris-Sakurai (and potentially other width schemes)
// for the resonances in
// METOOLS++/PS_Library/Resonance.[C,H]
//
// Named kf-codes (aliases) in ATOOLS/Phys/Flavour_Tags.H - essentially the
// PDG codes.
//
// Check overall norm in the end by calculating partial widths.
// 
///////////////////////////////////////////////////////////////////////////

VA_0_PiPi::VA_0_PiPi(const ATOOLS::Flavour_Vector& flavs,
		     const std::vector<int>& indices,
		     const std::string& name) :
  Current_Base(flavs, indices, name),
  m_norm(1.), m_deltaM2(sqr(m_flavs[p_i[1]].Mass(true))-
			sqr(m_flavs[p_i[0]].Mass(true))),
  p_fplus(NULL), p_fzero(NULL)
{ }

VA_0_PiPi::~VA_0_PiPi() {
  if (p_fplus) { delete p_fplus; p_fplus = NULL; }
  if (p_fzero) { delete p_fzero; p_fzero = NULL; }
}


void VA_0_PiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D   q = moms[p_i[1]]+moms[p_i[0]];
  double Q2 = q.Abs2();
  Vec4D   v = (moms[p_i[1]]-moms[p_i[0]]) - m_deltaM2/Q2*q;
  Complex V = (*p_fplus)(moms), S = (*p_fzero)(moms);
  Insert(m_norm * ( V * v + S * q ), 0);
}

void VA_0_PiPi::SetModelParameters(struct GeneralModel model) {
  FF_Parameters params(ff_model(model[string("FORM_FACTOR")]),
		       m_flavs,p_i,"F+_0_PP",&model);
  p_fplus = FF_Getter::GetObject("FF_0_PP",params);
  params.m_name = "F0_0_PP";
  p_fzero = FF_Getter::GetObject("FF_0_PP",params);
}

// need to update the references once we have the tau's under control
DEFINE_CURRENT_GETTER(METOOLS::VA_0_PiPi,"VA_0_PiPi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_PiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
