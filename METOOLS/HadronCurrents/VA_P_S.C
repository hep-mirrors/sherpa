#include "METOOLS/HadronCurrents/VA_P_S.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_S::VA_P_S(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_fplus(NULL), p_fzero(NULL), m_tiefzero(true), m_phase(0.5*M_PI)
{}

void VA_P_S::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  const Vec4D p = moms[p_i[0]], pX = moms[p_i[1]];
  const Complex phase(cos(m_phase),sin(m_phase));
  const Vec4C J = m_norm*phase*ScalarKernel(p,pX,p_fplus,p_fzero);
  Insert(anti?conj(J):J,0);
}

void VA_P_S::SetModelParameters(struct GeneralModel model) {
  CheckArity(1);
  ReadCommonParameters(model);
  m_phase = model("SCALAR_PHASE",0.5*M_PI);
  p_fplus = MakeFF("Fplus",model);
  p_fzero = MakeFF("Fzero",model);

  m_tiefzero = (model("TIE_FZERO",1.)>0.5);
  if (m_tiefzero) {
    const double fp0 = p_fplus->Value(0.,m_mX);
    const double f00 = p_fzero->Value(0.,m_mX);
    if (std::abs(f00)>1.e-300)
      p_fzero->SetF0(p_fzero->GetF0()*fp0/f00);
  }
  PrintCommonParameters();
  msg_Tracking()<<"###   f+^S(0) = "<<p_fplus->Value(0.,m_mX)
                <<", f0^S(0) = "<<p_fzero->Value(0.,m_mX)
                <<", global phase = "<<m_phase<<" rad\n"
                <<"###   NOTE: electron-mode fits constrain f+^S only; "
                <<"the f0^S shape here is a model choice.\n";
}

// Tag is the plain "VA_P_S". The legacy HADRONS++ Current_Library
// class of the same name has been retired (its whole directory is
// commented out of the build), so there is exactly one getter
// answering to this tag. Both libraries registered into the same
//   Getter<METOOLS::Current_Base, METOOLS::ME_Parameters>
// registry, which is why they could not coexist: the table resolved
// to whichever registered first, every key in the YAML missed, and
// the legacy current fell back to ISGW without an error. If the old
// library is ever re-enabled, this tag must be suffixed again.
DEFINE_CURRENT_GETTER(METOOLS::VA_P_S,"VA_P_S")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_S>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D \\rightarrow S \\ell \\nu $, "
    <<"$S=K^*_0(1430),f_0(980),f_0(500),a_0(980)$ \n\n"
    <<"Order: 0 = decaying $D$, 1 = recoiling scalar \n\n"
    <<"Form factors: {\\tt Fplus}, {\\tt Fzero} (i.e. $f^S_+,f^S_0$), \n"
    <<"same shape/parameter interface as {\\tt VA\\_P\\_P}. \n"
    <<"{\\tt SCALAR\\_PHASE} sets the global phase (default $\\pi/2$, \n"
    <<"the $+i$ of Eq.(Scurrent)); it is a convention, but must be \n"
    <<"kept consistent across every channel using the same state. \n\n"
    <<"Reference: charm\\_semileptonic.tex, Eq.(Scurrent). \n"
    <<std::endl;
}
