#include "METOOLS/HadronCurrents/VA_P_P.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_P::VA_P_P(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_fplus(NULL), p_fzero(NULL), m_tiefzero(true)
{}

void VA_P_P::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // No Levi-Civita structure, so the charge-conjugate current differs
  // only by a global phase; conj() is applied for uniformity.
  const Vec4D p = moms[p_i[0]], pX = moms[p_i[1]];
  const Vec4C J = m_norm*ScalarKernel(p,pX,p_fplus,p_fzero);
  Insert(anti?conj(J):J,0);
}

void VA_P_P::SetModelParameters(struct GeneralModel model) {
  CheckArity(1);
  ReadCommonParameters(model);
  p_fplus = MakeFF("Fplus",model);
  p_fzero = MakeFF("Fzero",model);

  // Impose f+(0)=f0(0) exactly rather than hoping two independently
  // supplied normalisations agree. The note is explicit that most
  // experimental fits determine f+ only and that one must NOT silently
  // set f0(q^2)=f+(q^2) - so the SHAPE of f0 stays whatever the YAML
  // selected (a 0^+ pole by default), only its normalisation is tied.
  m_tiefzero = (model("TIE_FZERO",1.)>0.5);
  if (m_tiefzero) {
    const double fp0 = p_fplus->Value(0.,m_mX);
    const double f00 = p_fzero->Value(0.,m_mX);
    if (std::abs(f00)>1.e-300)
      p_fzero->SetF0(p_fzero->GetF0()*fp0/f00);
    else
      msg_Error()<<"Error in "<<METHOD<<": f0 vanishes at q^2=0 for "
                 <<m_name<<", cannot impose f+(0)=f0(0).\n";
  }
  PrintCommonParameters();
  msg_Tracking()<<"###   f+(0) = "<<p_fplus->Value(0.,m_mX)
                <<", f0(0) = "<<p_fzero->Value(0.,m_mX)
                <<" (tied: "<<(m_tiefzero?"yes":"no")<<")\n";
}

// Tag is the plain "VA_P_P". The legacy HADRONS++ Current_Library
// class of the same name has been retired (its whole directory is
// commented out of the build), so there is exactly one getter
// answering to this tag. Both libraries registered into the same
//   Getter<METOOLS::Current_Base, METOOLS::ME_Parameters>
// registry, which is why they could not coexist: the table resolved
// to whichever registered first, every key in the YAML missed, and
// the legacy current fell back to ISGW without an error. If the old
// library is ever re-enabled, this tag must be suffixed again.
DEFINE_CURRENT_GETTER(METOOLS::VA_P_P,"VA_P_P")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_P>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D \\rightarrow P \\ell \\nu $, $P=K,\\pi,\\eta,\\eta'$ \n\n"
    <<"Order: 0 = decaying $D$, 1 = recoiling pseudoscalar \n\n"
    <<"Form factors: {\\tt Fplus}, {\\tt Fzero}. Each takes its own \n"
    <<"{\\tt <name>\\_SHAPE} from the {\\tt ffq2\\_shape} enum in \n"
    <<"{\\tt FF\\_P\\_X.H} plus the parameters of that shape. \n"
    <<"$f_0(0)=f_+(0)$ is imposed unless {\\tt TIE\\_FZERO: 0}. \n\n"
    <<"Common keys: {\\tt Vcq}, {\\tt ISOSPIN}, {\\tt CKM\\_DRESSED}, \n"
    <<"{\\tt DYNAMIC\\_MASS}, {\\tt Q2MIN}. \n\n"
    <<"Reference: charm\\_semileptonic.tex, Eq.(Pcurrent). \n"
    <<std::endl;
}
