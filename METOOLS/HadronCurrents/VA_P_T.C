#include "METOOLS/HadronCurrents/VA_P_T.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_T::VA_P_T(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_VT(NULL), p_A0T(NULL), p_A1T(NULL), p_A2T(NULL),
  m_A1T_0(1.), m_rVT(1.), m_r2T(1.), m_useratios(true)
{}

void VA_P_T::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // Identical to VA_P_V up to eps* -> eps_T = eps*^{mu nu} p_nu / M,
  // which is precisely what the note says, and now literally what the
  // code does: same kernel, different polarization object.
  const Vec4D p = moms[p_i[0]], pX = moms[p_i[1]];
  std::vector<Vec4C> epsT;
  EffectiveTensorPolarizations(pX,p,epsT);
  for (size_t h=0;h<5;h++) {
    const Vec4C J = m_norm*VectorKernel(p,pX,epsT[h],
                                        p_VT,p_A0T,p_A1T,p_A2T,m_epssign);
    Insert(anti?conj(J):J,h);
  }
}

void VA_P_T::SetModelParameters(struct GeneralModel model) {
  CheckArity(1);
  ReadCommonParameters(model);

  p_VT  = MakeFF("VT", model);
  p_A0T = MakeFF("A0T",model);
  p_A1T = MakeFF("A1T",model);
  p_A2T = MakeFF("A2T",model);

  m_useratios = (model("USE_RATIOS",1.)>0.5);
  m_A1T_0 = model("A1T_0",1.);
  m_rVT   = model("rVT",  1.);   // MODEL ASSUMPTION, not a measurement
  m_r2T   = model("r2T",  1.);   // MODEL ASSUMPTION, not a measurement
  if (m_useratios) {
    const double v0  = p_VT ->Value(0.,m_mX);
    const double a10 = p_A1T->Value(0.,m_mX);
    const double a20 = p_A2T->Value(0.,m_mX);
    if (std::abs(a10)>1.e-300) p_A1T->SetF0(p_A1T->GetF0()*m_A1T_0/a10);
    if (std::abs(v0) >1.e-300) p_VT ->SetF0(p_VT ->GetF0()*m_rVT*m_A1T_0/v0);
    if (std::abs(a20)>1.e-300) p_A2T->SetF0(p_A2T->GetF0()*m_r2T*m_A1T_0/a20);
  }
  PrintCommonParameters();
  msg_Tracking()<<"###   rVT = "<<m_rVT<<", r2T = "<<m_r2T
                <<" - MODEL ASSUMPTIONS (BESIII fixed both to 1), "
                <<"not fitted form-factor ratios.\n";
}

// Tag is the plain "VA_P_T". The legacy HADRONS++ Current_Library
// class of the same name has been retired (its whole directory is
// commented out of the build), so there is exactly one getter
// answering to this tag. Both libraries registered into the same
//   Getter<METOOLS::Current_Base, METOOLS::ME_Parameters>
// registry, which is why they could not coexist: the table resolved
// to whichever registered first, every key in the YAML missed, and
// the legacy current fell back to ISGW without an error. If the old
// library is ever re-enabled, this tag must be suffixed again.
DEFINE_CURRENT_GETTER(METOOLS::VA_P_T,"VA_P_T")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_T>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D \\rightarrow T \\ell \\nu $, $T=K^*_2(1430)$ \n\n"
    <<"Order: 0 = decaying $D$, 1 = recoiling tensor \n\n"
    <<"Form factors: {\\tt VT}, {\\tt A0T}, {\\tt A1T}, {\\tt A2T}, in \n"
    <<"the $D\\to V$ basis with $\\eps^*\\to\\eps_T=\\eps^{*\\mu\\nu}p_\\nu/M$. \n"
    <<"Normalisation: {\\tt A1T\\_0}, {\\tt rVT}, {\\tt r2T} - the \n"
    <<"latter two default to 1, reproducing the BESIII 2026 D-wave \n"
    <<"MODEL ASSUMPTION rather than any measurement. \n\n"
    <<"Reference: charm\\_semileptonic.tex, Sec. D to T. \n"
    <<std::endl;
}
