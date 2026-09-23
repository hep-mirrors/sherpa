#include "METOOLS/HadronCurrents/VA_P_V.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_V::VA_P_V(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_V(NULL), p_A0(NULL), p_A1(NULL), p_A2(NULL),
  m_A1_0(1.), m_rV(1.), m_r2(1.), m_useratios(true)
{}

void VA_P_V::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  // The whole Lorentz structure now lives in VA_P_X_Base::VectorKernel,
  // which VA_P_T and VA_P_PP call with a different "polarization".
  const Vec4D p = moms[p_i[0]], pX = moms[p_i[1]];
  std::vector<Vec4C> eps;
  PolarizationVectors(pX,eps);
  for (size_t h=0;h<3;h++) {
    const Vec4C J = m_norm*VectorKernel(p,pX,eps[h],
                                        p_V,p_A0,p_A1,p_A2,m_epssign);
    Insert(anti?conj(J):J,h);
  }
}

void VA_P_V::SetModelParameters(struct GeneralModel model) {
  CheckArity(1);
  ReadCommonParameters(model);

  p_V  = MakeFF("V", model);
  p_A0 = MakeFF("A0",model);
  p_A1 = MakeFF("A1",model);
  p_A2 = MakeFF("A2",model);

  // Experimental charm analyses quote {A_1(0), r_V, r_2} rather than
  // three independent normalisations, so that is the default input
  // route. Setting USE_RATIOS: 0 leaves whatever <name>_F0 the YAML
  // gave each form factor untouched (useful for a lattice/BGL input).
  m_useratios = (model("USE_RATIOS",1.)>0.5);
  m_A1_0 = model("A1_0",1.);
  m_rV   = model("rV",  1.);
  m_r2   = model("r2",  1.);
  if (m_useratios) {
    // Rescale each shape so that its value at q^2=0 is the requested
    // one; every shape in FF_P_X is linear in its stored F(0), so a
    // single multiplicative rescale is exact.
    const double v0 = p_V ->Value(0.,m_mX);
    const double a10= p_A1->Value(0.,m_mX);
    const double a20= p_A2->Value(0.,m_mX);
    if (std::abs(a10)>1.e-300) p_A1->SetF0(p_A1->GetF0()*m_A1_0/a10);
    if (std::abs(v0) >1.e-300) p_V ->SetF0(p_V ->GetF0()*m_rV*m_A1_0/v0);
    if (std::abs(a20)>1.e-300) p_A2->SetF0(p_A2->GetF0()*m_r2*m_A1_0/a20);
  }
  // A_0's stored normalisation is irrelevant - only its SHAPE is used,
  // the normalisation being fixed event by event by Eq.(A0constraint).
  PrintCommonParameters();
  const double A1_0 = p_A1->Value(0.,m_mX), A2_0 = p_A2->Value(0.,m_mX);
  msg_Tracking()<<"###   A1(0) = "<<A1_0<<", V(0) = "<<p_V->Value(0.,m_mX)
                <<", A2(0) = "<<A2_0<<"\n"
                <<"###   A0(0) [from Eq.(A0constraint), nominal m] = "
                <<((m_M+m_mX)*A1_0-(m_M-m_mX)*A2_0)/(2.*m_mX)
                <<" - recomputed per event with the generated s_X\n";
}

// Tag is the plain "VA_P_V". The legacy HADRONS++ Current_Library
// class of the same name has been retired (its whole directory is
// commented out of the build), so there is exactly one getter
// answering to this tag. Both libraries registered into the same
//   Getter<METOOLS::Current_Base, METOOLS::ME_Parameters>
// registry, which is why they could not coexist: the table resolved
// to whichever registered first, every key in the YAML missed, and
// the legacy current fell back to ISGW without an error. If the old
// library is ever re-enabled, this tag must be suffixed again.
DEFINE_CURRENT_GETTER(METOOLS::VA_P_V,"VA_P_V")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_V>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D \\rightarrow V \\ell \\nu $, "
    <<"$V=K^*(892),\\rho,\\omega,\\phi$ \n\n"
    <<"Order: 0 = decaying $D$, 1 = recoiling vector \n\n"
    <<"Form factors: {\\tt V}, {\\tt A0}, {\\tt A1}, {\\tt A2}, each \n"
    <<"with its own {\\tt <name>\\_SHAPE} and parameters. \n\n"
    <<"Normalisation: {\\tt A1\\_0}, {\\tt rV}$=V(0)/A_1(0)$, \n"
    <<"{\\tt r2}$=A_2(0)/A_1(0)$ (set {\\tt USE\\_RATIOS: 0} to use \n"
    <<"per-form-factor {\\tt \\_F0} values instead). $A_0(0)$ is never \n"
    <<"read: it is fixed by Eq.(A0constraint) event by event, using \n"
    <<"the generated $\\sqrt{s_X}$, so the $1/q^2$ structure stays \n"
    <<"exactly regular. Only the $A_0$ SHAPE is taken from the YAML. \n\n"
    <<"The strong decay and line shape of $V$ are NOT part of this \n"
    <<"current - they are handled by the ordinary decay chain. \n\n"
    <<"Reference: charm\\_semileptonic.tex, Eqs.(PVvector),(PVaxial). \n"
    <<std::endl;
}
