#include "METOOLS/HadronCurrents/VA_P_A.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

VA_P_A::VA_P_A(const ATOOLS::Flavour_Vector& flavs,
               const std::vector<int>& indices,
               const std::string& name) :
  VA_P_X_Base(flavs,indices,name),
  p_A(NULL), p_V0(NULL), p_V1(NULL), p_V2(NULL),
  m_V1_0(1.), m_rA(-0.112), m_rV2(-0.043),
  m_useratios(true), m_phase(-0.5*M_PI)
{}

void VA_P_A::Calc(const ATOOLS::Vec4D_Vector& moms, bool anti)
{
  const Vec4D p = moms[p_i[0]], pX = moms[p_i[1]];
  std::vector<Vec4C> eps;
  PolarizationVectors(pX,eps);
  for (size_t h=0;h<3;h++) {
    const Vec4C J = m_norm*AxialKernel(p,pX,eps[h],p_A,p_V0,p_V1,p_V2,
                                       m_epssign,m_phase);
    Insert(anti?conj(J):J,h);
  }
}

void VA_P_A::SetModelParameters(struct GeneralModel model) {
  CheckArity(1);
  ReadCommonParameters(model);
  m_phase = model("AXIAL_PHASE",-0.5*M_PI);

  p_A  = MakeFF("A", model);
  p_V0 = MakeFF("V0",model);
  p_V1 = MakeFF("V1",model);
  p_V2 = MakeFF("V2",model);

  m_useratios = (model("USE_RATIOS",1.)>0.5);
  m_V1_0 = model("V1_0",1.);
  m_rA   = model("rA",  -0.112);
  m_rV2  = model("rV2", -0.043);
  if (m_useratios) {
    const double v10 = p_V1->Value(0.,m_mX);
    const double v20 = p_V2->Value(0.,m_mX);
    const double a0  = p_A ->Value(0.,m_mX);
    if (std::abs(v10)>1.e-300) p_V1->SetF0(p_V1->GetF0()*m_V1_0/v10);
    if (std::abs(v20)>1.e-300) p_V2->SetF0(p_V2->GetF0()*m_rV2*m_V1_0/v20);
    if (std::abs(a0) >1.e-300) p_A ->SetF0(p_A ->GetF0()*m_rA *m_V1_0/a0);
  }
  PrintCommonParameters();
  msg_Tracking()<<"###   V1(0) = "<<p_V1->Value(0.,m_mX)
                <<", V2(0) = "<<p_V2->Value(0.,m_mX)
                <<", A(0) = "<<p_A->Value(0.,m_mX)<<"\n"
                <<"###   V0(0) fixed by Eq.(V0constraint) per event.\n"
                <<"###   SIGN WARNING: the measured rA, rV2 signs are "
                <<"convention-dependent - see VA_P_A.H.\n";
}

// Tag is the plain "VA_P_A". The legacy HADRONS++ Current_Library
// class of the same name has been retired (its whole directory is
// commented out of the build), so there is exactly one getter
// answering to this tag. Both libraries registered into the same
//   Getter<METOOLS::Current_Base, METOOLS::ME_Parameters>
// registry, which is why they could not coexist: the table resolved
// to whichever registered first, every key in the YAML missed, and
// the legacy current fell back to ISGW without an error. If the old
// library is ever re-enabled, this tag must be suffixed again.
DEFINE_CURRENT_GETTER(METOOLS::VA_P_A,"VA_P_A")

void ATOOLS::Getter<METOOLS::Current_Base,
                    METOOLS::ME_Parameters,METOOLS::VA_P_A>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ D \\rightarrow A \\ell \\nu $, $A=K_1(1270),K_1(1400)$ \n\n"
    <<"Order: 0 = decaying $D$, 1 = recoiling axial vector \n\n"
    <<"Form factors: {\\tt A}, {\\tt V0}, {\\tt V1}, {\\tt V2}. \n"
    <<"Normalisation: {\\tt V1\\_0}, {\\tt rA}$=A(0)/V_1(0)$, \n"
    <<"{\\tt rV2}$=V_2(0)/V_1(0)$. $V_0(0)$ follows from \n"
    <<"Eq.(V0constraint) event by event. \n\n"
    <<"Note that {\\tt rV2} is NOT the $r_V$ of {\\tt VA\\_P\\_V}; the \n"
    <<"BESIII paper calls it $r_V$, which is why it is renamed here. \n\n"
    <<"Reference: charm\\_semileptonic.tex, Eqs.(PAvector),(PAaxial); \n"
    <<"BESIII, Phys.Rev.Lett. 135, 091801 (2025), arXiv:2503.02196. \n"
    <<std::endl;
}
