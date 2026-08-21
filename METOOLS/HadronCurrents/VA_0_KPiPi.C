#include "METOOLS/HadronCurrents/VA_0_KPiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// KpiPi isospin family (Finkemeier-Mirkes hep-ph/9503474): K^-pi^-pi^+,
// pi^-K0bar pi^0. Same momentum convention and Calc() structure as
// VA_0_KPiK - see that file's header comment.
//
///////////////////////////////////////////////////////////////////////////

VA_0_KPiPi::VA_0_KPiPi(const ATOOLS::Flavour_Vector& flavs,
			const std::vector<int>& indices,
			const std::string& name) :
  Current_Base(flavs, indices, name),
  m_norm(1.),
  p_f1(NULL), p_f2(NULL), p_f3(NULL), p_fS(NULL)
{
  msg_Out()<<METHOD<<"(N_f = "<<m_flavs.size()<<"):\n";
  for (size_t i=0;i<p_i.size();i++) {
    msg_Out()<<"    *  i = "<<i<<": "<<p_i[i]<<"  --> "
	     <<m_flavs[p_i[i]]<<".\n";
  }
}

VA_0_KPiPi::~VA_0_KPiPi() {
  if (p_f1) { delete p_f1; p_f1 = NULL; }
  if (p_f2) { delete p_f2; p_f2 = NULL; }
  if (p_f3) { delete p_f3; p_f3 = NULL; }
  if (p_fS) { delete p_fS; p_fS = NULL; }
}

void VA_0_KPiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p1    = moms[p_i[0]],  p2 = moms[p_i[1]],  p3 = moms[p_i[2]];
  Vec4D q     = p1+p2+p3,    dq13 = p1-p3,       dq23 = p2-p3;
  double s123 = q.Abs2();
  double Qq13 = q*dq13/s123, Qq23 = q*dq23/s123;
  Vec4D  v1   = dq13-Qq13*q,   v2 = dq23-Qq23*q;
  Vec4C  v4   = Vec4C(cross(p1,p2,p3));
  Complex F1  = (*p_f1)(moms), F2 = (*p_f2)(moms);
  Complex F3  = (*p_f3)(moms), FS = (*p_fS)(moms);
  Insert( m_norm * (F1*v1 + F2*v2 + F3*q + FS*v4), 0);
}

void VA_0_KPiPi::SetModelParameters(struct GeneralModel model) {
  map<string,double> pmap;
  FF_Parameters params(ff_model(model[string("FORM_FACTOR")]),
		       m_flavs,p_i,pmap,"",&model);
  params.m_name = "F1_0_KPiPi";
  p_f1 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F2_0_KPiPi";
  p_f2 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F3_0_KPiPi";
  p_f3 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "FS_0_KPiPi";
  p_fS = FF_Getter::GetObject("FF_0_PPP",params);
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_KPiPi,"VA_0_KPiPi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_KPiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow K \\pi \\pi $ (KpiPi isospin family) \n\n"
    <<"Order: 0 = q1, 1 = q2, 2 = q3 in the per-channel convention of \n"
    <<"Finkemeier & Mirkes, hep-ph/9503474: \n"
    <<"  K^- pi^- pi^+,   pi^- K0bar pi^0 \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 100 :} Kuehn-Santamaria-style \n"
    <<"      (Finkemeier-Mirkes hep-ph/9503474, Tabs.I/II) \n"
    <<"  \\end{itemize} \n"
    <<"Status: needs K1(1270)/K1(1400)/K*(1714) lineshapes - see \n"
    <<"FF_0_PPP.C for details on what is missing. \n\n"
    <<"Reference: https://arxiv.org/abs/hep-ph/9503474 \n"
    <<std::endl;
}
