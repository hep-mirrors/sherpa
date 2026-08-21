#include "METOOLS/HadronCurrents/VA_0_PiPiPi.H"
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

VA_0_PiPiPi::VA_0_PiPiPi(const ATOOLS::Flavour_Vector& flavs,
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

VA_0_PiPiPi::~VA_0_PiPiPi() {
  if (p_f1) { delete p_f1; p_f1 = NULL; }
  if (p_f2) { delete p_f2; p_f2 = NULL; }
  if (p_f3) { delete p_f3; p_f3 = NULL; }
  if (p_fS) { delete p_fS; p_fS = NULL; }
}


void VA_0_PiPiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p1    = moms[p_i[0]],  p2 = moms[p_i[1]],  p3 = moms[p_i[2]];
  // v1, v2 are built from (p1-p3) and (p2-p3) - i.e. the "odd" pion
  // p3 subtracted from each of the two Bose-symmetric pions p1,p2 -
  // matching the V1,V2 structures of Kuhn-Santamaria Z.Phys.C 48
  // (1990) 445, Eq.(3.3)/(3.4): V1=(q1-q3)_T, V2=(q2-q3)_T (see also
  // FF_0_PPP.H's index-convention comment: p1=q1,p2=q2,p3=q3). The
  // form factors F1,F2 computed below assume exactly this pairing -
  // F1 multiplies V1 and depends on s2=(p1+p3)^2 (the KS Brho(s2)),
  // F2 multiplies V2 and depends on the genuine s1=(p2+p3)^2 (the KS
  // Brho(s1)) - see FF_0_PPP.C::FF_KS. An earlier version of this
  // line built the vectors as (p2-p1)/(p3-p1) instead, which pairs
  // the form factors with the wrong momentum combinations and does
  // not reproduce the paper's amplitude (confirmed by explicit
  // numerical check against Eq.(3.3)-(3.5) using generic propagators).
  Vec4D q     = p1+p2+p3,    dq13 = p1-p3,       dq23 = p2-p3;
  double s123 = q.Abs2(),      s1 = (q-p1).Abs2(), s2 = (q-p2).Abs2();
  double Qq13 = q*dq13/s123, Qq23 = q*dq23/s123;
  Vec4D  v1   = dq13-Qq13*q,   v2 = dq23-Qq23*q;
  Vec4C  v4   = Vec4C(cross(p1,p2,p3));
  msg_Out()<<"\n"<<"*** "<<METHOD<<"("<<p_f1<<"/"<<p_f2<<"/"<<p_f3<<"/"<<p_fS<<"): "
  	   <<m_flavs[p_i[0]]<<" + "<<m_flavs[p_i[1]]<<" + "
  	   <<m_flavs[p_i[2]]<<": "
  	   <<"Qq13 = "<<Qq13<<", Qq23 = "<<Qq23<<" from s = "<<s123<<"\n"
  	   <<p1<<"("<<p1.Abs2()<<") + "<<p2<<"("<<p2.Abs2()<<") + "
  	   <<p3<<"("<<p3.Abs2()<<")\n";
  Complex F1  = (*p_f1)(moms), F2 = (*p_f2)(moms);
  Complex F3  = (*p_f3)(moms), FS = (*p_fS)(moms);
  msg_Out()<<F1<<", "<<F2<<", "<<F3<<","<<FS<<"\n"; 
  Insert( m_norm * (F1*v1 + F2*v2 + F3*q + FS*v4), 0);
}

void VA_0_PiPiPi::SetModelParameters(struct GeneralModel model) {
  // I think we should add the parameters here, maybe from some
  // generic parameter.yaml structure or so?  Will resolve during
  // the MCnet school while chatting with Chris & Eno
  map<string,double> pmap;
  FF_Parameters params(ff_model(model[string("FORM_FACTOR")]),
		       m_flavs,p_i,pmap,"",&model);
  msg_Out()<<METHOD<<"("<<m_flavs[p_i[0]]<<"/"<<m_flavs[p_i[1]]<<"/"<<m_flavs[p_i[2]]<<")\n";
  params.m_name = "F1_0_PPP";
  p_f1 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F2_0_PPP";
  p_f2 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F3_0_PPP";
  p_f3 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "FS_0_PPP";
  p_fS = FF_Getter::GetObject("FF_0_PPP",params);
}

// need to update the references once we have the tau's under control
DEFINE_CURRENT_GETTER(METOOLS::VA_0_PiPiPi,"VA_0_PPP")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_PiPiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^0$, 2 = $\\pi^\\pm$ or "
    <<"0 = $\\pi^\\pm$, 1 = $\\pi^\\pm$, 2 = $\\pi^\\mp$\n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 100 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 200 :} Resonance Chiral Theory \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 201 :} RChL2012 (arXiv:1203.3955 "
    <<"Sec.2.1, Eqs.4-11, with the sigma-meson extension of "
    <<"arXiv:1310.1053; a1 width uses the polynomial fit Eq.(17) of "
    <<"that paper, not the full double phase-space integral) \n"
    <<"  \\end{itemize} \n"
    <<"References: \n"
    <<"  https://arxiv.org/pdf/1203.3955 (Sec.2.1, RChL2012 3pi current) \n"
    <<"  https://arxiv.org/abs/1310.1053 (sigma-meson extension) \n"
    <<std::endl;
}
