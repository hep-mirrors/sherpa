#include "METOOLS/HadronCurrents/VA_0_KKPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// Comment to be deleted once we have it all checked ....
// This is the KKpi sibling of VA_0_PiPiPi, added while harvesting
// - https://arxiv.org/pdf/1509.09140 (Sec.3.3)
// Structurally identical to VA_0_PiPiPi::Calc() except that the F5
// (Levi-Civita / vector) term is not assumed to vanish here.
// The actual CPC/RChL form factors (Sec.3.3) are NOT implemented yet -
// see FF_0_PPP.C, "F1_0_KKPi" et al. This file exists to fix the
// current's Lorentz structure and the KKpi ordering convention so we
// have something concrete to iterate on.
//
///////////////////////////////////////////////////////////////////////////

VA_0_KKPi::VA_0_KKPi(const ATOOLS::Flavour_Vector& flavs,
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

VA_0_KKPi::~VA_0_KKPi() {
  if (p_f1) { delete p_f1; p_f1 = NULL; }
  if (p_f2) { delete p_f2; p_f2 = NULL; }
  if (p_f3) { delete p_f3; p_f3 = NULL; }
  if (p_fS) { delete p_fS; p_fS = NULL; }
}


void VA_0_KKPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  // Convention: 0 = K^+, 1 = K^-, 2 = pi^- (see header).
  Vec4D p1    = moms[p_i[0]],  p2 = moms[p_i[1]],  p3 = moms[p_i[2]];
  Vec4D q     = p1+p2+p3,    dq21 = p2-p1,       dq31 = p3-p1;
  double s123 = q.Abs2(),      s1 = (q-p1).Abs2(), s2 = (q-p2).Abs2();
  double Qq21 = q*dq21/s123, Qq31 = q*dq31/s123;
  Vec4D  v1   = dq21-Qq21*q,   v2 = dq31-Qq31*q;
  Vec4C  v4   = Vec4C(cross(p1,p2,p3));
  // per-event debug output (re-enabled - this is the K+K-pi- RChL2012
  // channel actively under investigation for the Q^2/K+K- peak-
  // position question). sKK=(p1+p2)^2=(K+K-)^2 is printed explicitly
  // since that is the specific sub-mass in question (propRho's own
  // argument, per F1_0_KKPi::FF_RChL2012's paper_s2 mapping) - s1/s2
  // above are this Current's own local (K-pi- same-sign)/(K+pi-)
  // pairs, not the same quantity.
  double sKK  = (p1+p2).Abs2();
  /* per-event debug output, commented out again (K1(1270) width
     investigation concluded - see conversation history) - uncomment
     to re-enable if needed for future debugging.
  msg_Out()<<"*** "<<METHOD<<": "<<m_flavs[p_i[0]]<<" + "<<m_flavs[p_i[1]]<<" + "
	   <<m_flavs[p_i[2]]<<":  Q^2 = "<<s123<<" (Q = "<<sqrt(s123)<<" GeV),  "
	   <<"s(K+K-) = "<<sKK<<" (sqrt = "<<sqrt(sKK)<<" GeV),  "
	   <<"s(K+pi-) = "<<s2<<" (sqrt = "<<sqrt(s2)<<" GeV),  "
	   <<"s(K-pi-) = "<<s1<<" (sqrt = "<<sqrt(s1)<<" GeV)\n";
  */
  Complex F1  = (p_f1!=NULL ? (*p_f1)(moms) : Complex(0.,0.));
  Complex F2  = (p_f2!=NULL ? (*p_f2)(moms) : Complex(0.,0.));
  Complex F3  = (p_f3!=NULL ? (*p_f3)(moms) : Complex(0.,0.));
  Complex FS  = (p_fS!=NULL ? (*p_fS)(moms) : Complex(0.,0.));
  // msg_Out()<<"    F1 = "<<F1<<",  F2 = "<<F2<<",  F3 = "<<F3<<",  FS = "<<FS<<"\n";
  // NULL guards above: the FF_0_PPP getter dispatch returns NULL for
  // any flavour combination it doesn't recognise (e.g. a KKpi isospin
  // variant with no FixMode() branch yet) - dereferencing an
  // unconditional (*p_f1)(moms) etc. in that case is a NULL-pointer
  // crash, not the "falls back to zero" behaviour the rest of this
  // codebase assumes. Confirmed by an actual crash log for exactly
  // this reason (K^- K_S pi^0, an intentionally-unrecognised channel).
  Insert( m_norm * (F1*v1 + F2*v2 + F3*q + FS*v4), 0);
}

void VA_0_KKPi::SetModelParameters(struct GeneralModel model) {
  map<string,double> pmap;
  FF_Parameters params(ff_model(model[string("FORM_FACTOR")]),
		       m_flavs,p_i,pmap,"",&model);
  params.m_name = "F1_0_KKPi";
  p_f1 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F2_0_KKPi";
  p_f2 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F3_0_KKPi";
  p_f3 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "FS_0_KKPi";
  p_fS = FF_Getter::GetObject("FF_0_PPP",params);
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_KKPi,"VA_0_KKPi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_KKPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow K^+ K^- \\pi^- $ \n\n"
    <<"Order: 0 = $K^+$, 1 = $K^-$, 2 = $\\pi^-$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 201 :} RChL2012 (arXiv:1203.3955 "
    <<"Sec.2.2, Eqs.12-15; F3 is identically 0 for this channel) \n"
    <<"  \\end{itemize} \n"
    <<"Status: form factors for FORM_FACTOR=100,101,200,202 not implemented "
    <<"(return 0). \n\n"
    <<"Reference: https://arxiv.org/pdf/1203.3955 (Sec.2.2) \n"
    <<std::endl;
}
