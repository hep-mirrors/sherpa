#include "METOOLS/HadronCurrents/VA_0_EtaPiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// eta pi^- pi^0 / eta' pi^- pi^0 - see VA_0_EtaPiPi.H for the design
// note. Same Calc() structure as VA_0_KPiK/VA_0_KPiPi (kept generic
// even though F1/F2/F3 are always 0 here, for consistency with the
// other 3-meson Current classes and in case a dedicated RChiPT current
// with nonzero F1/F2 is ever added - see the FIXME in FF_0_PPP.C's
// FS_0_EtaPiPi::Construct()).
//
///////////////////////////////////////////////////////////////////////////

VA_0_EtaPiPi::VA_0_EtaPiPi(const ATOOLS::Flavour_Vector& flavs,
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

VA_0_EtaPiPi::~VA_0_EtaPiPi() {
  if (p_f1) { delete p_f1; p_f1 = NULL; }
  if (p_f2) { delete p_f2; p_f2 = NULL; }
  if (p_f3) { delete p_f3; p_f3 = NULL; }
  if (p_fS) { delete p_fS; p_fS = NULL; }
}

void VA_0_EtaPiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p1    = moms[p_i[0]],  p2 = moms[p_i[1]],  p3 = moms[p_i[2]];
  Vec4D q     = p1+p2+p3,    dq13 = p1-p3,       dq23 = p2-p3;
  double s123 = q.Abs2();
  double Qq13 = q*dq13/s123, Qq23 = q*dq23/s123;
  Vec4D  v1   = dq13-Qq13*q,   v2 = dq23-Qq23*q;
  Vec4C  v4   = Vec4C(cross(p1,p2,p3));
  Complex F1  = (p_f1!=NULL ? (*p_f1)(moms) : Complex(0.,0.));
  Complex F2  = (p_f2!=NULL ? (*p_f2)(moms) : Complex(0.,0.));
  Complex F3  = (p_f3!=NULL ? (*p_f3)(moms) : Complex(0.,0.));
  Complex FS  = (p_fS!=NULL ? (*p_fS)(moms) : Complex(0.,0.));
  // NULL guards above: the FF_0_PPP getter dispatch returns NULL for
  // any flavour combination it doesn't recognise - dereferencing an
  // unconditional (*p_f1)(moms) etc. in that case is a NULL-pointer
  // crash, not a graceful "falls back to zero" (confirmed by an
  // actual crash log, VA_0_KKPi/K^- K_S pi^0 - same fix applied
  // uniformly across every VA_0_XXX 3-meson Current class).
  Insert( m_norm * (F1*v1 + F2*v2 + F3*q + FS*v4), 0);
}

void VA_0_EtaPiPi::SetModelParameters(struct GeneralModel model) {
  map<string,double> pmap;
  FF_Parameters params(ff_model(model[string("FORM_FACTOR")]),
		       m_flavs,p_i,pmap,"",&model);
  params.m_name = "F1_0_EtaPiPi";
  p_f1 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F2_0_EtaPiPi";
  p_f2 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "F3_0_EtaPiPi";
  p_f3 = FF_Getter::GetObject("FF_0_PPP",params);
  params.m_name = "FS_0_EtaPiPi";
  p_fS = FF_Getter::GetObject("FF_0_PPP",params);
}

DEFINE_CURRENT_GETTER(METOOLS::VA_0_EtaPiPi,"VA_0_EtaPiPi")

void ATOOLS::Getter<METOOLS::Current_Base,
		    METOOLS::ME_Parameters,METOOLS::VA_0_EtaPiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\eta \\pi \\pi $ (eta(') pi pi family) \n\n"
    <<"Order: 0 = eta (or eta'), 1 = pi^-, 2 = pi^0 - this exact order \n"
    <<"is required (not interchangeable): \n"
    <<"  eta pi^- pi^0,   eta' pi^- pi^0 \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 100 :} old TAUOLA/VMD current \n"
    <<"      (anomalous vector form factor only - G-parity forbids the \n"
    <<"      axial/scalar pieces in the isospin limit) \n"
    <<"  \\end{itemize} \n"
    <<"Status: overall WZW-anomaly/eta-eta' mixing normalisation \n"
    <<"(N_etapipi/N_etaprimepipi) is a placeholder - see FS_0_EtaPiPi's \n"
    <<"FIXME comment in FF_0_PPP.C. No dedicated RChiPT current \n"
    <<"(Gomez Dumm & Roig) is implemented - falls back to the constant \n"
    <<"form factor for FORM_FACTOR != 100. \n\n"
    <<"Reference: tau_two_meson_currents_KS_RChiT.tex, Sec. \n"
    <<"'eta pi pi and eta' pi pi' \n"
    <<std::endl;
}
