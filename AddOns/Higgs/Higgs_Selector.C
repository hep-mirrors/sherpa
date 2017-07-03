#include "PHASIC++/Selectors/Selector.H"

namespace HIGGS {

  class Higgs_Selector: public PHASIC::Selector_Base {
  private:

    double m_pt1, m_pt2, m_eta, m_mmin, m_mmax, m_dr, m_epspt;

    bool Trigger(ATOOLS::Vec4D &py1,ATOOLS::Vec4D &py2,
		 ATOOLS::Vec4D &pj);

  public:

    Higgs_Selector(PHASIC::Process_Base *const proc,
                   double pt1,double pt2,double eta,
		   double mmin,double mmax,
		   double dr,double epspt);

    ~Higgs_Selector();

    bool   Trigger(const ATOOLS::Vec4D_Vector &,
                   ATOOLS::NLO_subevt *const sub=NULL);

    void   BuildCuts(PHASIC::Cut_Data *cuts);

  };

}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"

using namespace HIGGS;
using namespace PHASIC;
using namespace ATOOLS;

Higgs_Selector::Higgs_Selector(Process_Base *const proc,
			       double pt1, double pt2, double eta,
			       double mmin, double mmax,
			       double dr,double epspt):
  Selector_Base("HiggsFinder"),
  m_pt1(pt1), m_pt2(pt2), m_eta(eta),
  m_mmin(mmin), m_mmax(mmax),
  m_dr(dr), m_epspt(epspt)
{
  m_smin       = sqr(m_mmin);
}

Higgs_Selector::~Higgs_Selector()
{
}

void Higgs_Selector::BuildCuts(PHASIC::Cut_Data *cuts)
{
  for (int i=m_nin;i<m_n;++i)
    if (p_fl[i].IsPhoton())
      for (int j=m_nin;j<m_n;++j)
        if (p_fl[j].IsPhoton()) {
	  cuts->scut[j][i] = 
	    cuts->scut[i][j] = Max(sqr(m_mmin),cuts->scut[i][j]);
	}
}

bool Higgs_Selector::Trigger(const Vec4D_Vector &p,
                             ATOOLS::NLO_subevt *const sub)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  size_t n(sub?sub->m_n:m_n);
  const Flavour *const fl(sub?sub->p_fl:p_fl);
  Vec4D py1, py2, pj;
  for (size_t i(m_nin);i<n;++i) {
    if (fl[i].IsPhoton()) {
      if (py1==Vec4D()) py1=p[i];
      else {
        if (py2!=Vec4D()) msg_Error()<<METHOD<<"(): Not a yy event."<<std::endl;
        py2=p[i];
      }
    }
    if (fl[i].Strong()) {
      if (pj!=Vec4D()) msg_Error()<<METHOD<<"(): Not a yy event."<<std::endl;
      pj=p[i];
    }
  }
  return Trigger(py1,py2,pj);
}

bool Higgs_Selector::Trigger(Vec4D &py1,Vec4D &py2,Vec4D &pj)
{
  if (py1==Vec4D() || py2==Vec4D())
    msg_Error()<<METHOD<<"(): Not a yy event."<<std::endl;
  if (py2.PPerp2()>py1.PPerp2()) std::swap<Vec4D>(py1,py2);
  bool trigger(true);
  if (py1.PPerp2()<sqr(m_pt1)) trigger=false;
  if (py2.PPerp2()<sqr(m_pt2)) trigger=false;
  if (dabs(py1.Eta())>m_eta) trigger=false;
  if (dabs(py2.Eta())>m_eta) trigger=false;
  double m2((py1+py2).Abs2());
  if (m2<sqr(m_mmin) || m2>sqr(m_mmax)) trigger=false;
  if (pj.PPerp()>m_epspt) {
    if (py1.DR(pj)<m_dr) trigger=false;
    if (py2.DR(pj)<m_dr) trigger=false;
  }
  return (1-m_sel_log->Hit(1-trigger));
}


DECLARE_ND_GETTER(Higgs_Selector,"HiggsFinder",Selector_Base,Selector_Key,true);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,Higgs_Selector>::
operator()(const Selector_Key &key) const
{
  if (key.empty() || key.front().size()<5) THROW(critical_error,"Invalid syntax");
 
  double dr=0.0, epspt=1.0e12;
  if (key.front().size()>6) {
    dr=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][5]));
    epspt=ToType<double>(key.p_read->Interpreter()->Interprete(key[0][6]));
  }

  Higgs_Selector *jf
    (new Higgs_Selector
     (key.p_proc,
      ToType<double>(key.p_read->Interpreter()->Interprete(key[0][0])),
      ToType<double>(key.p_read->Interpreter()->Interprete(key[0][1])),
      ToType<double>(key.p_read->Interpreter()->Interprete(key[0][2])),
      ToType<double>(key.p_read->Interpreter()->Interprete(key[0][3])),
      ToType<double>(key.p_read->Interpreter()->Interprete(key[0][4])),dr,epspt));
  jf->SetProcess(key.p_proc);
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Higgs_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"HiggsFinder pt1 pt2 eta mmin mmax [dr epspt]";
}
