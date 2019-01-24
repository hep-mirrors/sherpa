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

    bool   Trigger(ATOOLS::Selector_List &);

    void   BuildCuts(PHASIC::Cut_Data *cuts);

  };

}

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

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

bool Higgs_Selector::Trigger(Selector_List &sl)
{
  DEBUG_FUNC(m_on);
  if (!m_on) return true;
  Vec4D py1, py2, pj;
  for (size_t i(m_nin);i<sl.size();++i) {
    if (sl[i].Flavour().IsPhoton()) {
      if (py1==Vec4D()) py1=sl[i].Momentum();
      else {
        if (py2!=Vec4D()) msg_Error()<<METHOD<<"(): Not a yy event."<<std::endl;
        py2=sl[i].Momentum();
      }
    }
    if (sl[i].Flavour().Strong()) {
      if (pj!=Vec4D()) msg_Error()<<METHOD<<"(): Not a yy event."<<std::endl;
      pj=sl[i].Momentum();
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
  auto s = key.m_settings;
  const auto pt1 = s["PT1"].SetDefault(0.0).Get<double>();
  const auto pt2 = s["PT2"].SetDefault(0.0).Get<double>();
  const auto eta = s["Eta"].SetDefault(0.0).Get<double>();
  const auto mrange = s["MassRange"]
    .SetDefault({0.0, std::numeric_limits<double>::max()})
    .GetVector<double>();
  if (mrange.size() != 2)
    THROW(fatal_error, "Invalid syntax, MassRange needs two entries");
  const auto dr = s["DR"].SetDefault(0.0).Get<double>();
  const auto epspt = s["EpsPT"].SetDefault(1.0e12).Get<double>();
  Higgs_Selector *jf(new Higgs_Selector(key.p_proc,
                                        pt1, pt2,
                                        eta,
                                        mrange[0], mrange[1],
                                        dr, epspt));
  jf->SetProcess(key.p_proc);
  return jf;
}

void ATOOLS::Getter<Selector_Base,Selector_Key,Higgs_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"{Type: HiggsFinder, PT1: pt1, PT2: pt2, Eta: eta, MassRange: [mmin, mmax], DR: dR, EpsPT: pt}";
}
