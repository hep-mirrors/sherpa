#include "AddOns/Analysis/Analyses/Analysis_Base.H"

namespace ANALYSIS {

  class xQ2_Analysis: public Analysis_Base {  
  public:

    xQ2_Analysis(const std::string &listname);

    void Evaluate(double weight,double ncount,int mode);
    Primitive_Observable_Base *Copy() const;

  };// end of class xQ2_Analysis

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

xQ2_Analysis::xQ2_Analysis
(const std::string &listname): Analysis_Base(listname)
{
  m_name="xQ2__"+m_listname;
  m_histos.resize(20,NULL);
  for (int i(0);i<4;++i) {
    m_histos[5*i+0] = new Histogram(11,1.0,1.0e5,100,"sigr");
    m_histos[5*i+1] = new Histogram(11,1.0,1.0e5,100,"sigr");
    m_histos[5*i+2] = new Histogram(11,1.0,1.0e5,100,"sigr");
    m_histos[5*i+3] = new Histogram(11,1.0,1.0e5,100,"sigr");
    m_histos[5*i+4] = new Histogram(11,1.0,1.0e5,100,"sigr");
  }
}

void xQ2_Analysis::Evaluate(double weight, double ncount,int mode)
{
  // compute observables
  const Blob_List *p_bl(p_ana->GetBlobList());
  Vec4D l(rpa->gen.PBeam(0)), lp, pp(rpa->gen.PBeam(1));
  Blob *me(p_bl->FindFirst(btp::Signal_Process));
  for (int i(0);i<me->NOutP();++i)
    if (me->OutParticle(i)->Flav().IsLepton()) {
      lp=me->OutParticle(i)->Momentum();
      break;
    }
  Vec4D qq(l-lp);
  double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), y((pp*qq)/(pp*l));
  double alpha(1./137.), w(Q2*Q2*x/(2.*M_PI*sqr(alpha)*(1.+sqr(1.-y)))*weight);
  w/=Q2*log(10.0)*rpa->Picobarn(); // Q2*log(10) arises from loarithmic plotting
  for (int i(0);i<4;++i) {
    int o(5*i);
    double s(pow(10.,i));
    FillHisto(o+0,Q2,(x>4.e-5*s && x<6.e-5*s)?w/(2.e-5*s):0.0,ncount,mode);
    FillHisto(o+1,Q2,(x>6.e-5*s && x<1.e-4*s)?w/(4.e-5*s):0.0,ncount,mode);
    FillHisto(o+2,Q2,(x>1.e-4*s && x<1.6e-4*s)?w/(6.e-5*s):0.0,ncount,mode);
    FillHisto(o+3,Q2,(x>1.6e-4*s && x<2.6e-4*s)?w/(1.e-4*s):0.0,ncount,mode);
    FillHisto(o+4,Q2,(x>2.6e-4*s && x<4.e-4*s)?w/(1.4e-4*s):0.0,ncount,mode);
  }
}

Primitive_Observable_Base *xQ2_Analysis::Copy() const 
{
  return new xQ2_Analysis(m_listname);
}

DECLARE_ND_GETTER(xQ2_Analysis,"xQ2",
		  Primitive_Observable_Base,Analysis_Key,true);

Primitive_Observable_Base *ATOOLS::Getter
<Primitive_Observable_Base,Analysis_Key,xQ2_Analysis>::
operator()(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  s.SetDefault("");
  const auto listname = s.Get<std::string>();
  if (listname.empty())
    THROW(missing_input, "Missing list.");
  return new xQ2_Analysis(listname);
}

void ATOOLS::Getter
<Primitive_Observable_Base,Analysis_Key,xQ2_Analysis>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list";
}

