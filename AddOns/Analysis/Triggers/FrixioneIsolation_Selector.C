#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>

using namespace ATOOLS;

namespace ANALYSIS {

  struct edr {
    double E, dr;
    edr(double _E,double _dr) : E(_E), dr(_dr) {}
  };
  class Order_edr {
  public:
    int operator()(const edr a, const edr b);
  };
  int Order_edr::operator()(const edr a, const edr b) {
    if (a.dr<b.dr) return 1;
    return 0;
  }

  class FrixioneIsolation_Selector: public Two_List_Trigger_Base {  
  protected:

    double m_d0, m_n, m_eps;

    double Chi(double eg,double dr)
    {
      if (m_n==0) return m_eps;
      if (m_n<0) return 0.;
      return eg*pow((1.-cos(dr))/(1.-cos(m_d0)),m_n);
    }

    double DR(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2)
    {
      return  sqrt(sqr(DEta12(p1,p2)) + sqr(DPhi12(p1,p2)));
    }

    double DEta12(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2)
    {
      // eta1,2 = -log(tan(theta_1,2)/2)   
      double c1=p1[3]/Vec3D(p1).Abs();
      double c2=p2[3]/Vec3D(p2).Abs();
      return  0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
    }

    double DPhi12(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2)
    {
      double pt1=sqrt(p1[1]*p1[1]+p1[2]*p1[2]);
      double pt2=sqrt(p2[1]*p2[1]+p2[2]*p2[2]);
      return acos(Min(1.0,Max(-1.0,((p1[1]*p2[1]+p1[2]*p2[2])/(pt1*pt2)))));
    }

  public:

    inline FrixioneIsolation_Selector
    (const double d0,const double n,const double eps,
     const std::string &inlist,const std::string &reflist,
     const std::string &outlist):
      Two_List_Trigger_Base(inlist,reflist,outlist),
      m_d0(d0), m_n(n), m_eps(eps) {}
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,double value,double ncount)
    {
      for (size_t j=0;j<reflist.size();++j) {
	double egamma=reflist[j]->Momentum().PPerp();
	std::vector<edr> edrlist;
	for (size_t i=0;i<inlist.size();++i) {
	  double dr=DR(reflist[j]->Momentum(),inlist[i]->Momentum());
	  if (dr<m_d0) edrlist.push_back(edr(inlist[i]->Momentum().PPerp(),dr));
	}
	if (edrlist.size()>0) {
	  stable_sort(edrlist.begin(),edrlist.end(),Order_edr());
	  double etot=0.;
	  for (size_t i=0;i<edrlist.size();i++) {
	    etot+=edrlist[i].E;
	    if (etot<Chi(egamma,edrlist[i].dr)) return;
	  }
	  edrlist.clear();
	}
      }
      outlist.resize(inlist.size());
      for (size_t i=0;i<inlist.size();++i) 
	outlist[i] = new ATOOLS::Particle(*inlist[i]);
    }

    Analysis_Object *GetCopy() const 
    {
      return new FrixioneIsolation_Selector
	(m_d0,m_n,m_eps,m_inlist,m_reflist,m_outlist);
    }

  };// end of class FrixioneIsolation_Selector

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(FrixioneIsolation_Selector,"IsolationCutSel",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Argument_Matrix,FrixioneIsolation_Selector>::
operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()<1 || parameters[0].size()<6) return NULL;
  int eps(parameters[0].size()>6?ATOOLS::ToType<double>(parameters[0][6]):1.0);
  int crit1=ToType<int>(parameters[0][0]);
  Flavour flav=Flavour((kf_code)abs(crit1));
  if (crit1<0) flav=flav.Bar();
  return new FrixioneIsolation_Selector
    (ATOOLS::ToType<double>(parameters[0][1]),
     ATOOLS::ToType<double>(parameters[0][2]),eps,
     parameters[0][3],parameters[0][4],parameters[0][5]);
}									

void ATOOLS::Getter
<Analysis_Object,Argument_Matrix,FrixioneIsolation_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"kf delta_0 n inlist reflist outlist [epsilon]";
}
