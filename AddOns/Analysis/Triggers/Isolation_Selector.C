#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "ATOOLS/Org/Message.H"

#include <algorithm>

using namespace ATOOLS;

namespace ANALYSIS {

  class Isolation_Selector: public Two_List_Trigger_Base {  
  protected:

    double m_dr, m_emax;

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

    double DR(const ATOOLS::Vec4D & p1,const ATOOLS::Vec4D & p2)
    {
      return  sqrt(sqr(DEta12(p1,p2)) + sqr(DPhi12(p1,p2)));
    }

  public:

    inline Isolation_Selector
    (const double dr,const double emax,
     const std::string &inlist,const std::string &reflist,
     const std::string &outlist):
      Two_List_Trigger_Base(inlist,reflist,outlist),
      m_dr(dr), m_emax(emax) {}
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,double value,double ncount)
    {
      for (size_t j=0;j<reflist.size();++j) {
	double et=0.0;
	for (size_t i=0;i<inlist.size();++i) {
	  double dr=DR(reflist[j]->Momentum(),inlist[i]->Momentum());
	  if (dr<m_dr) et+=inlist[i]->Momentum().EPerp();
	}
	if (et>m_emax) return;
      }
      outlist.resize(inlist.size());
      for (size_t i=0;i<inlist.size();++i) 
	outlist[i] = new ATOOLS::Particle(*inlist[i]);
    }

    Analysis_Object *GetCopy() const 
    {
      return new Isolation_Selector
	(m_dr,m_emax,m_inlist,m_reflist,m_outlist);
    }

  };// end of class Isolation_Selector

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Isolation_Selector,"PhotonIsolation",
	       Analysis_Object,Analysis_Key);

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Analysis_Key,Isolation_Selector>::
operator()(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 5)
    THROW(missing_input, "Missing parameter values.");
  return new Isolation_Selector
    (s.Interprete<double>(parameters[0]),
     s.Interprete<double>(parameters[1]),
     parameters[2],parameters[3],parameters[4]);
}									

void ATOOLS::Getter
<Analysis_Object,Analysis_Key,Isolation_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[DR, E_max, inlist, reflist, outlist]";
}
