#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "ATOOLS/Org/Message.H"
namespace ANALYSIS {

  class HT_Selector: public Trigger_Base {  
  protected:

    double m_xmin, m_xmax;
    int    m_mode;

  public:

    inline HT_Selector(const double min,const double max,const int mode,
		       const std::string &inlist,const std::string &outlist):
      Trigger_Base(inlist,outlist), m_xmin(min), m_xmax(max), m_mode(mode) {}
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,double value,double ncount)
    {
      double HT=0.0;
      for (ATOOLS::Particle_List::const_iterator
	     pit(inlist.begin());pit!=inlist.end();++pit)
	HT+=m_mode?(*pit)->Momentum().EPerp():(*pit)->Momentum().PPerp();
      if (HT>=m_xmin && HT<=m_xmax) {
	outlist.resize(inlist.size());
	for (size_t i=0;i<inlist.size();++i) 
	  outlist[i] = new ATOOLS::Particle(*inlist[i]);
      }
    }

    Analysis_Object *GetCopy() const 
    {
      return new HT_Selector(m_xmin,m_xmax,m_mode,m_inlist,m_outlist);
    }

  };// end of class HT_Selector

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

DECLARE_GETTER(HT_Selector,"HTSel",Analysis_Object,Analysis_Key);

Analysis_Object *ATOOLS::Getter
<Analysis_Object,Analysis_Key,HT_Selector>::operator()
(const Analysis_Key& key) const
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  if (parameters.size() < 4)
    THROW(missing_input, "Missing parameter values.");
  const int mode{
    parameters.size() > 4 ? s.Interprete<int>(parameters[4]) : 0 };
  return new HT_Selector(s.Interprete<double>(parameters[0]),
                         s.Interprete<double>(parameters[1]),mode,
                         parameters[2],parameters[3]);
}

void ATOOLS::Getter
<Analysis_Object,Analysis_Key,HT_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"[min, max, inlist, outlist, mode]  ... mode is optional and defaults to 0";
}
