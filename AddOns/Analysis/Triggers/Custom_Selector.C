#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

namespace ANALYSIS {

  class Custom_Selector_Base: public Trigger_Base {  
  protected:

    double m_xmin, m_xmax;

  public:

    Custom_Selector_Base
    (const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);

    virtual bool Select(const ATOOLS::Particle_List &inlist) const = 0;
    
  };// end of class Custom_Selector_Base

  class SHT_Selector: public Custom_Selector_Base {  
  public:

    SHT_Selector(const double min,const double max,
		 const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle_List &inlist) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class SHT_Selector

  // hHTF: a selector for hadronic HT by Fernando:P:) OJO: based on PT!!!
  // based on SHTSel!
  class hHTF_Selector: public Custom_Selector_Base {  
  public:

    hHTF_Selector(const double min,const double max,
		 const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle_List &inlist) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class SHF_Selector

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetCustomSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("LeadJets").Get<std::string>();
  return new Class(min,max,inlist,outlist);
}

#define DEFINE_CUSTOM_GETTER_METHOD(CLASS)		\
  Analysis_Object *ATOOLS::Getter			\
  <Analysis_Object,Analysis_Key,CLASS>::		\
  operator()(const Analysis_Key& key) const	\
  { return GetCustomSelector<CLASS>(key); }

#define DEFINE_CUSTOM_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter						\
  <Analysis_Object,Analysis_Key,CLASS>::		\
  PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 30, Max: 70, InList: Jets, OutList: LeadJets}"; }

#define DEFINE_CUSTOM_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_CUSTOM_GETTER_METHOD(CLASS)			\
  DEFINE_CUSTOM_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Custom_Selector_Base::
Custom_Selector_Base(const double min,const double max,
		     const std::string &inlist,const std::string &outlist):
  Trigger_Base(inlist,outlist)
{
  m_xmin=min;
  m_xmax=max;
}

void Custom_Selector_Base::Evaluate(const ATOOLS::Particle_List &inlist,
					  ATOOLS::Particle_List &outlist,
					  double weight,double ncount)
{
  if (!Select(inlist)) return;
  outlist.resize(inlist.size());
  for (size_t i=0;i<inlist.size();++i) 
    outlist[i] = new ATOOLS::Particle(*inlist[i]);
}

DEFINE_CUSTOM_GETTER(SHT_Selector,"SHTSel")

SHT_Selector::SHT_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Custom_Selector_Base(min,max,inlist,outlist) {}

bool SHT_Selector::Select(const ATOOLS::Particle_List &inlist) const
{
  size_t item=0;
  double ht=0.;
  ATOOLS::Vec4D mom(0.,0.,0.,0.);
  for (size_t i=0;i<inlist.size();++i) {
    if (inlist[i]->Flav()==ATOOLS::Flavour(kf_jet)) {
      if (item>0) ht+=inlist[i]->Momentum().PPerp();
      item++;
    }
    else mom+=inlist[i]->Momentum();
  }

  ht+=mom.PPerp();

  return ht>=m_xmin && ht<=m_xmax;
}

Analysis_Object *SHT_Selector::GetCopy() const
{
  return new SHT_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_CUSTOM_GETTER(hHTF_Selector,"hHTFSel")

hHTF_Selector::hHTF_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Custom_Selector_Base(min,max,inlist,outlist) {}

bool hHTF_Selector::Select(const ATOOLS::Particle_List &inlist) const
{
//  size_t item=0;
  double ht=0.;
  for (size_t i=0;i<inlist.size();++i) {
    if (inlist[i]->Flav()==ATOOLS::Flavour(kf_jet)) {
      ht+=inlist[i]->Momentum().PPerp();
    }
//std::cout<<"i: "<<i<<" ht: "<<ht<<std::endl;
//      if (item>0) ht+=inlist[i]->Momentum().PPerp();
//      item++;
  }

  return ht>=m_xmin && ht<=m_xmax;
}

Analysis_Object *hHTF_Selector::GetCopy() const
{
  return new hHTF_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

