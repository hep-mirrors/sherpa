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
GetCustomSelector(const Argument_Matrix &parameters) 
{				
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    return new Class(ATOOLS::ToType<int>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     parameters[0][2],parameters[0][3]);
  }
  if (parameters.size()<4) return NULL;
  double min=30.0, max=70.0; 
  std::string inlist="Jets", outlist="LeadJets";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
  }
  return new Class(min,max,inlist,outlist);
}									

#define DEFINE_CUSTOM_GETTER_METHOD(CLASS,NAME)	\
  Analysis_Object *				\
  NAME::operator()(const Argument_Matrix &parameters) const	\
  { return GetCustomSelector<CLASS>(parameters); }

#define DEFINE_CUSTOM_PRINT_METHOD(NAME)			\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max inlist outlist"; }

#define DEFINE_CUSTOM_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_CUSTOM_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_CUSTOM_PRINT_METHOD(NAME)

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

DEFINE_CUSTOM_GETTER(SHT_Selector,
		     SHT_Selector_Getter,"SHTSel")

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

DEFINE_CUSTOM_GETTER(hHTF_Selector,
		     hHTF_Selector_Getter,"hHTFSel")

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

