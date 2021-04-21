#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Math/Variable.H"
#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {

  class One_Particle_Extractor: public Two_List_Trigger_Base {  
  private:

    ATOOLS::Variable_Base<double> *p_variable;

    ATOOLS::Flavour m_flavour;

    double m_xmin, m_xmax;
    size_t m_item;

  public:

    One_Particle_Extractor(const std::string &type,const ATOOLS::Flavour flav,
			   const size_t item,const double min,const double max,
			   const std::string &inlist,const std::string &reflist,
			   const std::string &outlist);
    
    ~One_Particle_Extractor();

    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,
		  double weight=1.,double ncount=1);
    
    Analysis_Object *GetCopy() const;
    
  };// end of class One_Particle_Extractor

}// end of namespace ANALYSIS

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetOneParticleSelector(const Analysis_Key& key)
{									
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto type = s["Type"].SetDefault("p_\\perp").Get<double>();
  const auto inlist = s["InList"].SetDefault("Jets").Get<std::string>();
  const auto reflist = s["RefList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("LeadJets").Get<std::string>();
  const auto item = s["Item"].SetDefault(0).Get<size_t>();
  const auto kf = s["Flav"].SetDefault(kf_jet).Get<int>();
  const auto flav = Flavour{kf};
  return new Class(type,flav,item,min,max,inlist,reflist,outlist);
}									

#define DEFINE_ONE_EXTRACTOR_GETTER_METHOD(CLASS,NAME)	\
  Analysis_Object *				\
  NAME::operator()(const Analysis_Key& key) const	\
  { return GetOneParticleSelector<CLASS>(parameters); }

#define DEFINE_ONE_EXTRACTOR_PRINT_METHOD(NAME)			\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"e.g. {Min: 30, Max: 70, Type: "p_\\\\perp", InList: Jets, RefList: Jets, OutList: LeadJets, Item: 0, Flav: 93}"; }

#define DEFINE_ONE_EXTRACTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_ONE_EXTRACTOR_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_ONE_EXTRACTOR_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

DEFINE_ONE_EXTRACTOR_GETTER(One_Particle_Extractor,
			    One_Particle_Extractor_Getter,"OnePartExt")

One_Particle_Extractor::
One_Particle_Extractor(const std::string &type,const ATOOLS::Flavour flav,
		       const size_t item,const double min,const double max,
		       const std::string &inlist,const std::string &reflist,
		       const std::string &outlist):
  Two_List_Trigger_Base(inlist,reflist,outlist),
  p_variable(ATOOLS::Variable_Getter::GetObject(type,"")),
  m_flavour(flav), m_item(item)
{
  m_xmin=min;
  m_xmax=max;
}

One_Particle_Extractor::~One_Particle_Extractor()
{
  delete p_variable;
}

void One_Particle_Extractor::Evaluate(const ATOOLS::Particle_List &inlist,
				      const ATOOLS::Particle_List &reflist,
				      ATOOLS::Particle_List &outlist,
				      double weight,double ncount)
{
  int no=-1; 
  for (size_t i=0;i<reflist.size();++i) {
    if (reflist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==kf_none) {
      if (++no==(int)m_item) {
	double pt=(*p_variable)(&reflist[i]->Momentum());
	for (size_t j=0;j<inlist.size();++j) {
	  if (!(*inlist[j]==*reflist[i]) || pt<m_xmin || pt>m_xmax) 
	    outlist.push_back(new ATOOLS::Particle(*inlist[j]));
	}
	break;
      }
    }
  }
}

Analysis_Object *One_Particle_Extractor::GetCopy() const
{
  return new One_Particle_Extractor
    (p_variable->Name(),m_flavour,m_item,
     m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

