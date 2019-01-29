#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Two_Particle_X_Selector_Base: public Two_List_Trigger_Base {  
  protected:

    ATOOLS::Flavour m_flavour;

    double m_xmin, m_xmax;
    int m_item;

  public:

    Two_Particle_X_Selector_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const double min,const double max,
     const std::string &inlist,const std::string &reflist,
     const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);
    
    virtual bool Select(const Vec4D& p1,const Vec4D& p2) const = 0;

  };// end of class Two_Particle_X_Selector_Base

  class Two_DPhiL_Selector: public Two_Particle_X_Selector_Base {  
  public:

    Two_DPhiL_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DPhi_Selector

  class Two_DEtaL_Selector: public Two_Particle_X_Selector_Base {  
  public:

    Two_DEtaL_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DEta_Selector

  class Two_DYL_Selector: public Two_Particle_X_Selector_Base {  
  public:

    Two_DYL_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DY_Selector

  class Two_PTL_Selector: public Two_Particle_X_Selector_Base {  
  public:

    Two_PTL_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_PT_Selector

  class Two_DRL_Selector: public Two_Particle_X_Selector_Base {  
  public:

    Two_DRL_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DR_Selector

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetTwoParticleLDeltaSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"].SetDefault("Jets").Get<std::string>();
  const auto reflist = s["RefList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("LeadJets").Get<std::string>();
  const auto item = s["Item"].SetDefault(0).Get<size_t>();
  const auto kf = s["Flav"].SetDefault(kf_jet).Get<int>();
  ATOOLS::Flavour flav{ ATOOLS::Flavour((kf_code)(std::abs(kf))) };
  if (kf<0) flav=flav.Bar();
  return new Class(flav,item,min,max,inlist,reflist,outlist);
}

#define DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter					\
  <Analysis_Object,Analysis_Key,CLASS>::				\
  operator()(const Analysis_Key& key) const			\
  { return GetTwoParticleLDeltaSelector<CLASS>(key); }

#define DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(CLASS)	\
  void ATOOLS::Getter					\
  <Analysis_Object,Analysis_Key,CLASS>::		\
  PrintInfo(std::ostream &str,const size_t width) const \
  { str<<"e.g. {Min: 30, Max: 70, InList: Jets, RefList: Jets, OutList: LeadJets, Item: 0, Flav: 93}"; }

#define DEFINE_TWO_SELECTOR_DELTA_GETTER(CLASS,TAG)		\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS)		\
  DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Two_Particle_X_Selector_Base::
Two_Particle_X_Selector_Base(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_List_Trigger_Base(inlist,reflist,outlist),
  m_flavour(flav), 
  m_item(item) 
{
  m_xmin=min;
  m_xmax=max;
}

void Two_Particle_X_Selector_Base::Evaluate
(const ATOOLS::Particle_List &inlist,const ATOOLS::Particle_List &reflist,
 ATOOLS::Particle_List &outlist,double value,double ncount)
{

  Vec4D refmom(0.,0.,0.,0.);
  for (size_t i=0;i<reflist.size();++i) {
    refmom+=reflist[i]->Momentum();
  }
  if (m_item>=0) {
    int no=-1; 
    size_t pos=std::string::npos;
    for (size_t i=0;i<inlist.size();++i) {
      if (inlist[i]->Flav()==m_flavour || 
	  m_flavour.Kfcode()==kf_none) {
	++no;
	if (no==(int)m_item) {
	  pos=i;
	  break;
	}
      }
    }
    if (pos==std::string::npos) return;
    if (!Select(inlist[pos]->Momentum(),refmom)) return;
  }
  else {
    for (size_t i=0;i<inlist.size();++i) {
      if (inlist[i]->Flav()==m_flavour) {
	if (!Select(inlist[i]->Momentum(),refmom)) return;
      }
    }
  }

  outlist.resize(inlist.size());
  for (size_t i=0;i<inlist.size();++i) 
    outlist[i] = new ATOOLS::Particle(*inlist[i]);

}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DPhiL_Selector,"TwoDPhiXSel")

Two_DPhiL_Selector::
Two_DPhiL_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_X_Selector_Base(flav,item,min,max,
			     inlist,reflist,outlist) {}

bool Two_DPhiL_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dphi(p1.DPhi(p2)/M_PI*180.0);
  if (dphi<m_xmin || dphi>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DPhiL_Selector::GetCopy() const
{
  return new Two_DPhiL_Selector(m_flavour,m_item,
			       m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DEtaL_Selector,"TwoDEtaXSel")

Two_DEtaL_Selector::
Two_DEtaL_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_X_Selector_Base(flav,item,min,max,
			     inlist,reflist,outlist) {}

bool Two_DEtaL_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double deta=dabs(p1.Eta()-p2.Eta());
  if (deta<m_xmin || deta>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DEtaL_Selector::GetCopy() const
{
  return new Two_DEtaL_Selector(m_flavour,m_item,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DYL_Selector,"TwoDYXSel")

Two_DYL_Selector::
Two_DYL_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_X_Selector_Base(flav,item,min,max,
			     inlist,reflist,outlist) {}

bool Two_DYL_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dy=dabs(p1.Y()-p2.Y());
  if (dy<m_xmin || dy>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DYL_Selector::GetCopy() const
{
  return new Two_DYL_Selector(m_flavour,m_item,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}


DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PTL_Selector,"TwoPTXSel")

Two_PTL_Selector::
Two_PTL_Selector(const ATOOLS::Flavour flav,const size_t item,
		const double min,const double max,
		const std::string &inlist,const std::string &reflist,
		const std::string &outlist):
  Two_Particle_X_Selector_Base(flav,item,min,max,
			     inlist,reflist,outlist) {}

bool Two_PTL_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double pt=(p1+p2).PPerp();
  if (pt<m_xmin || pt>m_xmax) return false;
  return true;
}

Analysis_Object *Two_PTL_Selector::GetCopy() const
{
  return new Two_PTL_Selector(m_flavour,m_item,
			     m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}


DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DRL_Selector,"TwoDRXSel")

Two_DRL_Selector::
Two_DRL_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_X_Selector_Base(flav,item,min,max,
			     inlist,reflist,outlist) {}

bool Two_DRL_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dr=sqrt(sqr(p1.Eta()-p2.Eta())+
		 sqr(p1.DPhi(p2)));
  if (dr<m_xmin || dr>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DRL_Selector::GetCopy() const
{
  return new Two_DRL_Selector(m_flavour,m_item,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}


