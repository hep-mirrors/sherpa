#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Jet_Particle_Selector_Base: public Trigger_Base {  
  protected:

    ATOOLS::Flavour m_flavour;

    double m_xmin, m_xmax;
    size_t m_item;

  public:

    Jet_Particle_Selector_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);
    
    virtual bool Select(const Vec4D& p1,const Vec4D& p2) const = 0;

  };// end of class Jet_Particle_Selector_Base

  class Jet_Particle_DPhi_Selector: public Jet_Particle_Selector_Base {  
  public:

    Jet_Particle_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
			       const double min,const double max,
			       const std::string &inlist,const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Jet_Particle_DPhi_Selector

  class Jet_Particle_DEta_Selector: public Jet_Particle_Selector_Base {  
  public:

    Jet_Particle_DEta_Selector(const ATOOLS::Flavour flav,const size_t item,
			       const double min,const double max,
			       const std::string &inlist,
			       const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Jet_Particle_DEta_Selector

  class Jet_Particle_DY_Selector: public Jet_Particle_Selector_Base {  
  public:

    Jet_Particle_DY_Selector(const ATOOLS::Flavour flav,const size_t item,
			     const double min,const double max,
			     const std::string &inlist,
			     const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Jet_Particle_DY_Selector

  class Jet_Particle_DR_Selector: public Jet_Particle_Selector_Base {  
  public:

    Jet_Particle_DR_Selector(const ATOOLS::Flavour flav,const size_t item,
			     const double min,const double max,
			     const std::string &inlist,
			     const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Jet_Particle_DR_Selector

  class Jet_Particle_DRY_Selector: public Jet_Particle_Selector_Base {  
  public:

    Jet_Particle_DRY_Selector(const ATOOLS::Flavour flav,const size_t item,
			      const double min,const double max,
			      const std::string &inlist,
			      const std::string &outlist);
    
    bool Select(const Vec4D& p1,const Vec4D& p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Jet_Particle_DRY_Selector

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetJetParticleDeltaSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("LeadJets").Get<std::string>();
  const auto item = s["Item"].SetDefault(0).Get<size_t>();
  const auto kf = s["Flav"].SetDefault(kf_jet).Get<int>();
  ATOOLS::Flavour flav{ ATOOLS::Flavour((kf_code)(std::abs(kf))) };
  if (kf<0) flav=flav.Bar();
  return new Class(flav,item,min,max,inlist,outlist);
}

#define DEFINE_SELECTOR_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter				\
  <Analysis_Object,Analysis_Key,CLASS>::			\
  operator()(const Analysis_Key& key) const		\
  { return GetJetParticleDeltaSelector<CLASS>(key); }

#define DEFINE_SELECTOR_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Analysis_Key,CLASS>::	\
  PrintInfo(std::ostream &str,const size_t width) const		\
  { str<<"e.g. {Min: 30, Max: 70, InList: Jets, OutList: LeadJets, Item: 0, Flav: 93}"; }

#define DEFINE_JET_SELECTOR_DELTA_GETTER(CLASS,TAG)		\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_SELECTOR_GETTER_METHOD(CLASS)				\
  DEFINE_SELECTOR_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Jet_Particle_Selector_Base::
Jet_Particle_Selector_Base(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &outlist):
  Trigger_Base(inlist,outlist),
  m_flavour(flav), 
  m_item(item) 
{
  m_xmin=min;
  m_xmax=max;
}

void Jet_Particle_Selector_Base::Evaluate
(const ATOOLS::Particle_List &inlist,
 ATOOLS::Particle_List &outlist,double value,double ncount)
{
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
  for (size_t i=0;i<inlist.size();++i) if (pos!=i) {
    if (inlist[i]->Flav().Kfcode()==kf_jet) {
      if (!Select(inlist[pos]->Momentum(),inlist[i]->Momentum())) return;
    }
  }

  outlist.resize(inlist.size());
  for (size_t i=0;i<inlist.size();++i) 
    outlist[i] = new ATOOLS::Particle(*inlist[i]);
}

DEFINE_JET_SELECTOR_DELTA_GETTER(Jet_Particle_DPhi_Selector,"JetDPhiSel")

Jet_Particle_DPhi_Selector::
Jet_Particle_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
			   const double min,const double max,
			   const std::string &inlist,const std::string &outlist):
  Jet_Particle_Selector_Base(flav,item,min,max,inlist,outlist) {}

bool Jet_Particle_DPhi_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dphi(p1.DPhi(p2)/M_PI*180.0);
  if (dphi<m_xmin || dphi>m_xmax) return false;
  return true;
}

Analysis_Object *Jet_Particle_DPhi_Selector::GetCopy() const
{
  return new Jet_Particle_DPhi_Selector(m_flavour,m_item,
					m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_JET_SELECTOR_DELTA_GETTER(Jet_Particle_DEta_Selector,"JetDEtaSel")

Jet_Particle_DEta_Selector::
Jet_Particle_DEta_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &outlist):
  Jet_Particle_Selector_Base(flav,item,min,max,inlist,outlist) {}

bool Jet_Particle_DEta_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double deta=dabs(p1.Eta()-p2.Eta());
  if (deta<m_xmin || deta>m_xmax) return false;
  return true;
}

Analysis_Object *Jet_Particle_DEta_Selector::GetCopy() const
{
  return new Jet_Particle_DEta_Selector(m_flavour,m_item,
					m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_JET_SELECTOR_DELTA_GETTER(Jet_Particle_DY_Selector,"JetDYSel")

Jet_Particle_DY_Selector::
Jet_Particle_DY_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &outlist):
  Jet_Particle_Selector_Base(flav,item,min,max,inlist,outlist) {}

bool Jet_Particle_DY_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dy=dabs(p1.Y()-p2.Y());
  if (dy<m_xmin || dy>m_xmax) return false;
  return true;
}

Analysis_Object *Jet_Particle_DY_Selector::GetCopy() const
{
  return new Jet_Particle_DY_Selector(m_flavour,m_item,
				      m_xmin,m_xmax,m_inlist,m_outlist);
}


DEFINE_JET_SELECTOR_DELTA_GETTER(Jet_Particle_DR_Selector,"JetDRSel")

Jet_Particle_DR_Selector::
Jet_Particle_DR_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &outlist):
  Jet_Particle_Selector_Base(flav,item,min,max,inlist,outlist) {}

bool Jet_Particle_DR_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dr=sqrt(sqr(p1.Eta()-p2.Eta())+
		 sqr(p1.DPhi(p2)));
  if (dr<m_xmin || dr>m_xmax) return false;
  return true;
}

Analysis_Object *Jet_Particle_DR_Selector::GetCopy() const
{
  return new Jet_Particle_DR_Selector(m_flavour,m_item,
				      m_xmin,m_xmax,m_inlist,m_outlist);
}



DEFINE_JET_SELECTOR_DELTA_GETTER(Jet_Particle_DRY_Selector,"JetDRYSel")

Jet_Particle_DRY_Selector::
Jet_Particle_DRY_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const double min,const double max,
		  const std::string &inlist,const std::string &outlist):
  Jet_Particle_Selector_Base(flav,item,min,max,inlist,outlist) {}

bool Jet_Particle_DRY_Selector::Select(const Vec4D &p1,const Vec4D &p2) const
{
  double dr=sqrt(sqr(p1.Y()-p2.Y())+
		 sqr(p1.DPhi(p2)));
  if (dr<m_xmin || dr>m_xmax) return false;
  return true;
}

Analysis_Object *Jet_Particle_DRY_Selector::GetCopy() const
{
  return new Jet_Particle_DRY_Selector(m_flavour,m_item,
				       m_xmin,m_xmax,m_inlist,m_outlist);
}


