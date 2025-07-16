#include "AddOns/Analysis/Triggers/Trigger_Base.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"

#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {

  class Particle_Selector_Base: public Trigger_Base {  
  protected:

    double m_xmin, m_xmax;

  public:

    Particle_Selector_Base(const double min,const double max,
			   const std::string &inlist,
			   const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);
    
    virtual bool Select(const ATOOLS::Particle *p) const = 0;

  };// end of class Particle_Selector_Base

  class PT_Selector: public Particle_Selector_Base {  
  public:

    PT_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class PT_Selector

  class ET_Selector: public Particle_Selector_Base {  
  public:

    ET_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class ET_Selector

  class Eta_Selector: public Particle_Selector_Base {  
  public:

    Eta_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Eta_Selector

  class Abs_Eta_Selector: public Particle_Selector_Base {  
  public:

    Abs_Eta_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Abs_Eta_Selector

  class Y_Selector: public Particle_Selector_Base {  
  public:

    Y_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Y_Selector

  class Abs_Y_Selector: public Particle_Selector_Base {  
  public:

    Abs_Y_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Abs_Y_Selector

  class Phi_Selector: public Particle_Selector_Base {  
  public:

    Phi_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class Phi_Selector

  class DPhi_Selector: public Two_List_Trigger_Base {  
  private:

    size_t m_item;
    double m_xmin, m_xmax;

    ATOOLS::Flavour m_flavour;

  public:

    DPhi_Selector(const double min,const double max,
		  const ATOOLS::Flavour flav,
		  const size_t item,const std::string &reflist,
		  const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);
    
    Analysis_Object *GetCopy() const;
    
  };// end of class DPhi_Selector

}// end of namespace ANALYSIS

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ANALYSIS;
using namespace ATOOLS;

template <class Class>
Analysis_Object *
GetParticleSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"]
    .SetDefault(std::string(finalstate_list))
    .Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("Selected").Get<std::string>();
  return new Class(min,max,inlist,outlist);
}

#define DEFINE_SELECTOR_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter				\
  <Analysis_Object,Analysis_Key,CLASS>::			\
  operator()(const Analysis_Key& key) const		\
  { return GetParticleSelector<CLASS>(key); }

#define DEFINE_SELECTOR_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Analysis_Key,CLASS>::	\
  PrintInfo(std::ostream &str,const size_t width) const		\
  { str<<"e.g. {Min: 30, Max: 70, InList: FinalState, OutList: Selected}"; }

#define DEFINE_SELECTOR_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_SELECTOR_GETTER_METHOD(CLASS)				\
  DEFINE_SELECTOR_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

template <class Class>
Analysis_Object *
GetParticleDSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"]
    .SetDefault(std::string(finalstate_list))
    .Get<std::string>();
  const auto reflist = s["RefList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("Selected").Get<std::string>();
  const auto item = s["Item"].SetDefault(0).Get<size_t>();
  const auto kf = s["Flav"].SetDefault(kf_jet).Get<int>();
  ATOOLS::Flavour flav{ ATOOLS::Flavour((kf_code)(std::abs(kf))) };
  if (kf<0) flav=flav.Bar();

  return new Class(min,max,flav,item,reflist,inlist,outlist);
}

#define DEFINE_SELECTOR_D_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter				\
  <Analysis_Object,Analysis_Key,CLASS>::			\
  operator()(const Analysis_Key& key) const		\
  { return GetParticleDSelector<CLASS>(key); }

#define DEFINE_SELECTOR_D_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Analysis_Key,CLASS>::	\
  PrintInfo(std::ostream &str,const size_t width) const		\
  { str<<"e.g. {Min: 30, Max: 70, InList: FinalState, RefList: Jets, OutList: Selected, Item: 0, Flav: 93}"; }

#define DEFINE_SELECTOR_D_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_SELECTOR_D_GETTER_METHOD(CLASS)			\
  DEFINE_SELECTOR_D_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Particle_Selector_Base::
Particle_Selector_Base(const double min,const double max,
		       const std::string &inlist,const std::string &outlist):
  Trigger_Base(inlist,outlist)
{
  m_xmin=min;
  m_xmax=max;
}

void Particle_Selector_Base::Evaluate(const ATOOLS::Particle_List &inlist,
				      ATOOLS::Particle_List &outlist,
				      double value,double ncount)
{
  for (size_t i=0;i<inlist.size();++i) {
    if (Select(inlist[i])) 
      outlist.push_back(new ATOOLS::Particle(*inlist[i]));
  }
}

DEFINE_SELECTOR_GETTER(PT_Selector,"PTSel")

PT_Selector::PT_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool PT_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().PPerp());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *PT_Selector::GetCopy() const
{
  return new PT_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(ET_Selector,"ETSel")

ET_Selector::ET_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool ET_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().EPerp());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *ET_Selector::GetCopy() const
{
  return new ET_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(Eta_Selector,"EtaSel")

Eta_Selector::Eta_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool Eta_Selector::Select(const ATOOLS::Particle *p) const
{
  double eta(p->Momentum().Eta());
  return eta>=m_xmin && eta<=m_xmax;
}

Analysis_Object *Eta_Selector::GetCopy() const
{
  return new Eta_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(Abs_Eta_Selector,"AbsEtaSel")

Abs_Eta_Selector::Abs_Eta_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool Abs_Eta_Selector::Select(const ATOOLS::Particle *p) const
{
  double eta(p->Momentum().Eta());
  return dabs(eta)>=m_xmin && dabs(eta)<=m_xmax;
}

Analysis_Object *Abs_Eta_Selector::GetCopy() const
{
  return new Abs_Eta_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(Y_Selector,"YSel")

Y_Selector::Y_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool Y_Selector::Select(const ATOOLS::Particle *p) const
{
  double y(p->Momentum().Y());
  return y>=m_xmin && y<=m_xmax;
}

Analysis_Object *Y_Selector::GetCopy() const
{
  return new Y_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(Abs_Y_Selector,"AbsYSel")

Abs_Y_Selector::Abs_Y_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool Abs_Y_Selector::Select(const ATOOLS::Particle *p) const
{
  double y(p->Momentum().Y());
  return dabs(y)>=m_xmin && dabs(y)<=m_xmax;
}

Analysis_Object *Abs_Y_Selector::GetCopy() const
{
  return new Abs_Y_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_GETTER(Phi_Selector,"PhiSel")

Phi_Selector::Phi_Selector
(const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  Particle_Selector_Base(min,max,inlist,outlist) {}

bool Phi_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().Phi());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *Phi_Selector::GetCopy() const
{
  return new Phi_Selector(m_xmin,m_xmax,m_inlist,m_outlist);
}

DEFINE_SELECTOR_D_GETTER(DPhi_Selector,"DPhiSel")

DPhi_Selector::
DPhi_Selector(const double min,const double max,
	      const ATOOLS::Flavour flav,
	      const size_t item,const std::string &reflist,
	      const std::string &inlist,const std::string &outlist):
  Two_List_Trigger_Base(inlist,reflist,outlist),
  m_item(item), m_flavour(flav)
{
  m_xmin=min;
  m_xmax=max;
}

void DPhi_Selector::Evaluate(const ATOOLS::Particle_List &inlist,
			     const ATOOLS::Particle_List &reflist,
			     ATOOLS::Particle_List &outlist,
			     double value,double ncount)
{
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<reflist.size();++i) 
    if (reflist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==kf_none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  for (size_t i=0;i<inlist.size();++i) {
    double phi=
      ATOOLS::dabs(inlist[i]->Momentum().
		   DPhi(reflist[pos]->Momentum())/M_PI*180.0);
    if (phi>=m_xmin && phi<=m_xmax) 
      outlist.push_back(new ATOOLS::Particle(*inlist[i]));
  }
}

Analysis_Object *DPhi_Selector::GetCopy() const
{
  return new DPhi_Selector(m_xmin,m_xmax,m_flavour,m_item,
			   m_reflist,m_inlist,m_outlist);
}

