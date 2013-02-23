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
GetJetParticleDeltaSelector(const Argument_Matrix &parameters) 
{
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     parameters[0][4],parameters[0][5]);
  }
  if (parameters.size()<6) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", outlist="LeadJets";
  size_t item=0;
  ATOOLS::Flavour flav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour((kf_code)(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(flav,item,min,max,inlist,outlist);
}									

#define DEFINE_SELECTOR_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter				\
  <Analysis_Object,Argument_Matrix,CLASS>::			\
  operator()(const Argument_Matrix &parameters) const		\
  { return GetJetParticleDeltaSelector<CLASS>(parameters); }

#define DEFINE_SELECTOR_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Argument_Matrix,CLASS>::	\
  PrintInfo(std::ostream &str,const size_t width) const		\
  { str<<"flav item min max inlist outlist"; }

#define DEFINE_JET_SELECTOR_DELTA_GETTER(CLASS,TAG)		\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Argument_Matrix);	\
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


