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
Analysis_Object *const 
GetTwoParticleLDeltaSelector(const Argument_Matrix &parameters) 
{
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(flav,ATOOLS::ToType<int>(parameters[0][1]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     parameters[0][4],parameters[0][5],parameters[0][6]);
  }
  if (parameters.size()<7) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", reflist="Jets", outlist="LeadJets";
  size_t item=0;
  ATOOLS::Flavour flav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
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
  return new Class(flav,item,min,max,inlist,reflist,outlist);
}									

#define DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  Analysis_Object *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetTwoParticleLDeltaSelector<CLASS>(parameters); }

#define DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav item min max inlist reflist outlist"; }

#define DEFINE_TWO_SELECTOR_DELTA_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(NAME)

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DPhiL_Selector,
				 Two_DPhiL_Selector_Getter,"TwoDPhiXSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DEtaL_Selector,
				 Two_DEtaL_Selector_Getter,"TwoDEtaXSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DYL_Selector,
				 Two_DYL_Selector_Getter,"TwoDYXSel")

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


DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PTL_Selector,
				 Two_PTL_Selector_Getter,"TwoPTXSel")

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


DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DRL_Selector,
				 Two_DRL_Selector_Getter,"TwoDRXSel")

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


