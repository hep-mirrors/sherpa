#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

namespace ANALYSIS {

  class One_Particle_Selector_Base: public Trigger_Base {  
  protected:

    ATOOLS::Flavour m_flavour;
    
    double m_xmin, m_xmax;

    size_t       m_item;
    unsigned int m_mode;
    
  public:

    One_Particle_Selector_Base
    (const ATOOLS::Flavour flav,const size_t item,const int mode,
     const double min,const double max,
     const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);

    virtual bool Select(const ATOOLS::Particle *p) const = 0;
    
  };// end of class One_Particle_Selector_Base

  class One_PT_Selector: public One_Particle_Selector_Base {  
  public:

    One_PT_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,
		    const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class One_PT_Selector

  class One_ET_Selector: public One_Particle_Selector_Base {  
  public:

    One_ET_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,
		    const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class One_ET_Selector

  class One_Eta_Selector: public One_Particle_Selector_Base {  
  public:

    One_Eta_Selector(const ATOOLS::Flavour flav,
		     const size_t item,const int mode,
		     const double min,const double max,
		     const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class One_Eta_Selector

  class One_Y_Selector: public One_Particle_Selector_Base {  
  public:

    One_Y_Selector(const ATOOLS::Flavour flav,
		   const size_t item,const int mode,
		   const double min,const double max,
		   const std::string &inlist,const std::string &outlist);
    
    bool Select(const ATOOLS::Particle *p) const;
    
    Analysis_Object *GetCopy() const;
    
  };// end of class One_Y_Selector

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *
GetOneParticleSelector(const Argument_Matrix &parameters) 
{				
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     parameters[0][5],parameters[0][6]);
  }
  if (parameters.size()<7) return NULL;
  double min=30.0, max=70.0; 
  std::string inlist="Jets", outlist="LeadJets";
  size_t item=0;
  int mode=0;
  ATOOLS::Flavour flav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Mode") mode=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour((kf_code)(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(flav,item,mode,min,max,inlist,outlist);
}									

#define DEFINE_ONE_SELECTOR_GETTER_METHOD(CLASS)		\
  Analysis_Object *ATOOLS::Getter				\
  <Analysis_Object,Argument_Matrix,CLASS>::			\
  operator()(const Argument_Matrix &parameters) const		\
  { return GetOneParticleSelector<CLASS>(parameters); }

#define DEFINE_ONE_SELECTOR_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Argument_Matrix,CLASS>::	\
  PrintInfo(std::ostream &str,const size_t width) const		\
  { str<<"flav item mode min max inlist outlist"; }

#define DEFINE_ONE_SELECTOR_GETTER(CLASS,TAG)			\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_ONE_SELECTOR_GETTER_METHOD(CLASS)			\
  DEFINE_ONE_SELECTOR_PRINT_METHOD(CLASS)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

One_Particle_Selector_Base::
One_Particle_Selector_Base(const ATOOLS::Flavour flav,
			   const size_t item,const int mode,
			   const double min,const double max,
			   const std::string &inlist,
			   const std::string &outlist):
  Trigger_Base(inlist,outlist),
  m_flavour(flav), m_item(item), m_mode(mode)
{
  m_xmin=min;
  m_xmax=max;
}

void One_Particle_Selector_Base::Evaluate(const ATOOLS::Particle_List &inlist,
					  ATOOLS::Particle_List &outlist,
					  double weight,double ncount)
{
  int no=-1; 
  size_t pos=std::string::npos;
  for (size_t i=0;i<inlist.size();++i) 
    if (inlist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==kf_none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	break;
      }
    }
  if (pos==std::string::npos) return;
  if (!Select(inlist[pos])) return;
  if (m_mode==0) { 
    outlist.resize(inlist.size());
    for (size_t i=0;i<inlist.size();++i) 
      outlist[i] = new ATOOLS::Particle(*inlist[i]);
  }
  else {
    int diff_flavour=0;
    for (size_t i=0;i<inlist.size();++i) {
      if (inlist[i]->Flav()!=m_flavour) ++diff_flavour;
    }
    int size = diff_flavour+1; 
    outlist.resize(size);
    for (size_t i=0;i<inlist.size();++i) {
      if (inlist[i]->Flav()!=m_flavour || i==pos) {
	outlist[size-1] = new ATOOLS::Particle(*inlist[i]);
	--size;
      }
    }
  }
}

DEFINE_ONE_SELECTOR_GETTER(One_PT_Selector,"OnePTSel")

One_PT_Selector::One_PT_Selector
(const ATOOLS::Flavour flav,const size_t item,const int mode,
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  One_Particle_Selector_Base(flav,item,mode,min,max,inlist,outlist) {}

bool One_PT_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().PPerp());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *One_PT_Selector::GetCopy() const
{
  return new One_PT_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,
			     m_inlist,m_outlist);
}

DEFINE_ONE_SELECTOR_GETTER(One_ET_Selector,"OneETSel")

One_ET_Selector::One_ET_Selector
(const ATOOLS::Flavour flav,const size_t item,const int mode,
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  One_Particle_Selector_Base(flav,item,mode,min,max,inlist,outlist) {}

bool One_ET_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().EPerp());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *One_ET_Selector::GetCopy() const
{
  return new One_ET_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,
			     m_inlist,m_outlist);
}

DEFINE_ONE_SELECTOR_GETTER(One_Eta_Selector,"OneEtaSel")

One_Eta_Selector::One_Eta_Selector
(const ATOOLS::Flavour flav,const size_t item,const int mode,
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  One_Particle_Selector_Base(flav,item,mode,min,max,inlist,outlist) {}

bool One_Eta_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().Eta());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *One_Eta_Selector::GetCopy() const
{
  return new One_Eta_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,
			      m_inlist,m_outlist);
}

DEFINE_ONE_SELECTOR_GETTER(One_Y_Selector,"OneYSel")

One_Y_Selector::One_Y_Selector
(const ATOOLS::Flavour flav,const size_t item,const int mode,
 const double min,const double max,
 const std::string &inlist,const std::string &outlist):
  One_Particle_Selector_Base(flav,item,mode,min,max,inlist,outlist) {}

bool One_Y_Selector::Select(const ATOOLS::Particle *p) const
{
  double pt(p->Momentum().Y());
  return pt>=m_xmin && pt<=m_xmax;
}

Analysis_Object *One_Y_Selector::GetCopy() const
{
  return new One_Y_Selector(m_flavour,m_item,m_mode,m_xmin,m_xmax,
			    m_inlist,m_outlist);
}

