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
Analysis_Object *const 
GetParticleSelector(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    return new Class(ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     parameters[0][2],parameters[0][3]);
  }
  if (parameters.size()<4) return NULL;
  double min=30.0, max=70.0;
  std::string inlist=finalstate_list, outlist="Selected";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
  }
  return new Class(min,max,inlist,outlist);
}									

#define DEFINE_SELECTOR_GETTER_METHOD(CLASS,NAME)		\
  Analysis_Object *				\
  NAME::operator()(const Argument_Matrix &parameters) const	\
  { return GetParticleSelector<CLASS>(parameters); }

#define DEFINE_SELECTOR_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max inlist outlist"; }

#define DEFINE_SELECTOR_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_SELECTOR_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_SELECTOR_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

template <class Class>
Analysis_Object *const 
GetParticleDSelector(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][3]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     flav,ATOOLS::ToType<int>(parameters[0][3]),
		     parameters[0][5],parameters[0][4],parameters[0][6]);
  }
  if (parameters.size()<7) return NULL;
  double min=30.0, max=70.0;
  std::string inlist=finalstate_list, reflist="Jets", outlist="Selected";
  size_t item=0;
  ATOOLS::Flavour flav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour((kf_code)(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
  }
  return new Class(min,max,flav,item,reflist,inlist,outlist);
}									

#define DEFINE_SELECTOR_D_GETTER_METHOD(CLASS,NAME)		\
  Analysis_Object *				\
  NAME::operator()(const Argument_Matrix &parameters) const	\
  { return GetParticleDSelector<CLASS>(parameters); }

#define DEFINE_SELECTOR_D_PRINT_METHOD(NAME)				\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max flav item inlist reflist outlist"; }

#define DEFINE_SELECTOR_D_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_SELECTOR_D_GETTER_METHOD(CLASS,NAME)			\
  DEFINE_SELECTOR_D_PRINT_METHOD(NAME)

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

DEFINE_SELECTOR_GETTER(PT_Selector,PT_Selector_Getter,"PTSel")

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

DEFINE_SELECTOR_GETTER(ET_Selector,ET_Selector_Getter,"ETSel")

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

DEFINE_SELECTOR_GETTER(Eta_Selector,Eta_Selector_Getter,"EtaSel")

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

DEFINE_SELECTOR_GETTER(Abs_Eta_Selector,Abs_Eta_Selector_Getter,"AbsEtaSel")

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

DEFINE_SELECTOR_GETTER(Y_Selector,Y_Selector_Getter,"YSel")

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

DEFINE_SELECTOR_GETTER(Abs_Y_Selector,Abs_Y_Selector_Getter,"AbsYSel")

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

DEFINE_SELECTOR_GETTER(Phi_Selector,Phi_Selector_Getter,"PhiSel")

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

DEFINE_SELECTOR_D_GETTER(DPhi_Selector,DPhi_Selector_Getter,"DPhiSel")

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

