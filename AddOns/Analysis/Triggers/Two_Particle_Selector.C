#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class Two_Particle_Selector_Base: public Two_List_Trigger_Base {  
  protected:

    ATOOLS::Flavour m_flavour, m_refflavour;

    double m_xmin, m_xmax;
    size_t m_item, m_refitem;

  public:

    Two_Particle_Selector_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const ATOOLS::Flavour ref,const size_t refitem,
     const double min,const double max,
     const std::string &inlist,const std::string &reflist,
     const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &inlist,
		  const ATOOLS::Particle_List &reflist,
		  ATOOLS::Particle_List &outlist,
		  double value,double ncount);
    
    virtual bool Select(const Particle *p1,const Particle *p2) const = 0;

  };// end of class Two_Particle_Selector_Base

  class Two_DPhi_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DPhi_Selector

  class Two_DEta_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_DEta_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DEta_Selector

  class Two_PEta_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_PEta_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_PEta_Selector

  class Two_DY_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_DY_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DY_Selector

  class Two_PY_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_PY_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_PY_Selector

  class Two_CY_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_CY_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_CY_Selector

  class Two_Mass_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_Mass_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_Mass_Selector

  class Two_MT_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_MT_Selector(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour ref,const size_t refitem,
		    const double min,const double max,
		    const std::string &inlist,const std::string &reflist,
		    const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_MT_Selector

  class Two_PT_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_PT_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_PT_Selector

  class Two_DR_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_DR_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_DR_Selector

  class Two_ETFrac_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_ETFrac_Selector(const ATOOLS::Flavour flav,const size_t item,
		       const ATOOLS::Flavour ref,const size_t refitem,
		       const double min,const double max,
		       const std::string &inlist,const std::string &reflist,
		       const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_ETFrac_Selector

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Analysis_Object *const 
GetTwoParticleDeltaSelector(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<9) return NULL;
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    kf=ATOOLS::ToType<int>(parameters[0][2]);
    ATOOLS::Flavour refflav((kf_code)abs(kf));
    if (kf<0) refflav=refflav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     refflav,ATOOLS::ToType<size_t>(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][5]),
		     parameters[0][6],parameters[0][7],parameters[0][8]);
  }
  if (parameters.size()<9) return NULL;
  double min=30.0, max=70.0;
  std::string inlist="Jets", reflist="Jets", outlist="LeadJets";
  size_t item=0, refitem=1;
  ATOOLS::Flavour flav(kf_jet), refflav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="OutList") outlist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Item1") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Item2") refitem=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Flav1") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      flav=ATOOLS::Flavour((kf_code)(abs(kf)));
      if (kf<0) flav=flav.Bar();
    }
    else if (parameters[i][0]=="Flav2") {
      int kf=ATOOLS::ToType<int>(parameters[i][1]);
      refflav=ATOOLS::Flavour((kf_code)(abs(kf)));
      if (kf<0) refflav=refflav.Bar();
    }
  }
  return new Class(flav,item,refflav,refitem,min,max,inlist,reflist,outlist);
}									

#define DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  Analysis_Object *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetTwoParticleDeltaSelector<CLASS>(parameters); }

#define DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 flav2 item2 min max inlist reflist outlist"; }

#define DEFINE_TWO_SELECTOR_DELTA_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Analysis_Object,Argument_Matrix);	\
  DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

Two_Particle_Selector_Base::
Two_Particle_Selector_Base(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_List_Trigger_Base(inlist,reflist,outlist),
  m_flavour(flav), m_refflavour(refflav),
  m_item(item), m_refitem(refitem)
{
  m_xmin=min;
  m_xmax=max;
}

void Two_Particle_Selector_Base::Evaluate
(const ATOOLS::Particle_List &inlist,const ATOOLS::Particle_List &reflist,
 ATOOLS::Particle_List &outlist,double value,double ncount)
{
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<reflist.size();++i) {
    if (reflist[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==kf_none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	if (refpos!=std::string::npos) break;
      }
    }
    if (reflist[i]->Flav()==m_refflavour || 
	m_refflavour.Kfcode()==kf_none) {
      ++refno;
      if (refno==(int)m_refitem) {
	refpos=i;
	if (pos!=std::string::npos) break;
      }
    }
  }
  if (pos==std::string::npos || refpos==std::string::npos) return;
  if (Select(reflist[pos],reflist[refpos])) {
    outlist.resize(inlist.size());
    for (size_t i=0;i<inlist.size();++i) 
      outlist[i] = new ATOOLS::Particle(*inlist[i]);
  }
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DPhi_Selector,
				 Two_DPhi_Selector_Getter,"TwoDPhiSel")

Two_DPhi_Selector::
Two_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_DPhi_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double dphi(p1->Momentum().DPhi(p2->Momentum())/M_PI*180.0);
  if (dphi<m_xmin || dphi>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DPhi_Selector::GetCopy() const
{
  return new Two_DPhi_Selector(m_flavour,m_item,m_refflavour,m_refitem,
			       m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DEta_Selector,
				 Two_DEta_Selector_Getter,"TwoDEtaSel")

Two_DEta_Selector::
Two_DEta_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_DEta_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double deta=dabs(p1->Momentum().Eta()-p2->Momentum().Eta());
  if (deta<m_xmin || deta>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DEta_Selector::GetCopy() const
{
  return new Two_DEta_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PEta_Selector,
				 Two_PEta_Selector_Getter,"TwoPEtaSel")

Two_PEta_Selector::
Two_PEta_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_PEta_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double peta=p1->Momentum().Eta()*p2->Momentum().Eta();
  if (peta<m_xmin || peta>m_xmax) return false;
  return true;
}

Analysis_Object *Two_PEta_Selector::GetCopy() const
{
  return new Two_PEta_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DY_Selector,
				 Two_DY_Selector_Getter,"TwoDYSel")

Two_DY_Selector::
Two_DY_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_DY_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double dy=dabs(p1->Momentum().Y()-p2->Momentum().Y());
  if (dy<m_xmin || dy>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DY_Selector::GetCopy() const
{
  return new Two_DY_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PY_Selector,
				 Two_PY_Selector_Getter,"TwoPYSel")

Two_PY_Selector::
Two_PY_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_PY_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double py=p1->Momentum().Y()*p2->Momentum().Y();
  if (py<m_xmin || py>m_xmax) return false;
  return true;
}

Analysis_Object *Two_PY_Selector::GetCopy() const
{
  return new Two_PY_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_CY_Selector,
				 Two_CY_Selector_Getter,"TwoCYSel")

Two_CY_Selector::
Two_CY_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_CY_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double cy=dabs((p1->Momentum()+p2->Momentum()).Y());
  if (cy<m_xmin || cy>m_xmax) return false;
  return true;
}

Analysis_Object *Two_CY_Selector::GetCopy() const
{
  return new Two_CY_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_Mass_Selector,
				 Two_Mass_Selector_Getter,"TwoMassSel")

Two_Mass_Selector::
Two_Mass_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_Mass_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double mass=(p1->Momentum()+p2->Momentum()).Mass();
  if (mass<m_xmin || mass>m_xmax) return false;
  return true;
}

Analysis_Object *Two_Mass_Selector::GetCopy() const
{
  return new Two_Mass_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_MT_Selector,
				 Two_MT_Selector_Getter,"TwoMTSel")

Two_MT_Selector::
Two_MT_Selector(const ATOOLS::Flavour flav,const size_t item,
		const ATOOLS::Flavour refflav,const size_t refitem,
		const double min,const double max,
		const std::string &inlist,const std::string &reflist,
		const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_MT_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double mt=(p1->Momentum()+p2->Momentum()).MPerp();
  if (mt<m_xmin || mt>m_xmax) return false;
  return true;
}

Analysis_Object *Two_MT_Selector::GetCopy() const
{
  return new Two_MT_Selector(m_flavour,m_item,m_refflavour,m_refitem,
			     m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PT_Selector,
				 Two_PT_Selector_Getter,"TwoPTSel")

Two_PT_Selector::
Two_PT_Selector(const ATOOLS::Flavour flav,const size_t item,
		const ATOOLS::Flavour refflav,const size_t refitem,
		const double min,const double max,
		const std::string &inlist,const std::string &reflist,
		const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_PT_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double pt=(p1->Momentum()+p2->Momentum()).PPerp();
  if (pt<m_xmin || pt>m_xmax) return false;
  return true;
}

Analysis_Object *Two_PT_Selector::GetCopy() const
{
  return new Two_PT_Selector(m_flavour,m_item,m_refflavour,m_refitem,
			     m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DR_Selector,
				 Two_DR_Selector_Getter,"TwoDRSel")

Two_DR_Selector::
Two_DR_Selector(const ATOOLS::Flavour flav,const size_t item,
		  const ATOOLS::Flavour refflav,const size_t refitem,
		  const double min,const double max,
		  const std::string &inlist,const std::string &reflist,
		  const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_DR_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double dr=(sqr(p1->Momentum().Eta()-p2->Momentum().Eta())+
	     sqr(p1->Momentum().DPhi(p2->Momentum())));
  if (dr<m_xmin || dr>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DR_Selector::GetCopy() const
{
  return new Two_DR_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_ETFrac_Selector,
				 Two_ETFrac_Selector_Getter,"TwoETFracSel")

Two_ETFrac_Selector::
Two_ETFrac_Selector(const ATOOLS::Flavour flav,const size_t item,
		   const ATOOLS::Flavour refflav,const size_t refitem,
		   const double min,const double max,
		   const std::string &inlist,const std::string &reflist,
		   const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_ETFrac_Selector::Select(const Particle *p1,const Particle *p2) const
{
  double efrac=p1->Momentum().EPerp()/p2->Momentum().EPerp();
  if (efrac<m_xmin || efrac>m_xmax) return false;
  return true;
}

Analysis_Object *Two_ETFrac_Selector::GetCopy() const
{
  return new Two_ETFrac_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				 m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}


