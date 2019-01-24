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

  class Two_MT2_Selector: public Two_Particle_Selector_Base {  
  public:

    Two_MT2_Selector(const ATOOLS::Flavour flav,const size_t item,
		     const ATOOLS::Flavour ref,const size_t refitem,
		     const double min,const double max,
		     const std::string &inlist,const std::string &reflist,
		     const std::string &outlist);
    
    bool Select(const Particle *p1,const Particle *p2) const;

    Analysis_Object *GetCopy() const;
    
  };// end of class Two_MT2_Selector

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
Analysis_Object *
GetTwoParticleDeltaSelector(const Analysis_Key& key)
{
  ATOOLS::Scoped_Settings s{ key.m_settings };
  const auto min = s["Min"].SetDefault(30.0).Get<double>();
  const auto max = s["Max"].SetDefault(70.0).Get<double>();
  const auto inlist = s["InList"].SetDefault("Jets").Get<std::string>();
  const auto reflist = s["RefList"].SetDefault("Jets").Get<std::string>();
  const auto outlist = s["OutList"].SetDefault("LeadJets").Get<std::string>();
  const auto item = s["Item1"].SetDefault(0).Get<size_t>();
  const auto refitem = s["Item2"].SetDefault(1).Get<size_t>();
  std::vector<ATOOLS::Flavour> flavs;
  flavs.reserve(2);
  for (size_t i{ 0 }; i < 2; ++i) {
    const auto flavkey = "Flav" + ATOOLS::ToString(i + 1);
    const auto kf = s[flavkey].SetDefault(kf_jet).Get<int>();
    flavs.push_back(ATOOLS::Flavour((kf_code)std::abs(kf)));
    if (kf < 0)
      flavs.back() = flavs.back().Bar();
  }
  return new Class(flavs[0],item,flavs[2],refitem,min,max,inlist,reflist,outlist);
}

#define DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS)			\
  Analysis_Object *ATOOLS::Getter					\
  <Analysis_Object,Analysis_Key,CLASS>:: 				\
  operator()(const Analysis_Key& key) const			\
  { return GetTwoParticleDeltaSelector<CLASS>(key); }

#define DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(CLASS)			\
  void ATOOLS::Getter<Analysis_Object,Analysis_Key,CLASS>::		\
  PrintInfo(std::ostream &str,const size_t width) const			\
  { str<<"e.g. {Min: 30, Max: 70, InList: Jets, RefList: Jets, OutList: LeadJets, Item1: 0, Item2: 1, Flav1: 93, Flav1: 93}"; }

#define DEFINE_TWO_SELECTOR_DELTA_GETTER(CLASS,TAG)		\
  DECLARE_GETTER(CLASS,TAG,Analysis_Object,Analysis_Key);	\
  DEFINE_TWO_SELECTOR_DELTA_GETTER_METHOD(CLASS)		\
  DEFINE_TWO_SELECTOR_DELTA_PRINT_METHOD(CLASS)

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DPhi_Selector,"TwoDPhiSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DEta_Selector,"TwoDEtaSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PEta_Selector,"TwoPEtaSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DY_Selector,"TwoDYSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PY_Selector,"TwoPYSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_CY_Selector,"TwoCYSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_Mass_Selector,"TwoMassSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_MT_Selector,"TwoMTSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_MT2_Selector,"TwoMT2Sel")

Two_MT2_Selector::
Two_MT2_Selector(const ATOOLS::Flavour flav,const size_t item,
		 const ATOOLS::Flavour refflav,const size_t refitem,
		 const double min,const double max,
		 const std::string &inlist,const std::string &reflist,
		 const std::string &outlist):
  Two_Particle_Selector_Base(flav,item,refflav,refitem,min,max,
			     inlist,reflist,outlist) {}

bool Two_MT2_Selector::Select(const Particle *p1,const Particle *p2) const
{
  Vec4D mom1 = p1->Momentum(), mom2 = p2->Momentum();
  double mass = sqrt(2.*(mom1.PPerp()*mom2.PPerp()-mom1[1]*mom2[1]-mom1[2]*mom2[2]));
  if (mass<m_xmin || mass>m_xmax) return false;
  return true;
}

Analysis_Object *Two_MT2_Selector::GetCopy() const
{
  return new Two_MT2_Selector(m_flavour,m_item,m_refflavour,m_refitem,
			      m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_PT_Selector,"TwoPTSel")

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

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_DR_Selector,"TwoDRSel")

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
  double dr=sqrt(sqr(p1->Momentum().Eta()-p2->Momentum().Eta())+
	     sqr(p1->Momentum().DPhi(p2->Momentum())));
  if (dr<m_xmin || dr>m_xmax) return false;
  return true;
}

Analysis_Object *Two_DR_Selector::GetCopy() const
{
  return new Two_DR_Selector(m_flavour,m_item,m_refflavour,m_refitem,
				   m_xmin,m_xmax,m_inlist,m_reflist,m_outlist);
}

DEFINE_TWO_SELECTOR_DELTA_GETTER(Two_ETFrac_Selector,"TwoETFracSel")

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


