#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class STwo_Particle_Observable_Base: public Primitive_Observable_Base {  
  protected:

    std::string     m_reflist;
    ATOOLS::Flavour m_flavour, m_refflavour;

    size_t m_item, m_refitem;

  public:

    STwo_Particle_Observable_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const ATOOLS::Flavour ref,const size_t refitem,
     const int type,const double min,const double max,const int bins,
     const std::string &inlist,const std::string &reflist,
     const std::string &name);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    virtual bool Evaluate(const Particle *p1,const Particle *p2,
			  double weight=1.,double ncount=1) const = 0;

  };// end of class STwo_Particle_Observable_Base

  class Two_DPhi_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DPhi_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DPhi_Distribution

  class Two_DEta_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DEta_Distribution

  class Two_PEta_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PEta_Distribution

  class Two_DY_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DY_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DY_Distribution

  class Two_PY_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PY_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PY_Distribution

  class Two_Mass_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_Mass_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_Mass_Distribution

  class Two_PT_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PT_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PT_Distribution

  class Two_DR_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DR_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DR_Distribution

  class Two_ETFrac_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_ETFrac_Distribution(const ATOOLS::Flavour flav,const size_t item,
			    const ATOOLS::Flavour ref,const size_t refitem,
			    const int type,const double min,const double max,const int bins,
			    const std::string &inlist,const std::string &reflist);
    
    bool Evaluate(const Particle *p1,const Particle *p2,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_ETFrac_Distribution

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *const 
GetSTwoParticleObservable(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<8) return NULL;
    std::string list(parameters[0].size()>8?parameters[0][8]:"FinalState");
    std::string rlist(parameters[0].size()>9?parameters[0][9]:list);
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    kf=ATOOLS::ToType<int>(parameters[0][2]);
    ATOOLS::Flavour refflav((kf_code)abs(kf));
    if (kf<0) refflav=refflav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     refflav,ATOOLS::ToType<size_t>(parameters[0][3]),
		     HistogramType(parameters[0][7]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][5]),
		     ATOOLS::ToType<int>(parameters[0][6]),list,rlist);
  }
  if (parameters.size()<9) return NULL;
  int bins=100, scale=0;
  double min=30.0, max=70.0;
  std::string inlist="Jets", reflist="Jets";
  size_t item=0, refitem=1;
  ATOOLS::Flavour flav(kf_jet), refflav(kf_jet);
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    else if (parameters[i][0]=="InList") inlist=parameters[i][1];
    else if (parameters[i][0]=="RefList") reflist=parameters[i][1];
    else if (parameters[i][0]=="Min") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Max") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="Bins") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Item1") item=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Item2") refitem=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="Scale") scale=HistogramType(parameters[i][1]);
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
  return new Class(flav,item,refflav,refitem,scale,min,max,bins,inlist,reflist);
}									

#define DEFINE_TWO_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetSTwoParticleObservable<CLASS>(parameters); }

#define DEFINE_TWO_OBSERVABLE_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav1 item1 flav2 item2 min max bins Lin|LinErr|Log|LogErr [inlist [reflist]]"; }

#define DEFINE_TWO_OBSERVABLE_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_TWO_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_TWO_OBSERVABLE_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

STwo_Particle_Observable_Base::
STwo_Particle_Observable_Base(const ATOOLS::Flavour flav,const size_t item,
			      const ATOOLS::Flavour refflav,const size_t refitem,
			      const int type,const double min,const double max,const int bins,
			      const std::string &inlist,const std::string &reflist,
			      const std::string &name):
  Primitive_Observable_Base(type,min,max,bins), 
  m_reflist(reflist),
  m_flavour(flav),
  m_refflavour(refflav),
  m_item(item),
  m_refitem(refitem)
{
  m_listname=inlist;
  m_name=name+"_"+ToString(m_flavour)+"-"+ToString(m_item)+"_"
    +ToString(m_refflavour)+"-"+ToString(m_refitem)+".dat";
}

void STwo_Particle_Observable_Base::Evaluate(const ATOOLS::Particle_List &list,
					     double weight,double ncount)
{
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<list.size();++i) {
    if (list[i]->Flav()==m_flavour || 
	m_flavour.Kfcode()==kf_none) {
      ++no;
      if (no==(int)m_item) {
	pos=i;
	if (refpos!=std::string::npos) break;
      }
    }
  }
  for (size_t i=0;i<reflist->size();++i) {
    if ((*reflist)[i]->Flav()==m_refflavour || 
	m_refflavour.Kfcode()==kf_none) {
      ++refno;
      if (refno==(int)m_refitem) {
	refpos=i;
	if (pos!=std::string::npos) break;
      }
    }
  }
  if (pos==std::string::npos || refpos==std::string::npos) return;
  Evaluate(list[pos],(*reflist)[refpos],weight,ncount);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_DPhi_Distribution,
			     Two_DPhi_Distribution_Getter,"TwoDPhi")

  Two_DPhi_Distribution::
Two_DPhi_Distribution(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour refflav,const size_t refitem,
		      const int type,const double min,const double max,const int bins,
		      const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoDPhi") {}

bool Two_DPhi_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				     double weight,double ncount) const
{
  p_histo->Insert(p1->Momentum().DPhi(p2->Momentum()),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_DPhi_Distribution::Copy() const
{
  return new Two_DPhi_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				   m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_DEta_Distribution,
			     Two_DEta_Distribution_Getter,"TwoDEta")

  Two_DEta_Distribution::
Two_DEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour refflav,const size_t refitem,
		      const int type,const double min,const double max,const int bins,
		      const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoDEta") {}

bool Two_DEta_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				     double weight,double ncount) const
{
  p_histo->Insert(dabs(p1->Momentum().Eta()-p2->Momentum().Eta()),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_DEta_Distribution::Copy() const
{
  return new Two_DEta_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				   m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_PEta_Distribution,
			     Two_PEta_Distribution_Getter,"TwoPEta")
  
  Two_PEta_Distribution::
Two_PEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour refflav,const size_t refitem,
		      const int type,const double min,const double max,const int bins,
		      const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoPEta") {}

bool Two_PEta_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				     double weight,double ncount) const
{
  p_histo->Insert(p1->Momentum().Eta()*p2->Momentum().Eta(),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_PEta_Distribution::Copy() const
{
  return new Two_PEta_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				   m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_DY_Distribution,
			     Two_DY_Distribution_Getter,"TwoDY")

  Two_DY_Distribution::
Two_DY_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour refflav,const size_t refitem,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoDY") {}

bool Two_DY_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				   double weight,double ncount) const
{
  p_histo->Insert(dabs(p1->Momentum().Y()-p2->Momentum().Y()),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_DY_Distribution::Copy() const
{
  return new Two_DY_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_PY_Distribution,
			     Two_PY_Distribution_Getter,"TwoPY")

  Two_PY_Distribution::
Two_PY_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour refflav,const size_t refitem,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoPY") {}

bool Two_PY_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				   double weight,double ncount) const
{
  p_histo->Insert(p1->Momentum().Y()*p2->Momentum().Y(),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_PY_Distribution::Copy() const
{
  return new Two_PY_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_Mass_Distribution,
			     Two_Mass_Distribution_Getter,"TwoMass")

  Two_Mass_Distribution::
Two_Mass_Distribution(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour refflav,const size_t refitem,
		      const int type,const double min,const double max,const int bins,
		      const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoMass") {}

bool Two_Mass_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				     double weight,double ncount) const
{
  p_histo->Insert((p1->Momentum()+p2->Momentum()).Mass(),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_Mass_Distribution::Copy() const
{
  return new Two_Mass_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				   m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_PT_Distribution,
			     Two_PT_Distribution_Getter,"TwoPT")

  Two_PT_Distribution::
Two_PT_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour refflav,const size_t refitem,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoPT") {}

bool Two_PT_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				   double weight,double ncount) const
{
  p_histo->Insert((p1->Momentum()+p2->Momentum()).PPerp(),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_PT_Distribution::Copy() const
{
  return new Two_PT_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_DR_Distribution,
			     Two_DR_Distribution_Getter,"TwoDR")

  Two_DR_Distribution::
Two_DR_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour refflav,const size_t refitem,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoDR") {}

bool Two_DR_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				   double weight,double ncount) const
{
  p_histo->Insert(sqrt(sqr(p1->Momentum().Eta()-p2->Momentum().Eta())+
		       sqr(p1->Momentum().DPhi(p2->Momentum()))),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_DR_Distribution::Copy() const
{
  return new Two_DR_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_ETFrac_Distribution,
			     Two_ETFrac_Distribution_Getter,"TwoETFrac")

  Two_ETFrac_Distribution::
Two_ETFrac_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour refflav,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoETFrac") {}

bool Two_ETFrac_Distribution::Evaluate(const Particle *p1,const Particle *p2,
				       double weight,double ncount) const
{
  p_histo->Insert(p1->Momentum().EPerp()/p2->Momentum().EPerp(),weight,ncount);
  return true;
}

Primitive_Observable_Base *Two_ETFrac_Distribution::Copy() const
{
  return new Two_ETFrac_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				     m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}


	  
