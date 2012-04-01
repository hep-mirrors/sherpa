#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Histogram_2D.H"
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
		  double weight=1.,double ncount=1.);
    
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    
    virtual double Calc(const Particle *p1,const Particle *p2) = 0;

  };// end of class STwo_Particle_Observable_Base

  class STwo2D_Particle_Observable_Base: public Primitive_Observable_Base {  
  protected:

    std::string     m_reflist;
    ATOOLS::Flavour m_flavour, m_refflavour;

    size_t m_item, m_refitem;
    ATOOLS::Histogram_2D* p_2dhisto;

  public:

    STwo2D_Particle_Observable_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const ATOOLS::Flavour ref,const size_t refitem,
     const int type,const double min,const double max,const int bins,
     const std::string &inlist,const std::string &reflist,
     const std::string &name);
    ~STwo2D_Particle_Observable_Base();

    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();

    void EndEvaluation(double scale=1.);
    void Reset();
    void Restore(double scale=1.0);
    void Output(const std::string & pname);
    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base & ob);
    
    virtual double Calc1(const Particle *p) = 0;
    virtual double Calc2(const Particle *p) = 0;

  };// end of class STwo2D_Particle_Observable_Base

  class Two_DPhi_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DPhi_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DPhi_Distribution

  class Two_DEta_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DEta_Distribution

  class Two_PEta_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PEta_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PEta_Distribution

  class Two_DY_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DY_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DY_Distribution

  class Two_PY_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PY_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PY_Distribution

  class Two_Mass_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_Mass_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const ATOOLS::Flavour ref,const size_t refitem,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_Mass_Distribution

  class Two_PT_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_PT_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PT_Distribution

  class Two_DPT_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DPT_Distribution(const ATOOLS::Flavour flav,const size_t item,
			 const ATOOLS::Flavour ref,const size_t refitem,
			 const int type,const double min,const double max,const int bins,
			 const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DPT_Distribution

  class Two_RPT_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_RPT_Distribution(const ATOOLS::Flavour flav,const size_t item,
			 const ATOOLS::Flavour ref,const size_t refitem,
			 const int type,const double min,const double max,const int bins,
			 const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_RPT_Distribution

  class Two_DR_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_DR_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_DR_Distribution

  class Two_ETFrac_Distribution: public STwo_Particle_Observable_Base {  
  public:

    Two_ETFrac_Distribution(const ATOOLS::Flavour flav,const size_t item,
			    const ATOOLS::Flavour ref,const size_t refitem,
			    const int type,const double min,const double max,const int bins,
			    const std::string &inlist,const std::string &reflist);
    
    double Calc(const Particle *p1,const Particle *p2);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_ETFrac_Distribution

  class Two_PT2D_Distribution: public STwo2D_Particle_Observable_Base {  
  public:

    Two_PT2D_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const ATOOLS::Flavour ref,const size_t refitem,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist,const std::string &reflist);
    
    double Calc1(const Particle *p);
    double Calc2(const Particle *p);

    Primitive_Observable_Base *Copy() const;
    
  };// end of class Two_PT_Distribution

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *
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
#include "ATOOLS/Org/Shell_Tools.H"

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
  if (pos==std::string::npos || refpos==std::string::npos) {
    p_histo->Insert(0.,0.,ncount);
    return;
  }
  p_histo->Insert(Calc(list[pos],(*reflist)[refpos]),weight,ncount);
}

void STwo_Particle_Observable_Base::EvaluateNLOcontrib(double weight,double ncount)
{
  ATOOLS::Particle_List * plist =p_ana->GetParticleList(m_listname);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<plist->size();++i) {
    if ((*plist)[i]->Flav()==m_flavour || 
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
  if (pos==std::string::npos || refpos==std::string::npos) {
    p_histo->InsertMCB(0.,0.,ncount);
    return;
  }
  p_histo->InsertMCB(Calc((*plist)[pos],(*reflist)[refpos]),weight,ncount);
}

void STwo_Particle_Observable_Base::EvaluateNLOevt()
{
  p_histo->FinishMCB();
}


STwo2D_Particle_Observable_Base::
STwo2D_Particle_Observable_Base(const ATOOLS::Flavour flav,const size_t item,
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
  int type2d(m_type);
  if ((type2d%1000)/100==1) type2d+=900;
  if ((type2d%100)/10==1) type2d+=100*((type2d%100)/10);
  p_2dhisto = new Histogram_2D(type2d,m_xmin,m_xmax,m_nbins,m_xmin,m_xmax,m_nbins);
  m_listname=inlist;
  m_name=name+"_"+ToString(m_flavour)+"-"+ToString(m_item)+"_"
    +ToString(m_refflavour)+"-"+ToString(m_refitem)+".dat";
}

STwo2D_Particle_Observable_Base::~STwo2D_Particle_Observable_Base()
{
  delete p_2dhisto;
}

void STwo2D_Particle_Observable_Base::Evaluate(const ATOOLS::Particle_List &list,
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
  if (pos==std::string::npos || refpos==std::string::npos) {
    p_2dhisto->Insert(0.,0.,0.,ncount);
    return;
  }
  p_2dhisto->Insert(Calc1(list[pos]),Calc2((*reflist)[refpos]),weight,ncount);
}

void STwo2D_Particle_Observable_Base::EvaluateNLOcontrib(double weight,double ncount)
{
  ATOOLS::Particle_List * plist =p_ana->GetParticleList(m_listname);
  ATOOLS::Particle_List *reflist=p_ana->GetParticleList(m_reflist);
  int no=-1, refno=-1; 
  size_t pos=std::string::npos, refpos=std::string::npos;
  for (size_t i=0;i<plist->size();++i) {
    if ((*plist)[i]->Flav()==m_flavour || 
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
  if (pos==std::string::npos || refpos==std::string::npos) {
    p_2dhisto->InsertMCB(0.,0.,0.,ncount);
    return;
  }
  p_2dhisto->InsertMCB(Calc1((*plist)[pos]),Calc2((*reflist)[refpos]),weight,ncount);
}

void STwo2D_Particle_Observable_Base::EvaluateNLOevt()
{
  p_2dhisto->FinishMCB();
}

Primitive_Observable_Base & STwo2D_Particle_Observable_Base::operator+=(const Primitive_Observable_Base & ob)
{
  STwo2D_Particle_Observable_Base* ob2d = (STwo2D_Particle_Observable_Base*)(&ob);
  if (p_2dhisto) {
    (*p_2dhisto)+=(*ob2d->p_2dhisto);
  }
  else {
    abort();
  }
  return *this;
}

void STwo2D_Particle_Observable_Base::EndEvaluation(double scale) {
  if (p_2dhisto) {
    p_2dhisto->Finalize();
    if (scale!=1.) p_2dhisto->Scale(scale);
    p_2dhisto->Output();
  }
}

void STwo2D_Particle_Observable_Base::Restore(double scale) 
{
  if (p_2dhisto) {
    if (scale!=1.) p_2dhisto->Scale(scale);
    p_2dhisto->Restore();
  }
}

void STwo2D_Particle_Observable_Base::Output(const std::string & pname) {
  if (p_2dhisto) {
    ATOOLS::MakeDir(pname); 
    p_2dhisto->Output((pname+std::string("/")+m_name).c_str());
  }
}

void STwo2D_Particle_Observable_Base::Reset()
{
  if (p_2dhisto) p_2dhisto->Reset();
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

double Two_DPhi_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().DPhi(p2->Momentum());
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

double Two_DEta_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return dabs(p1->Momentum().Eta()-p2->Momentum().Eta());
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

double Two_PEta_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().Eta()*p2->Momentum().Eta();
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

double Two_DY_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return dabs(p1->Momentum().Y()-p2->Momentum().Y());
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

double Two_PY_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().Y()*p2->Momentum().Y();
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

double Two_Mass_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return (p1->Momentum()+p2->Momentum()).Mass();
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

double Two_PT_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return (p1->Momentum()+p2->Momentum()).PPerp();
}

Primitive_Observable_Base *Two_PT_Distribution::Copy() const
{
  return new Two_PT_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_DPT_Distribution,
			     Two_DPT_Distribution_Getter,"TwoDPT")

  Two_DPT_Distribution::
Two_DPT_Distribution(const ATOOLS::Flavour flav,const size_t item,
		     const ATOOLS::Flavour refflav,const size_t refitem,
		     const int type,const double min,const double max,const int bins,
		     const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoDPT") {}

double Two_DPT_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().PPerp()-p2->Momentum().PPerp();
}

Primitive_Observable_Base *Two_DPT_Distribution::Copy() const
{
  return new Two_DPT_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				  m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_RPT_Distribution,
			     Two_RPT_Distribution_Getter,"TwoRPT")

  Two_RPT_Distribution::
Two_RPT_Distribution(const ATOOLS::Flavour flav,const size_t item,
		     const ATOOLS::Flavour refflav,const size_t refitem,
		     const int type,const double min,const double max,const int bins,
		     const std::string &inlist,const std::string &reflist):
  STwo_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				inlist,reflist,"TwoRPT") {}

double Two_RPT_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().PPerp()/p2->Momentum().PPerp();
}

Primitive_Observable_Base *Two_RPT_Distribution::Copy() const
{
  return new Two_RPT_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
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

double Two_DR_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return sqrt(sqr(p1->Momentum().Eta()-p2->Momentum().Eta())+
	      sqr(p1->Momentum().DPhi(p2->Momentum())));
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

double Two_ETFrac_Distribution::Calc(const Particle *p1,const Particle *p2)
{
  return p1->Momentum().EPerp()/p2->Momentum().EPerp();
}

Primitive_Observable_Base *Two_ETFrac_Distribution::Copy() const
{
  return new Two_ETFrac_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				     m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

DEFINE_TWO_OBSERVABLE_GETTER(Two_PT2D_Distribution,
			     Two_PT2D_Distribution_Getter,"TwoPT2D")
	  
  Two_PT2D_Distribution::
Two_PT2D_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const ATOOLS::Flavour refflav,const size_t refitem,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist,const std::string &reflist):
  STwo2D_Particle_Observable_Base(flav,item,refflav,refitem,type,min,max,bins,
				  inlist,reflist,"TwoPT2D") {}

double Two_PT2D_Distribution::Calc1(const Particle *p)
{
  return (p->Momentum()).PPerp();
}

double Two_PT2D_Distribution::Calc2(const Particle *p)
{
  return (p->Momentum()).PPerp();
}

Primitive_Observable_Base *Two_PT2D_Distribution::Copy() const
{
  return new Two_PT2D_Distribution(m_flavour,m_item,m_refflavour,m_refitem,
				 m_type,m_xmin,m_xmax,m_nbins,m_listname,m_reflist);
}

	  
