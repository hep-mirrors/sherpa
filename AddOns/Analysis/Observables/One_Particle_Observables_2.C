#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ATOOLS;

namespace ANALYSIS {

  class SOne_Particle_Observable_Base: public Primitive_Observable_Base {  
  protected:

    ATOOLS::Flavour m_flavour;
    size_t          m_item;

  public:

    SOne_Particle_Observable_Base
    (const ATOOLS::Flavour flav,const size_t item,
     const int type,const double min,const double max,const int bins,
     const std::string &inlist,const std::string &name);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    virtual bool Evaluate(const Particle &p1,
			  double weight=1.,double ncount=1) const = 0;

  };// end of class SOne_Particle_Observable_Base

  class One_Phi_Distribution: public SOne_Particle_Observable_Base {  
  public:

    One_Phi_Distribution(const ATOOLS::Flavour flav,const size_t item,
			 const int type,const double min,const double max,const int bins,
			 const std::string &inlist);
    
    bool Evaluate(const Particle &p1,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class One_Phi_Distribution

  class One_Eta_Distribution: public SOne_Particle_Observable_Base {  
  public:

    One_Eta_Distribution(const ATOOLS::Flavour flav,const size_t item,
			 const int type,const double min,const double max,const int bins,
			 const std::string &inlist);
    
    bool Evaluate(const Particle &p1,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class One_Eta_Distribution

  class One_Y_Distribution: public SOne_Particle_Observable_Base {  
  public:

    One_Y_Distribution(const ATOOLS::Flavour flav,const size_t item,
		       const int type,const double min,const double max,const int bins,
		       const std::string &inlist);
    
    bool Evaluate(const Particle &p1,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class One_Y_Distribution

  class One_Mass_Distribution: public SOne_Particle_Observable_Base {  
  public:

    One_Mass_Distribution(const ATOOLS::Flavour flav,const size_t item,
			  const int type,const double min,const double max,const int bins,
			  const std::string &inlist);
    
    bool Evaluate(const Particle &p1,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class One_Mass_Distribution

  class One_PT_Distribution: public SOne_Particle_Observable_Base {  
  public:

    One_PT_Distribution(const ATOOLS::Flavour flav,const size_t item,
			const int type,const double min,const double max,const int bins,
			const std::string &inlist);
    
    bool Evaluate(const Particle &p1,
		  double weight=1.,double ncount=1) const;

    Primitive_Observable_Base *Copy() const;
    
  };// end of class One_PT_Distribution

}// end of namespace ANALYSIS

using namespace ANALYSIS;

template <class Class>
Primitive_Observable_Base *
GetSOneParticleObservable(const Argument_Matrix &parameters) 
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<6) return NULL;
    std::string list(parameters[0].size()>6?parameters[0][6]:"FinalState");
    int kf=ATOOLS::ToType<int>(parameters[0][0]);
    ATOOLS::Flavour flav((kf_code)abs(kf));
    if (kf<0) flav=flav.Bar();
    return new Class(flav,ATOOLS::ToType<size_t>(parameters[0][1]),
		     HistogramType(parameters[0][5]),
		     ATOOLS::ToType<double>(parameters[0][2]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<int>(parameters[0][4]),list);
  }
  return NULL;
}									

#define DEFINE_ONE_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetSOneParticleObservable<CLASS>(parameters); }

#define DEFINE_ONE_OBSERVABLE_PRINT_METHOD(NAME)		\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"flav item min max bins Lin|LinErr|Log|LogErr [inlist]"; }

#define DEFINE_ONE_OBSERVABLE_GETTER(CLASS,NAME,TAG)		\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_ONE_OBSERVABLE_GETTER_METHOD(CLASS,NAME)		\
  DEFINE_ONE_OBSERVABLE_PRINT_METHOD(NAME)

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

SOne_Particle_Observable_Base::
SOne_Particle_Observable_Base(const ATOOLS::Flavour flav,const size_t item,
			      const int type,const double min,const double max,const int bins,
			      const std::string &inlist,const std::string &name):
  Primitive_Observable_Base(type,min,max,bins), 
  m_flavour(flav),
  m_item(item)
{
  m_listname=inlist;
  m_name=name+"_"+ToString(m_flavour)+"-"+ToString(m_item)+".dat";
}

void SOne_Particle_Observable_Base::Evaluate(const ATOOLS::Particle_List &inlist,
					     double weight,double ncount)
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
  if (pos==std::string::npos) {
    p_histo->Insert(0,0,ncount);
    return;
  }
  Evaluate(*inlist[pos],weight,ncount);
}

DEFINE_ONE_OBSERVABLE_GETTER(One_Phi_Distribution,
			     One_Phi_Distribution_Getter,"OnePhi")

  One_Phi_Distribution::
One_Phi_Distribution(const ATOOLS::Flavour flav,const size_t item,
		     const int type,const double min,const double max,const int bins,
		     const std::string &inlist):
  SOne_Particle_Observable_Base(flav,item,type,min,max,bins,
				inlist,"OnePhi") {}

bool One_Phi_Distribution::Evaluate(const Particle &p1,
				    double weight,double ncount) const
{
  p_histo->Insert(p1.Momentum().Phi(),weight,ncount);
  return true;
}

Primitive_Observable_Base *One_Phi_Distribution::Copy() const
{
  return new One_Phi_Distribution(m_flavour,m_item,m_type,
				  m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_ONE_OBSERVABLE_GETTER(One_Eta_Distribution,
			     One_Eta_Distribution_Getter,"OneEta")

  One_Eta_Distribution::
One_Eta_Distribution(const ATOOLS::Flavour flav,const size_t item,
		     const int type,const double min,const double max,const int bins,
		     const std::string &inlist):
  SOne_Particle_Observable_Base(flav,item,type,min,max,bins,
				inlist,"OneEta") {}

bool One_Eta_Distribution::Evaluate(const Particle &p1,
				    double weight,double ncount) const
{
  p_histo->Insert(p1.Momentum().Eta(),weight,ncount);
  return true;
}

Primitive_Observable_Base *One_Eta_Distribution::Copy() const
{
  return new One_Eta_Distribution(m_flavour,m_item,m_type,
				  m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_ONE_OBSERVABLE_GETTER(One_Y_Distribution,
			     One_Y_Distribution_Getter,"OneY")

  One_Y_Distribution::
One_Y_Distribution(const ATOOLS::Flavour flav,const size_t item,
		   const int type,const double min,const double max,const int bins,
		   const std::string &inlist):
  SOne_Particle_Observable_Base(flav,item,type,min,max,bins,
				inlist,"OneY") {}

bool One_Y_Distribution::Evaluate(const Particle &p1,
				  double weight,double ncount) const
{
  p_histo->Insert(p1.Momentum().Y(),weight,ncount);
  return true;
}

Primitive_Observable_Base *One_Y_Distribution::Copy() const
{
  return new One_Y_Distribution(m_flavour,m_item,m_type,
				m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_ONE_OBSERVABLE_GETTER(One_Mass_Distribution,
			     One_Mass_Distribution_Getter,"OneMass")

  One_Mass_Distribution::
One_Mass_Distribution(const ATOOLS::Flavour flav,const size_t item,
		      const int type,const double min,const double max,const int bins,
		      const std::string &inlist):
  SOne_Particle_Observable_Base(flav,item,type,min,max,bins,
				inlist,"OneMass") {}

bool One_Mass_Distribution::Evaluate(const Particle &p1,
				     double weight,double ncount) const
{
  p_histo->Insert(p1.Momentum().Mass(),weight,ncount);
  return true;
}

Primitive_Observable_Base *One_Mass_Distribution::Copy() const
{
  return new One_Mass_Distribution(m_flavour,m_item,m_type,
				   m_xmin,m_xmax,m_nbins,m_listname);
}

DEFINE_ONE_OBSERVABLE_GETTER(One_PT_Distribution,
			     One_PT_Distribution_Getter,"OnePT")

  One_PT_Distribution::
One_PT_Distribution(const ATOOLS::Flavour flav,const size_t item,
		    const int type,const double min,const double max,const int bins,
		    const std::string &inlist):
  SOne_Particle_Observable_Base(flav,item,type,min,max,bins,
				inlist,"OnePT") {}

bool One_PT_Distribution::Evaluate(const Particle &p1,
				   double weight,double ncount) const
{
  p_histo->Insert(p1.Momentum().PPerp(),weight,ncount);
  return true;
}

Primitive_Observable_Base *One_PT_Distribution::Copy() const
{
  return new One_PT_Distribution(m_flavour,m_item,m_type,
				 m_xmin,m_xmax,m_nbins,m_listname);
}

