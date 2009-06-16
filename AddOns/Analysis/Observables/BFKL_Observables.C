#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class Decorrelation_vs_DEta: public Normalized_Observable {  
  protected:

    double m_ptmin;
    
  public:

    Decorrelation_vs_DEta(const int type,const double &ptmin,
			  const double etamin,const double etamax,
			  const int nbins,const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Decorrelation_vs_DEta

  class BFKL_DEta: public Primitive_Observable_Base {  
  protected:

    double m_ptmin;
    
  public:

    BFKL_DEta(const int type,const double &ptmin,
	 const double etamin,const double etamax,
	 const int nbins,const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class BFKL_DEta

  class BFKL_Two_Eta: public Primitive_Observable_Base {  
  protected:

    double m_ptmin;
    
  public:

    BFKL_Two_Eta(const int type,const double &ptmin,
	 const double etamin,const double etamax,
	 const int nbins,const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class BFKL_Two_Eta

  class BFKL_DPhi: public Primitive_Observable_Base {  
  protected:

    double m_ptmin, m_detamin, m_detamax;
    
  public:

    BFKL_DPhi(const int type,const double &ptmin,
	      const double detamin,const double detamax,
	      const double phimin,const double phimax,
	      const int nbins,const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,double ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class BFKL_DPhi

}// end of namespace ANALYSIS

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include <algorithm>

DECLARE_GETTER(Decorrelation_vs_DEta_Getter,"CosDPhivsDEta",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *Decorrelation_vs_DEta_Getter::
operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:"FinalState";
    return new Decorrelation_vs_DEta
      (HistogramType(parameters[0][4]),
       ATOOLS::ToType<double>(parameters[0][0]),
       ATOOLS::ToType<double>(parameters[0][1]),
       ATOOLS::ToType<double>(parameters[0][2]),
       ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  return NULL;
}									

void Decorrelation_vs_DEta_Getter::
PrintInfo(std::ostream &str,const size_t width) const 
{ 
  str<<"etmin min max bins Lin|LinErr|Log|LogErr [list]"; 
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Poincare.H"

#include <fstream>

using namespace ATOOLS;

double DPhi(const Vec4D &a,Vec4D b)
{
  Poincare rot(Vec4D(0.0,a[1],a[2],0.0),Vec4D::XVEC);
  rot.Rotate(b);
  double dp(b.Phi());
  return dp<0.0?dp+2.0*M_PI:dp;
}

class Order_Eta {
public:
  bool operator()(const Particle *a,const Particle *b)
  {
    return a->Momentum().Eta()>b->Momentum().Eta();
  }
};

Decorrelation_vs_DEta::
Decorrelation_vs_DEta(const int type,const double &ptmin,
		      const double etamin,const double etamax,
		      const int nbins,const std::string &listname):
  Normalized_Observable(type,etamin,etamax,nbins)
{
  m_listname=listname;
  m_ptmin=ptmin;
  m_name="CosDPhi_vs_DEta_"+m_listname+".dat";
}
    
void Decorrelation_vs_DEta::Evaluate(const ATOOLS::Particle_List &plist,
				     double weight,double ncount)
{
  Particle_List list(plist);
  std::stable_sort(list.begin(),list.end(),Order_Eta());
  if (list.size()<2 || (list.front()->Momentum().EPerp()<m_ptmin &&
			list.back()->Momentum().EPerp()<m_ptmin)) {
    p_obs->Insert(0.0,0.0,ncount);
    p_norm->Insert(0.0,0.0,ncount);
    return;
  }
  double dphijet(DPhi(list.front()->Momentum(),list.back()->Momentum()));
  double detajet(dabs(list.front()->Momentum().
		      DEta(list.back()->Momentum())));
  p_obs->Insert(detajet,cos(M_PI-dphijet)*weight,ncount);
  p_norm->Insert(detajet,weight,ncount);
}
    
Primitive_Observable_Base *Decorrelation_vs_DEta::Copy() const
{
  return new 
    Decorrelation_vs_DEta(m_type,m_ptmin,m_xmin,m_xmax,m_nbins,m_listname);
}

DECLARE_GETTER(BFKL_DEta_Getter,"BFKLDEta",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *BFKL_DEta_Getter::
operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:"FinalState";
    return new BFKL_DEta
      (HistogramType(parameters[0][4]),
       ATOOLS::ToType<double>(parameters[0][0]),
       ATOOLS::ToType<double>(parameters[0][1]),
       ATOOLS::ToType<double>(parameters[0][2]),
       ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  return NULL;
}									

void BFKL_DEta_Getter::
PrintInfo(std::ostream &str,const size_t width) const 
{ 
  str<<"etmin min max bins Lin|LinErr|Log|LogErr [list]"; 
}

BFKL_DEta::BFKL_DEta(const int type,const double &ptmin,
		     const double etamin,const double etamax,
		     const int nbins,const std::string &listname):
  Primitive_Observable_Base(type,etamin,etamax,nbins)
{
  m_listname=listname;
  m_ptmin=ptmin;
  m_name="BFKL_DEta_"+m_listname+".dat";
}
    
void BFKL_DEta::Evaluate(const ATOOLS::Particle_List &plist,
				     double weight,double ncount)
{
  Particle_List list(plist);
  std::stable_sort(list.begin(),list.end(),Order_Eta());
  if (list.size()<2 || (list.front()->Momentum().EPerp()<m_ptmin &&
			list.back()->Momentum().EPerp()<m_ptmin)) {
    p_histo->Insert(0.0,0.0,ncount);
    return;
  }
  double detajet(dabs(list.front()->Momentum().DEta(list.back()->Momentum())));
  p_histo->Insert(detajet,weight,ncount);
}
    
Primitive_Observable_Base *BFKL_DEta::Copy() const
{
  return new BFKL_DEta(m_type,m_ptmin,m_xmin,m_xmax,m_nbins,m_listname);
}

DECLARE_GETTER(BFKL_Two_Eta_Getter,"BFKLTwoEta",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *BFKL_Two_Eta_Getter::
operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<5) return NULL;
    std::string list=parameters[0].size()>5?parameters[0][5]:"FinalState";
    return new BFKL_Two_Eta
      (HistogramType(parameters[0][4]),
       ATOOLS::ToType<double>(parameters[0][0]),
       ATOOLS::ToType<double>(parameters[0][1]),
       ATOOLS::ToType<double>(parameters[0][2]),
       ATOOLS::ToType<int>(parameters[0][3]),list);
  }
  return NULL;
}									

void BFKL_Two_Eta_Getter::
PrintInfo(std::ostream &str,const size_t width) const 
{ 
  str<<"etmin min max bins Lin|LinErr|Log|LogErr [list]"; 
}

BFKL_Two_Eta::BFKL_Two_Eta(const int type,const double &ptmin,
		     const double etamin,const double etamax,
		     const int nbins,const std::string &listname):
  Primitive_Observable_Base(type,etamin,etamax,nbins)
{
  m_listname=listname;
  m_ptmin=ptmin;
  m_name="BFKL_Two_Eta_"+m_listname+".dat";
}
    
void BFKL_Two_Eta::Evaluate(const ATOOLS::Particle_List &plist,
				     double weight,double ncount)
{
  Particle_List list(plist);
  std::stable_sort(list.begin(),list.end(),Order_Eta());
  if (list.size()<2 || (list.front()->Momentum().EPerp()<m_ptmin &&
			list.back()->Momentum().EPerp()<m_ptmin)) {
    p_histo->Insert(0.0,0.0,ncount);
    return;
  }
  double detajet((list.front()->Momentum().Eta()+
		  list.back()->Momentum().Eta())/2.0);
  p_histo->Insert(detajet,weight,ncount);
}
    
Primitive_Observable_Base *BFKL_Two_Eta::Copy() const
{
  return new BFKL_Two_Eta(m_type,m_ptmin,m_xmin,m_xmax,m_nbins,m_listname);
}

DECLARE_GETTER(BFKL_DPhi_Getter,"BFKLDPhi",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *BFKL_DPhi_Getter::
operator()(const Argument_Matrix &parameters) const
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    std::string list=parameters[0].size()>7?parameters[0][7]:"FinalState";
    return new BFKL_DPhi
      (HistogramType(parameters[0][6]),
       ATOOLS::ToType<double>(parameters[0][0]),
       ATOOLS::ToType<double>(parameters[0][1]),
       ATOOLS::ToType<double>(parameters[0][2]),
       ATOOLS::ToType<double>(parameters[0][3]),
       ATOOLS::ToType<double>(parameters[0][4]),
       ATOOLS::ToType<int>(parameters[0][5]),list);
  }
  return NULL;
}									

void BFKL_DPhi_Getter::
PrintInfo(std::ostream &str,const size_t width) const 
{ 
  str<<"etmin detamin detamax min max bins Lin|LinErr|Log|LogErr [list]"; 
}

BFKL_DPhi::BFKL_DPhi(const int type,const double &ptmin,
		     const double detamin,const double detamax,
		     const double phimin,const double phimax,
		     const int nbins,const std::string &listname):
  Primitive_Observable_Base(type,phimin,phimax,nbins)
{
  m_listname=listname;
  m_ptmin=ptmin;
  m_detamin=detamin;
  m_detamax=detamax;
  m_name="BFKL_DPhi_"+ToString(m_detamin)+"-"+
    ToString(m_detamax)+"_"+m_listname+".dat";
}
    
void BFKL_DPhi::Evaluate(const ATOOLS::Particle_List &plist,
				     double weight,double ncount)
{
  Particle_List list(plist);
  std::stable_sort(list.begin(),list.end(),Order_Eta());
  double detajet(list.size()<2?0.0:
		 dabs(list.front()->Momentum().DEta(list.back()->Momentum())));
  if (list.size()<2 || detajet<m_detamin || detajet>m_detamax ||
      (list.front()->Momentum().EPerp()<m_ptmin &&
       list.back()->Momentum().EPerp()<m_ptmin)) {
    p_histo->Insert(0.0,0.0,ncount);
    return;
  }
  double dphijet(DPhi(list.front()->Momentum(),list.back()->Momentum()));
  p_histo->Insert(1.0-dphijet/M_PI,weight,ncount);
}
    
Primitive_Observable_Base *BFKL_DPhi::Copy() const
{
  return new BFKL_DPhi(m_type,m_ptmin,m_detamin,m_detamax,
		       m_xmin,m_xmax,m_nbins,m_listname);
}

