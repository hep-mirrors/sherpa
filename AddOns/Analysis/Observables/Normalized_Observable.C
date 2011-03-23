#include "AddOns/Analysis/Observables/Normalized_Observable.H"

#include "ATOOLS/Org/Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Normalized_Observable::Normalized_Observable():
  Primitive_Observable_Base()
{
}

Normalized_Observable::
Normalized_Observable(int type,double xmin,double xmax,int nbins,
		      const std::string &name,const int mode):
  Primitive_Observable_Base(type,xmin,xmax,nbins), m_mode(mode)
{
  if (!(type&4)) THROW(fatal_error,"Must be initialized with type X+4");
  p_obs = new Histogram(m_type,m_xmin,m_xmax,m_nbins,name);
  p_norm = new Histogram(m_type,m_xmin,m_xmax,m_nbins,name);
}

Normalized_Observable::
Normalized_Observable(const Normalized_Observable & old):
  Primitive_Observable_Base(old)
{
  if (old.p_histo) {
    p_obs = new Histogram(m_type,m_xmin,m_xmax,m_nbins,old.Name());
    p_norm = new Histogram(m_type,m_xmin,m_xmax,m_nbins,old.Name());
  }
}

Normalized_Observable::~Normalized_Observable()
{
  if (p_histo) {
    delete p_obs;
    delete p_norm;
  }
}

void Normalized_Observable::Reset()
{
  p_obs->Reset();
  p_norm->Reset();
}

Primitive_Observable_Base &
Normalized_Observable::operator+=(const Primitive_Observable_Base &obs)
{
  const Normalized_Observable *nobs((const Normalized_Observable*)&obs);
  if (p_histo) {
    (*p_obs)+=(*nobs->p_obs);
    (*p_norm)+=(*nobs->p_norm);
    if (!m_copied) {                                                 
      double n=ATOOLS::Max(1.0,double(p_obs->Fills()));                 
      for (int i=0;i<m_nbins+2;++i) {			   
	double val=p_obs->Bin(i)==0.0?0.0:
	  p_obs->Bin(i)/p_norm->Bin(i);
	p_histo->SetBin(i,val);
	double w2mo(n*p_obs->Bin2(i)), wmo(p_obs->Bin(i));
	double sigo(sqrt((w2mo-wmo*wmo)/Max(1.0,n-1.0)));
	double w2mi(n*p_norm->Bin2(i)), wmi(p_norm->Bin(i));
	double sigi(sqrt((w2mi-wmi*wmi)/Max(1.0,n-1.0)));
	double relerr=p_obs->Bin(i)==0.0?0.0:
	  (sigo/dabs(p_obs->Bin(i))+sigi/dabs(p_norm->Bin(i)));
	p_histo->SetBin2(i,sqr(val*relerr));
        p_histo->SetBinPS(i,p_norm->Bin(i));
	p_histo->SetBinPS2(i,sqr(sigi));
      }								   
      p_histo->SetFills((long int)p_obs->Fills());            
    }                                                                
  }
  else {
    THROW(critical_error,Name()+" has not overloaded the operator");
  }
  return *this;
}

void Normalized_Observable::EndEvaluation(double scale)
{
  double n=ATOOLS::Max(1.0,double(p_obs->Fills()));                 
  if (m_mode==0) p_obs->Scale(scale*m_nbins/(m_xmax-m_xmin)/n);
  else p_obs->Scale(scale/n);
  p_norm->Scale(scale/n);                                 
  if (!m_copied) {                                                 
    for (int i=0;i<m_nbins+2;++i) {			   
      double val=p_obs->Bin(i)==0.0?0.0:
	p_obs->Bin(i)/p_norm->Bin(i);
      p_histo->SetBin(i,val);
      double w2mo(n*p_obs->Bin2(i)), wmo(p_obs->Bin(i));
      double sigo(sqrt((w2mo-wmo*wmo)/Max(1.0,n-1.0)));
      double w2mi(n*p_norm->Bin2(i)), wmi(p_norm->Bin(i));
      double sigi(sqrt((w2mi-wmi*wmi)/Max(1.0,n-1.0)));
      double relerr=p_obs->Bin(i)==0.0?0.0:
	(sigo/dabs(p_obs->Bin(i))+sigi/dabs(p_norm->Bin(i)));
      p_histo->SetBin2(i,sqr(val*relerr));
      p_histo->SetBinPS(i,p_norm->Bin(i));
      p_histo->SetBinPS2(i,sqr(sigi));
    }								   
    p_histo->SetFills((long int)p_obs->Fills());            
  }                                                                
}

void Normalized_Observable::Restore(double scale)
{
  double n=ATOOLS::Max(1.0,double(p_obs->Fills()));                 
  if (m_mode==0) p_obs->Scale(scale*n*(m_xmax-m_xmin)/m_nbins);
  else p_obs->Scale(scale*n);
  p_norm->Scale(scale*n);                                 
}

void Normalized_Observable::Fill
(const double &x,const double &y,
 const double &weight,const double &ntrial)
{
  p_obs->Insert(x,y*weight,ntrial);
  p_norm->Insert(x,weight,ntrial);
}

void Normalized_Observable::FillMCB
(const double &x,const double &y,
 const double &weight,const double &ntrial)
{
  p_obs->InsertMCB(x,y*weight,ntrial);
  p_norm->InsertMCB(x,weight,ntrial);
}

void Normalized_Observable::FinishMCB()
{
  p_obs->FinishMCB();
  p_norm->FinishMCB();
}

Primitive_Observable_Base *Normalized_Observable::Copy() const
{
  return new Normalized_Observable(m_type,m_xmin,m_xmax,m_nbins,Name());
}
