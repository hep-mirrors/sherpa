#include "Scaled_Observables.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

Scaled_Observable_Base::Scaled_Observable_Base(int type,double xmin,double xmax,int nbins,
					       const std::string & listname, const std::string & name,
					       double ecms) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL), m_ecms(ecms)
{
  MyStrStream str;
  str<<name<<".dat";
  str>>m_name;

  if (listname!=std::string("")) m_listname = listname;
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Scaled_Observable_Base::Evaluate(double value,double weight, int ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void Scaled_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,
					    double weight, int ncount) 
{
  for (int i=0;i<nout;i++) Evaluate(moms[i],weight,ncount);
}


void Scaled_Observable_Base::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  for (Particle_Const_Iterator plit=plist.begin();plit!=plist.end();++plit) {
    Evaluate((*plit)->Momentum(),weight, ncount);
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Scaled_Momentum::Scaled_Momentum(int type,double xmin,double xmax,int nbins,
				 const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"ScaledMomentum",ecms) { }


void Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;

  p_histo->Insert(xp,weight,ncount); 
} 

Primitive_Observable_Base * Scaled_Momentum::Copy() const
{
  return new Scaled_Momentum(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Log_Scaled_Momentum::Log_Scaled_Momentum(int type,double xmin,double xmax,int nbins,
				 const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"LogScaledMomentum", ecms) { }


void Log_Scaled_Momentum::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double xp = 2.*Vec3D(mom).Abs()/m_ecms;
  double xi = - log(xp);

  p_histo->Insert(xi,weight,ncount); 
} 

Primitive_Observable_Base * Log_Scaled_Momentum::Copy() const
{
  return new Log_Scaled_Momentum(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Scaled_Energy::Scaled_Energy(int type,double xmin,double xmax,int nbins,
			     const std::string & listname, double ecms) :
  Scaled_Observable_Base(type,xmin,xmax,nbins,listname,"ScaledEnergy",ecms) { }


void Scaled_Energy::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double E = 2.*mom[0]/m_ecms;
  p_histo->Insert(E,weight,ncount); 
} 

Primitive_Observable_Base * Scaled_Energy::Copy() const
{
  return new Scaled_Energy(m_type,m_xmin,m_xmax,m_nbins,m_listname,m_ecms);
}

