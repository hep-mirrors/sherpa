#include "Two_Particle_Observables.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

Two_Particle_Observable_Base::Two_Particle_Observable_Base(const Flavour & _flav1,const Flavour & _flav2,
							   int _type,double _xmin,double _xmax,int _nbins,
							   const std::string & _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav1(_flav1), m_flav2(_flav2)
{
  MyStrStream str;
  str<<_name<<m_flav1<<m_flav2<<".dat";
  str>>m_name;
  //  m_name     = _name + m_flav1.Name() + m_flav2.Name() +std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Two_Particle_Observable_Base::Evaluate(double value,double weight, int ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void Two_Particle_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms, const ATOOLS::Flavour * flavs,
					    double weight, int ncount) 
{
  for (int i=0;i<nout;i++) { 
    if (flavs[i]==m_flav1) {
      for (int j=0;j<nout;j++) { 
	if (flavs[j]==m_flav2 && i!=j) Evaluate(moms[i],moms[j],weight,ncount); 
      }
    }
  }
}


void Two_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight, int ncount)
{
  for (Particle_Const_Iterator plit1=plist.begin();plit1!=plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_Const_Iterator plit2=plist.begin();plit2!=plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),weight,ncount);
	  return;
	}
      }
    }
  }
  Evaluate(Vec4D(1.,0,0,1.),Vec4D(1.,0,0,-1.),0, ncount);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Two_Particle_Mass::Two_Particle_Mass(const Flavour & _flav1,const Flavour & _flav2,
				     int _type,double _xmin,double _xmax,int _nbins,
				     const std::string & _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Mass") { }


void Two_Particle_Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double mass = sqrt((mom1+mom2).Abs2());
  p_histo->Insert(mass,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Mass::Copy() const
{
  return new Two_Particle_Mass(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Two_Particle_PT::Two_Particle_PT(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"PT") { }


void Two_Particle_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_PT::Copy() const 
{
  return new Two_Particle_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name);
}


Two_Particle_Scalar_PT::Two_Particle_Scalar_PT(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"SPT") { }


void Two_Particle_Scalar_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  p_histo->Insert((mom1.PPerp()+mom2.PPerp())/2.,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Scalar_PT::Copy() const 
{
  return new Two_Particle_Scalar_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

Two_Particle_Eta::Two_Particle_Eta(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Eta") { }


void Two_Particle_Eta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double eta = (mom1+mom2).Eta();
  p_histo->Insert(eta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Eta::Copy() const 
{
  return new Two_Particle_Eta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name);
}
