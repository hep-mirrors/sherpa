#include "Two_Particle_Observables.H"
#include "MyStrStream.H"

#ifdef ROOT_SUPPORT
#include "Scaling.H"
#include "My_Root.H"
#include "TH2D.h"
#endif 

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

Two_Particle_Observable_Base::Two_Particle_Observable_Base(const Flavour & _flav1,const Flavour & _flav2,
							   int _type,double _xmin,double _xmax,int _nbins,
							   const std::string & _name,const std::string & _lname) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav1(_flav1), m_flav2(_flav2)
{
  m_listname=_lname;
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
				     const std::string & _name, const std::string & _lname) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Mass",_lname) { }


void Two_Particle_Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double mass = sqrt((mom1+mom2).Abs2());
  p_histo->Insert(mass,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Mass::Copy() const
{
  return new Two_Particle_Mass(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Two_Particle_PT::Two_Particle_PT(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name, const std::string & _lname) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"PT",_lname) { }


void Two_Particle_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]) + sqr(mom1[2]+mom2[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_PT::Copy() const 
{
  return new Two_Particle_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

Two_Particle_Scalar_PT::Two_Particle_Scalar_PT(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name, const std::string & _lname) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"SPT",_lname) { }

void Two_Particle_Scalar_PT::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  p_histo->Insert((mom1.PPerp()+mom2.PPerp())/2.,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Scalar_PT::Copy() const 
{
  return new Two_Particle_Scalar_PT(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

#ifdef ROOT_SUPPORT
Two_Particle_Angles::Two_Particle_Angles(const Flavour & _flav1,const Flavour & _flav2,
					 int _type,double _xmin,double _xmax,int _nbins,
					 const std::string & _name, const std::string & _lname) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Angles",_lname) 
{ 
  (*MYROOT::myroot)(new TH2D(ATOOLS::ToString(this).c_str(),
			     (m_flav1.Name()+std::string("_")+m_flav2.Name()+
			      std::string("_Angles")).c_str(),
			     64,0.,M_PI,64,0.,M_PI),
		    m_flav1.Name()+std::string("_")+m_flav2.Name()+
		    std::string("_Angles"));
}

void Two_Particle_Angles::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  ((TH2D*)(*MYROOT::myroot)[m_flav1.Name()+std::string("_")+m_flav2.Name()+
			    std::string("_Angles")])->Fill(mom1.Theta(),mom2.Theta(),weight);
} 

Primitive_Observable_Base * Two_Particle_Angles::Copy() const 
{
  return new Two_Particle_Angles(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}
#endif

Two_Particle_Eta::Two_Particle_Eta(const Flavour & _flav1,const Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name,const std::string & _lname) :
    Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"Eta",_lname) 

{
}


void Two_Particle_Eta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{
  double eta = (mom1+mom2).Eta();
  p_histo->Insert(eta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_Eta::Copy() const 
{
  return new Two_Particle_Eta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Two_Particle_DEta::Two_Particle_DEta(const Flavour & _flav1,const Flavour & _flav2,
				     int _type,double _xmin,double _xmax,int _nbins,
				     const std::string & _name,const std::string & _lname) :
    Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"deta",_lname) 

{ 
}


void Two_Particle_DEta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{    
    double deta = abs((mom1.Eta()-mom2.Eta()));
    p_histo->Insert(deta,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DEta::Copy() const 
{
    return new Two_Particle_DEta(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}
    
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Two_Particle_DPhi::Two_Particle_DPhi(const Flavour & _flav1,const Flavour & _flav2,
				       int _type,double _xmin,double _xmax,int _nbins,
				     const std::string & _name,const std::string & _lname) :
    Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"dphi",_lname) 
{ 
}


void Two_Particle_DPhi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    p_histo->Insert(dphi,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DPhi::Copy() const 
{
    return new Two_Particle_DPhi(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Two_Particle_DR::Two_Particle_DR(const Flavour & _flav1,const Flavour & _flav2,
				       int _type,double _xmin,double _xmax,int _nbins,
				     const std::string & _name,const std::string & _lname) :
    Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,"dr",_lname) 

{
 }


void Two_Particle_DR::Evaluate(const Vec4D & mom1,const Vec4D & mom2,double weight, int ncount) 
{ 
    double pt1=sqrt(mom1[1]*mom1[1]+mom1[2]*mom1[2]);
    double pt2=sqrt(mom2[1]*mom2[1]+mom2[2]*mom2[2]);
    double dphi=acos((mom1[1]*mom2[1]+mom1[2]*mom2[2])/(pt1*pt2));
    double c1=mom1[3]/Vec3D(mom1).Abs();
    double c2=mom2[3]/Vec3D(mom2).Abs();
    double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
    double dr= sqrt(sqr(deta) + sqr(dphi)); 
    cout<<"Deat in DR "<<deta<<" DR is :  "<<dr<<endl;
    p_histo->Insert(dr,weight,ncount); 
} 

Primitive_Observable_Base * Two_Particle_DR::Copy() const 
{
    return new Two_Particle_DR(m_flav1,m_flav2,m_type,m_xmin,m_xmax,m_nbins,m_name,m_listname);
}

