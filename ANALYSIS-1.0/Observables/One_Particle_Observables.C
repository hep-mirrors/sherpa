#include "One_Particle_Observables.H"

using namespace ANALYSIS;
using namespace ATOOLS;
using namespace std;

One_Particle_Observable_Base::One_Particle_Observable_Base(const Flavour & _flav,
							   int _type,double _xmin,double _xmax,int _nbins,
							   const std::string & _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav(_flav)
{
  m_name     = _name + m_flav.Name()+std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void One_Particle_Observable_Base::Evaluate(double value,double weight, int ncount) 
{
  p_histo->Insert(value,weight,ncount); 
}

 
void One_Particle_Observable_Base::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, int ncount) 
{
  for (int i=0;i<nout;i++) { if (flavs[i]==m_flav) Evaluate(moms[i],weight,ncount); }
}


void One_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  for (Particle_Const_Iterator plit=plist.begin();plit!=plist.end();++plit) {
    if ((*plit)->Flav()==m_flav) {
      Evaluate((*plit)->Momentum(),weight, ncount);
      return;
    }
  }
  
  Evaluate(Vec4D(1.,0,0,1.),0, ncount);
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


One_Particle_ET::One_Particle_ET(const Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"ET") { }


void One_Particle_ET::Evaluate(const Vec4D & mom,double weight,int ncount) 
{
  double sintheta = sqrt(1.-sqr(mom[3])/(sqr(mom[1])+sqr(mom[2])));
  double et = mom[0]*sintheta;

  double pt2 = sqr(mom[1])+sqr(mom[2]);
  double p2  = sqr(mom[3])+pt2;
  double net = mom[0]*sqrt(pt2/p2);

  p_histo->Insert(net,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_ET::Copy() const
{
  return new One_Particle_ET(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_PT::One_Particle_PT(const Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"PT") { }


void One_Particle_PT::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double pt = sqrt(sqr(mom[1])+sqr(mom[2]));
  p_histo->Insert(pt,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_PT::Copy() const 
{
  return new One_Particle_PT(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_Eta::One_Particle_Eta(const Flavour & _flav,
				   int _type,double _xmin,double _xmax,int _nbins,
				   const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"Eta") { }


void One_Particle_Eta::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  //  double eta = -log(tan(sqrt(sqr(mom[1])+sqr(mom[2]))/(2.*mom[3])));;
  //  p_histo->Insert(eta,weight,ncount);

  double pt2=sqr(mom[1])+sqr(mom[2]);
  double pp =sqrt(pt2+sqr(mom[3]));
  double pz =dabs(mom[3]);
  double sn =mom[3]/pz;
  double value= sn*20.;
  if (pt2>1.e-10*pp*pp) {
    value = sn*0.5*log(sqr(pp+pz)/pt2);
  }
  p_histo->Insert(value,weight,ncount);
} 

Primitive_Observable_Base * One_Particle_Eta::Copy() const
{
  return new One_Particle_Eta(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_E::One_Particle_E(const Flavour & _flav,
			       int _type,double _xmin,double _xmax,int _nbins,
			       const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"E") { }


void One_Particle_E::Evaluate(const Vec4D & mom,double weight, int ncount) 
{
  double E = mom[0];
  p_histo->Insert(E,weight,ncount); 
} 

Primitive_Observable_Base * One_Particle_E::Copy() const
{
  return new One_Particle_E(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_EVis::One_Particle_EVis(const Flavour & _flav,
			       int _type,double _xmin,double _xmax,int _nbins,
			       const std::string & _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,"EVis") { }


void One_Particle_EVis::Evaluate(const Vec4D & mom,double weight, int ncount) 
{ } 

void One_Particle_EVis::Evaluate(int nout,const ATOOLS::Vec4D * moms,const ATOOLS::Flavour * flavs,
					    double weight, int ncount) 
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<nout;i++) {
    momsum += moms[i];
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}


void One_Particle_EVis::Evaluate(const Particle_List & plist,double weight,int ncount )
{
  ATOOLS::Vec4D momsum = Vec4D(0.,0.,0.,0.);
  for (Particle_Const_Iterator plit=plist.begin();plit!=plist.end();++plit) {
    momsum += (*plit)->Momentum();
  }
  p_histo->Insert(momsum.Abs(),weight,ncount); 
}
Primitive_Observable_Base * One_Particle_EVis::Copy() const
{
  return new One_Particle_EVis(m_flav,m_type,m_xmin,m_xmax,m_nbins,m_name);
}
