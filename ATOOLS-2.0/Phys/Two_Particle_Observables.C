#include "Two_Particle_Observables.H"

using namespace ATOOLS;
using namespace std;

Two_Particle_Observable_Base::Two_Particle_Observable_Base(Flavour & _flav1,Flavour & _flav2,
							   int _type,double _xmin,double _xmax,int _nbins,
							   std::string _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav1(_flav1), m_flav2(_flav2)
{
  name       = _name + m_flav1.Name() + m_flav2.Name() +std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Two_Particle_Observable_Base::Evaluate(double _value,double _weight) {
  histo->Insert(_value,_weight); 
}

 
void Two_Particle_Observable_Base::Evaluate(int _nout,ATOOLS::Vec4D * _moms,ATOOLS::Flavour * _flavs,
					    double _weight) 
{
  for (int i=0;i<_nout;i++) { 
    if (_flavs[i]==m_flav1) {
      for (int j=0;j<_nout;j++) { 
	if (_flavs[j]==m_flav2 && i!=j) Evaluate(_moms[i],_moms[j],_weight); 
      }
    }
  }
}


void Two_Particle_Observable_Base::Evaluate(const Particle_List & _plist,double _weight)
{
  size_t pos;
  bool   take;
  for (Particle_Const_Iterator plit1=_plist.begin();plit1!=_plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_Const_Iterator plit2=_plist.begin();plit2!=_plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),_weight);
	}
      }
    }
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Two_Particle_Mass::Two_Particle_Mass(Flavour & _flav1,Flavour & _flav2,
				     int _type,double _xmin,double _xmax,int _nbins,
				     std::string _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,_name) { }


void Two_Particle_Mass::Evaluate(Vec4D _mom1,Vec4D _mom2,double _weight) {
  double mass = sqrt((_mom1+_mom2).Abs2());
  histo->Insert(mass,_weight); 
} 

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Two_Particle_PT::Two_Particle_PT(Flavour & _flav1,Flavour & _flav2,
				 int _type,double _xmin,double _xmax,int _nbins,
				 std::string _name) :
  Two_Particle_Observable_Base(_flav1,_flav2,_type,_xmin,_xmax,_nbins,_name) { }


void Two_Particle_PT::Evaluate(Vec4D _mom1,Vec4D _mom2,double _weight) {
  double pt = sqrt(sqr(_mom1[1]+_mom2[1]) + sqr(_mom1[2]+_mom2[2]));
  histo->Insert(pt,_weight); 
} 

