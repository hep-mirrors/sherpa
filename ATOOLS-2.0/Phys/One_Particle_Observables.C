#include "One_Particle_Observables.H"

using namespace ATOOLS;
using namespace std;

One_Particle_Observable_Base::One_Particle_Observable_Base(Flavour & _flav,
							   int _type,double _xmin,double _xmax,int _nbins,
							   std::string _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL), 
  m_flav(_flav)
{
  name       = _name + std::string(m_flav.Name())+std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void One_Particle_Observable_Base::Evaluate(double _value,double _weight) {
  histo->Insert(_value,_weight); 
}

 
void One_Particle_Observable_Base::Evaluate(int _nout,ATOOLS::Vec4D * _moms,ATOOLS::Flavour * _flavs,
					    double _weight) 
{
  for (int i=0;i<_nout;i++) { if (_flavs[i]==m_flav) Evaluate(_moms[i],_weight); }
}


void One_Particle_Observable_Base::Evaluate(const Particle_List & _plist,double _weight)
{
  size_t pos;
  bool   take;
  for (Particle_Const_Iterator plit=_plist.begin();plit!=_plist.end();++plit) {
    if ((*plit)->Flav()==m_flav) Evaluate((*plit)->Momentum(),_weight);
  }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


One_Particle_ET::One_Particle_ET(Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 std::string _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,_name) { }


void One_Particle_ET::Evaluate(Vec4D _mom,double _weight) {
  double sintheta = sqrt(1.-sqr(_mom[3])/(sqr(_mom[1])+sqr(_mom[2])));
  histo->Insert(_mom[0]*sintheta,_weight); 
} 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_PT::One_Particle_PT(Flavour & _flav,
				 int _type,double _xmin,double _xmax,int _nbins,
				 std::string _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,_name) { }


void One_Particle_PT::Evaluate(Vec4D _mom,double _weight) {
  double pt = sqrt(sqr(_mom[1])+sqr(_mom[2]));
  histo->Insert(pt,_weight); 
} 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_Eta::One_Particle_Eta(Flavour & _flav,
				   int _type,double _xmin,double _xmax,int _nbins,
				   std::string _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,_name) { }


void One_Particle_Eta::Evaluate(Vec4D _mom,double _weight) {
  double eta = -log(tan(sqrt(sqr(_mom[1])+sqr(_mom[2]))/(2.*_mom[3])));;
  histo->Insert(eta,_weight); 
} 


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

One_Particle_E::One_Particle_E(Flavour & _flav,
			       int _type,double _xmin,double _xmax,int _nbins,
			       std::string _name) :
  One_Particle_Observable_Base(_flav,_type,_xmin,_xmax,_nbins,_name) { }


void One_Particle_E::Evaluate(Vec4D _mom,double _weight) {
  double E = _mom[0];
  histo->Insert(E,_weight); 
} 
