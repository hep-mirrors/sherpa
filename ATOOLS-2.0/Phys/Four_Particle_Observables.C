#include "Four_Particle_Observables.H"

using namespace ATOOLS;
using namespace std;

Four_Particle_Observable_Base::Four_Particle_Observable_Base(std::vector<Flavour> & _flavs,
							     int _type,double _xmin,double _xmax,
							     int _nbins,std::string _name) :
  Primitive_Observable_Base(_type,_xmin,_xmax,_nbins,NULL)
{
  if (_flavs.size()<4) {
    msg.Error()<<"Error in Four_Particle_Observable_Base:"<<std::endl
	       <<"   No four flavours specified, try to copy flavours."<<std::endl;
    do { _flavs.push_back(_flavs.back()); } while (_flavs.size()<4);
  }
  std::string help = std::string("");
  for (int i=0;i<4;i++) {
    m_flavs.push_back(_flavs[i]);
    help += _flavs[i].Name();
    if (i==1) help+=std::string("--");
  }
  name       = _name + std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Four_Particle_Observable_Base::Evaluate(double _value,double _weight) {
  histo->Insert(_value,_weight); 
}

 
void Four_Particle_Observable_Base::Evaluate(int _nout,ATOOLS::Vec4D * _moms,ATOOLS::Flavour * _flavs,
					    double _weight) 
{
  for (int i=0;i<_nout;i++) { 
    if (_flavs[i]==m_flavs[0]) {
      for (int j=0;j<_nout;j++) { 
	if (_flavs[j]==m_flavs[1] && i!=j) {
	  for (int k=0;k<_nout;k++) { 
	    if (_flavs[k]==m_flavs[2] && k!=j && k!=i) {
	      for (int l=0;l<_nout;l++) { 
		if (_flavs[l]==m_flavs[3] && l!=k && l!=j && l!=i) {
		  Evaluate(_moms[i],_moms[j],_moms[k],_moms[l],_weight);
		}
	      }
	    }
	  }
	} 
      }
    }
  }
}


void Four_Particle_Observable_Base::Evaluate(const Particle_List & _plist,double _weight)
{
  for (Particle_Const_Iterator plit1=_plist.begin();plit1!=_plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flavs[0]) {
      for (Particle_Const_Iterator plit2=_plist.begin();plit2!=_plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flavs[1] && plit1!=plit2) {
	  for (Particle_Const_Iterator plit3=_plist.begin();plit3!=_plist.end();++plit3) {
	    if ((*plit3)->Flav()==m_flavs[2] && plit3!=plit2 && plit3!=plit1) {
	      for (Particle_Const_Iterator plit4=_plist.begin();plit4!=_plist.end();++plit4) {
		if ((*plit4)->Flav()==m_flavs[3] && plit4!=plit3 && plit4!=plit2 && plit4!=plit1) {
		  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),
			   (*plit3)->Momentum(),(*plit4)->Momentum(),_weight);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void Four_Particle_PlaneAngle::Evaluate(Vec4D _mom1,Vec4D _mom2,Vec4D _mom3,Vec4D _mom4,
					double _weight)
{
  Vec3D normal1 = cross(Vec3D(_mom1),Vec3D(_mom2));
  Vec3D normal2 = cross(Vec3D(_mom3),Vec3D(_mom4));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs()); 
  histo->Insert(costh,_weight); 
}
 
Four_Particle_PlaneAngle::Four_Particle_PlaneAngle(std::vector<Flavour> & _flavs,
						   int _type,double _xmin,double _xmax,int _nbins,
						   std::string _name) :
  Four_Particle_Observable_Base(_flavs,_type,_xmin,_xmax,_nbins,_name) { }


