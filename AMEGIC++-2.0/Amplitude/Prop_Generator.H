#ifndef Prop_Generator_h
#define Prop_Generator_h

#include "Pfunc.H"
#include "Vector.H"
#include "Zfunc.H"

namespace AMEGIC {

  class Prop_Generator {
    std::list<Pfunc*> plist;
  public:
    Prop_Generator() {}
    void Convert(Point*);
    void Fill();
    void Kill(std::list<Zfunc*>&); 
    void Calculate(AMATOOLS::Vec4D*,int*);
    void Get(std::list<Pfunc*>& _plist) {
      for (std::list<Pfunc*>::iterator pit=plist.begin();pit!=plist.end();++pit) 
	_plist.push_back(*pit);
    }
  };
}
#endif
