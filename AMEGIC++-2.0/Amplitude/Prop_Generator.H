#ifndef Prop_Generator_h
#define Prop_Generator_h

#include "Zfunc.H"

namespace AMEGIC {
  class Point;
  
  class Prop_Generator {
    Pfunc_List plist;
  public:
    Prop_Generator() {}
    void Convert(Point*);
    void Fill();
    void Kill(Zfunc_List&); 
    void Get(Pfunc_List& _plist);
  };
}
#endif
