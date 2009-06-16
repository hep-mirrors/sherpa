#ifndef AMEGIC_Amplitude_Prop_Generator_H
#define AMEGIC_Amplitude_Prop_Generator_H

#include "AMEGIC++/Amplitude/Zfunc.H"

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
