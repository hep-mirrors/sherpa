#ifndef Color_Group_h
#define Color_Group_h

#include "Amplitude_Group.H"

namespace AMEGIC {
  class Color_Group: public Amplitude_Group {
  public:
    Color_Group() { groupname = std::string("Color"); }
    ~Color_Group() {
      for (size_t i=0;i<graphs.size();i++) if (graphs[i]->IsGroup()) delete graphs[i];
      graphs.clear();  
    }
  };
}
#endif












