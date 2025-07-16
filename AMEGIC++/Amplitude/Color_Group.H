#ifndef AMEGIC_Amplitude_Color_Group_H
#define AMEGIC_Amplitude_Color_Group_H

#include "AMEGIC++/Amplitude/Amplitude_Group.H"

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












