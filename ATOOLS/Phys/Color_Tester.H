#ifndef ATOOLS_Phys_Color_Tester_H
#define ATOOLS_Phys_Color_Tester_H

#include "ATOOLS/Phys/Parton_Finder.H"

namespace ATOOLS {
  
  class Color_Tester: public Parton_Tester {
  private:

    std::pair<unsigned int,unsigned int> m_color;
    
  public:

    // constructor
    Color_Tester(const unsigned int index,const unsigned int color);

    // member functions
    void Turn();
    bool Test(const Particle *parton) const;
    
  };

}// end of namespace ATOOLS

#endif
