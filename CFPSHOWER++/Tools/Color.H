#ifndef CFPShower_Tools_Color_H
#define CFPShower_Tools_Color_H

#include "ATOOLS/Phys/Flow.H"
#include <vector>

namespace CFPSHOWER {
  struct Color : public std::vector<unsigned int> {
  public:
    Color(const unsigned int & trip=0, const unsigned int & anti=0) {
      resize(2);
      (*this)[0] = trip; (*this)[1] = anti;
    }
    Color(const Color & color) {
      resize(2);
      for (size_t i=0;i<2;i++) (*this)[i] = color[i];
    }
  };
  std::ostream &operator<<(std::ostream & s,const Color & color);
  typedef std::vector<Color> Color_Vector;
}

#endif
