#ifndef DIRE__Tools__Color_H
#define DIRE__Tools__Color_H

#include <iostream>
#include <vector>

namespace DIRE {

  struct Color {
    int m_i, m_j, m_n;
    double m_w;
    inline Color(const int i=-1,const int j=-1,
		 const double &w=1.0,const int n=0):
      m_i(i), m_j(j), m_w(w), m_n(n) {}
  };// end of struct Color

  std::ostream &operator<<(std::ostream &s,const Color &c);

  typedef std::vector<Color> Color_Vector;

}// end of namespace DIRE

#endif
