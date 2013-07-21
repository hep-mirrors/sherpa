#ifndef ATOOLS__Math__Marsaglia_H
#define ATOOLS__Math__Marsaglia_H

#include <iostream>

namespace ATOOLS {

  /*
    Combined KISS and SWB generator by George Marsaglia

    http://groups.google.co.uk/group/sci.stat.math/msg/edcb117233979602
    http://groups.google.co.uk/group/sci.math.num-analysis/msg/eb4ddde782b17051
   */

  class Marsaglia {
  public:

    typedef unsigned int UL;

  private:

    UL m_z, m_w, m_jsr, m_jcong, m_a, m_b;
    UL m_t[256], m_x, m_y, m_bro;

    unsigned char m_c;

  public:

    Marsaglia();

    void Init(UL i1,UL i2,UL i3,UL i4,UL i5,UL i6);
  
    double Get();

    void WriteStatus(std::ostream &str);
    bool ReadStatus(std::istream &str);

  };// end of class Marsaglia

}// end of namespace ATOOLS

#endif
