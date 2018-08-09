#ifndef CFPSHOWER_Tools_Weight_H
#define CFPSHOWER_Tools_Weight_H

#include <vector>

namespace CFPSHOWER {
  class Parton;
  
  struct Weight {
    double m_f, m_g, m_h;
  public:
    inline Weight(const double &f,const double &g,const double &h):
      m_f(f), m_g(g), m_h(h) {}
    inline const double operator()() const { return m_f/m_g; }
    inline const double Accept()     const { return m_g/m_h; }
    inline const double Reject()     const { return m_g/m_h*(m_h-m_f)/(m_g-m_f); }
  };// end of struct Weight

  inline std::ostream &operator<<(std::ostream &s,const Weight &w)
  { return s<<"W(f="<<w.m_f<<",g="<<w.m_g<<",h="<<w.m_h<<")"; }
}

#endif
