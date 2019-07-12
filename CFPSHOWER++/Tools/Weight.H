#ifndef CFPSHOWER_Tools_Weight_H
#define CFPSHOWER_Tools_Weight_H

#include <vector>

namespace CFPSHOWER {
  class Parton;
  
  struct Weight {
    double m_value, m_estover, m_realover;
  public:
    inline Weight(const double &value,const double &estover,const double &realover):
      m_value(value), m_estover(estover), m_realover(realover) {}
    inline const double operator()() const { return m_value/m_estover; }
    inline const double Accept()     const { return m_estover/m_realover; }
    inline const double Reject()     const {
      return m_estover/m_realover*(m_realover-m_value)/(m_estover-m_value);
    }
  };// end of struct Weight

  inline std::ostream &operator<<(std::ostream &s,const Weight &w)
  { return s<<"W(value = "<<w.m_value<<", "
	    <<"estimated over = "<<w.m_estover<<", "
	    <<"real over = "<<w.m_realover<<")"; }
}

#endif
