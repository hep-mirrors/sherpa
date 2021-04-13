#ifndef SHRIMPS_Ladders_Ladder_TProp_H
#define SHRIMPS_Ladders_Ladder_TProp_H

#include "ATOOLS/Math/Vector.H"
#include <list>

namespace SHRIMPS {
  struct colour_type {
    enum code {
      none    = 0,
      singlet = 1,
      triplet = 3,
      octet   = 8
    };
  };
  std::ostream & operator<<(std::ostream & s, const colour_type::code & colour);
  
  class T_Prop {
  private:
    colour_type::code m_col;
    ATOOLS::Vec4D     m_q;
    double            m_q2, m_qt2, m_q02, m_qmax2;
  public:
    T_Prop(const colour_type::code & col=colour_type::octet,
	   const ATOOLS::Vec4D & q=ATOOLS::Vec4D(0.,0.,0.,0.),
	   const double & q02=0.);

    inline void SetCol(const colour_type::code & col) { m_col   = col; }
    inline void SetQ(const ATOOLS::Vec4D & q)         { m_q     = q; }
    inline void SetQ2(const double & q2)              { m_q2    = q2; }
    inline void SetQT2(const double & qt2)            { m_qt2   = qt2; }
    inline void SetQ02(const double & q02)            { m_q02   = q02; }
    inline void SetQMax2(const double & qmax2)        { m_qmax2 = qmax2; }
    
    inline const colour_type::code & Col() const { return m_col; }
    inline const ATOOLS::Vec4D & Q()       const { return m_q; }
    inline const double & Q2()             const { return m_q2; }
    inline const double & QT2()            const { return m_qt2; }
    inline const double & Q02()            const { return m_q02; }
    inline const double & QMax2()          const { return m_qmax2; }
  };
  
  typedef std::list<T_Prop> TPropList;
  std::ostream & operator<<(std::ostream & s, const T_Prop & tprop);
  std::ostream & operator<<(std::ostream & s, const TPropList & props);
}
#endif
