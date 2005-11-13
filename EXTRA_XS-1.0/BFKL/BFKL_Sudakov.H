#ifndef BFKL_Sudakov_H
#define BFKL_Sudakov_H

#include "LEV_Group.H"

namespace EXTRAXS {

  class BFKL_Sudakov {
  public:

    typedef std::vector<ATOOLS::Vec4D> Vector_Vector;

  private:

    LEV_Group *p_levs;

    double m_ya, m_yb, m_kt2min, m_kt2max, m_y, m_qt2, m_phi;
    double m_asbarmax, m_integral, m_weight;

    double AlphaSBar(const double &kt2) const;
    
  public:

    // constructor
    BFKL_Sudakov();

    // destructor
    ~BFKL_Sudakov();

    // member functions
    bool Initialize();

    bool Dice();
    bool CalculateWeight(const Vector_Vector &moms);

    // inline functions
    inline void SetYMax(const double &ymax) { m_ya=ymax; }
    inline void SetYMin(const double &ymin) { m_yb=ymin; }

    inline void SetKT2Min(const double &kt2min) { m_kt2min=kt2min; }
    inline void SetKT2Max(const double &kt2max) { m_kt2max=kt2max; }

    inline double GetY() const   { return m_y;   }
    inline double GetQT2() const { return m_qt2; }
    inline double GetPhi() const { return m_phi; }

    inline double Weight() const { return m_weight; }

  };// end of class BFKL_Sudakov

}// end of namespace EXTRAXS

#endif
