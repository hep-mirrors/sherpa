#ifndef BFKL_Sudakov_H
#define BFKL_Sudakov_H

#include "LEV_Group.H"

namespace EXTRAXS {

  class BFKL_Sudakov {
  public:

    typedef std::vector<ATOOLS::Vec4D> Vector_Vector;

  private:

    LEV_Group *p_levs;

    double m_ya, m_yb, m_kt2min, m_kt2max, m_y, m_kt2, m_phi;
    double m_asmax, m_integral, m_gamma, m_oldy;

  public:

    // constructor
    BFKL_Sudakov();

    // destructor
    ~BFKL_Sudakov();

    // member functions
    bool Initialize();

    bool Dice();
    bool Approve(const Vector_Vector &k,
		 const ATOOLS::Vec4D &q1,const ATOOLS::Vec4D &q2);

    // inline functions
    inline void SetYMax(const double &ymax) { m_ya=ymax; }
    inline void SetYMin(const double &ymin) { m_yb=ymin; }

    inline void SetKT2Min(const double &kt2min) { m_kt2min=kt2min; }
    inline void SetKT2Max(const double &kt2max) { m_kt2max=kt2max; }

    inline double GetY() const   { return m_y;   }
    inline double GetKT2() const { return m_kt2; }
    inline double GetPhi() const { return m_phi; }

    inline LEV_Base *Selected() const { return p_levs->Selected(); }

  };// end of class BFKL_Sudakov

}// end of namespace EXTRAXS

#endif
