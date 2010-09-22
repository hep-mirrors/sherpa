#ifndef ATOOLS_Math_ZAlign_H
#define ATOOLS_Math_ZAlign_H

#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS {

  class ZAlign {
  private:

    Poincare m_cmso, m_zrot, m_cmsn;
    Vec4D    m_pao, m_pan, m_pb;

    int m_stat;

  public:

    inline ZAlign() {}

    ZAlign(const Vec4D &pa,const Vec4D &pb,
	   const double &ma2,const double &mb2,
	   const int mode=0);

    void  Align(Vec4D &p) const;
    Vec4D Align(const Vec4D &p) const;

    inline const Poincare &OldCMS() const { return m_cmso; }
    inline const Poincare &NewCMS() const { return m_cmsn; }
    inline const Poincare &Rotate() const { return m_zrot; }

    inline const Vec4D &PaOld() const { return m_pao; }
    inline const Vec4D &PaNew() const { return m_pan; }
    inline const Vec4D &Pb() const    { return m_pb;  }

    inline int Status() const { return m_stat; }

  };// end of class ZAlign

}// end of namespace ATOOLS

#endif
