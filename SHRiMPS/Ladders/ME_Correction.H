#ifndef SHRIMPS_Ladders_ME_Correction_H
#define SHRIMPS_Ladders_ME_Correction_H

#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"

namespace SHRIMPS {
  class ME_Correction {
  private:
    Continued_PDF * p_pdf[2];
    double          m_tmin, m_Ecms, m_S, m_Ebeam[2];

    void CalcPDFs(const double & x1,const double & x2,const double & mu2);
  public:
    ME_Correction(const double & tmin);
    ~ME_Correction() {};

    double operator()(Ladder * ladder);
    inline void SetPDFs(Continued_PDF * pdf0,Continued_PDF * pdf1) { p_pdf[0]=pdf0; p_pdf[1]=pdf1; }
  };
}
#endif
