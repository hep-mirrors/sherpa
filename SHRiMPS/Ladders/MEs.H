#ifndef SHRIMPS_Ladders_ME_Correction_H
#define SHRIMPS_Ladders_ME_Correction_H

#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"

namespace SHRIMPS {
  class MEs {
  private:
    Sigma_Partonic * p_sigma;
    double m_shatmin, m_thatmin;
  public:
    MEs(const double & smin=-1., const double & tmin=0.);
    ~MEs() {};
    void   SetPartonic(Sigma_Partonic * sigma);
    double operator()(Ladder * ladder,const double & qt2min);
    double operator()(Ladder * ladder);
    double PDFratio(const ATOOLS::Vec4D & qprev,const ATOOLS::Flavour & fprev,
		    const ATOOLS::Vec4D & qact,const ATOOLS::Flavour & fact,
		    const size_t & dir);
  };
}
#endif
