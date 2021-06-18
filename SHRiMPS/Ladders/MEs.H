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
    MEs(Sigma_Partonic * sigma,
	const double & smin=-1., const double & tmin=0.);
    ~MEs() {};
    double operator()(Ladder * ladder,const double & qt2min);
  };
}
#endif
