#ifndef BEAM_Main_Lepton_H
#define BEAM_Main_Lepton_H

#include "BEAM/Main/Beam_Base.H"


namespace BEAM {
  class Lepton_Beam : public Beam_Base {
  private:
    const ATOOLS::Flavour m_none;
  public:
    Lepton_Beam(const ATOOLS::Flavour,const double,const double,const int);
    Beam_Base * Copy();
    ATOOLS::Flavour Remnant();
    // Void  Gaussian();
    bool        CalculateWeight(const double,const double) { return true; }
  };
  /*!
    \file
    \brief contains the class BEAM::Lepton_Beam.
  */
  /*!
    \class Lepton_Beam
    \brief This class is for beam spectra including beam energy spread
  */
}
#endif
