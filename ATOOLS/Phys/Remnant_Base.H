#ifndef ATOOLS_Phys_Remnant_Base_H
#define ATOOLS_Phys_Remnant_Base_H

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vec4.H"

namespace BEAM { class Beam_Base; }

namespace ATOOLS {
  class Basic_Remnant_Base {
  protected:
    BEAM::Beam_Base * p_beam;
    ATOOLS::Blob    * p_beamblob;
    ATOOLS::Vec4D     m_pbeam;

    unsigned int m_beam;
  public:
    Basic_Remnant_Base(const unsigned int beam) :
      m_beam(beam), p_beam(nullptr), p_beamblob(nullptr) {};
    ~Basic_Remnant_Base() {}
    virtual bool TestExtract(const ATOOLS::Flavour &flav,
			     const ATOOLS::Vec4D &mom) { return true; }
    virtual void Reset() {}
    inline virtual void SetBeam(BEAM::Beam_Base * beam) { p_beam = beam; }
    inline virtual BEAM::Beam_Base * GetBeam() const    { return p_beam; }
  };
}

#endif
