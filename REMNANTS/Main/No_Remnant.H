#ifndef REMNANTS_Main_No_Remnant_H
#define REMNANTS_Main_No_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "PDF/Main/ISR_Handler.H"

namespace REMNANTS {
  class No_Remnant final : public Remnant_Base {
  public:
    No_Remnant(const size_t & beam,const size_t & tag=0);
    bool FillBlob(Colour_Generator* colours, ATOOLS::ParticleMomMap* ktmap = NULL,
                  const bool& copy = true) override;
    bool TestExtract(const ATOOLS::Flavour& flav, const ATOOLS::Vec4D& mom) override;
  };
}

#endif
