#ifndef REMNANTS_Main_No_Remnant_H
#define REMNANTS_Main_No_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "PDF/Main/ISR_Handler.H"

namespace REMNANTS {
  class No_Remnant: public Remnant_Base {
  public:    
    No_Remnant(const unsigned int _m_beam);
    bool FillBlob(ATOOLS::ParticleMomMap *ktmap=NULL,const bool & copy=true);
    bool TestExtract(const ATOOLS::Flavour &flav,const ATOOLS::Vec4D &mom);
  };
}

#endif
