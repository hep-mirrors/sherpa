#ifndef Electron_Remnant_H
#define Electron_Remnant_H

#include "Remnant_Base.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class Electron_Remnant: public Remnant_Base {
  private:

    PDF::PDF_Base *p_pdfbase;

  public:
    
    // constructor
    Electron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);
    bool AdjustKinematics();

    ATOOLS::Flavour ConstituentType(const ATOOLS::Flavour &flavour);

  };// end of class Remnant_Base

}// end of namespace SHERPA

#endif
