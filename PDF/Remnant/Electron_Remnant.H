#ifndef PDF_Remnant_Electron_Remnant_H
#define PDF_Remnant_Electron_Remnant_H

#include "PDF/Remnant/Remnant_Base.H"
#include "PDF/Main/ISR_Handler.H"

namespace PDF {

  class Electron_Remnant: public Remnant_Base {
  private:

    PDF::PDF_Base *p_pdfbase;

  public:
    
    // constructor
    Electron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam);

    // member functions
    bool FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist);
    bool AdjustKinematics();

  };// end of class Remnant_Base

}// end of namespace PDF

#endif
