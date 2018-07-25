#ifndef REMNANTS_Main_Electron_Remnant_H
#define REMNANTS_Main_Electron_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "PDF/Main/PDF_Base.H"

namespace REMNANTS {

  class Electron_Remnant: public Remnant_Base {
  private:
    FlavourList     m_constituents;
    PDF::PDF_Base * p_pdfbase;
  public:
    Electron_Remnant(PDF::PDF_Base * pdf,const unsigned int _m_beam);

    bool FillBlob(ATOOLS::ParticleMomMap *ktmap=NULL,const bool & copy=true);
    bool TestExtract(const ATOOLS::Flavour &flav,const ATOOLS::Vec4D &mom);
  };
}

#endif
