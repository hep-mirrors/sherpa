#ifndef REMNANTS_Main_Electron_Remnant_H
#define REMNANTS_Main_Electron_Remnant_H

#include "REMNANTS/Main/Remnant_Base.H"
#include "PDF/Main/PDF_Base.H"
#include "YFS/Main/YFS_Handler.H"

namespace REMNANTS {

  class Electron_Remnant final : public Remnant_Base {
  private:
    FlavourList     m_constituents;
    PDF::PDF_Base * p_pdfbase;
    YFS::YFS_Handler * p_yfs;
  public:
    Electron_Remnant(PDF::PDF_Base * pdf,const unsigned int & beam,
		     const unsigned int & tag=0);

    Electron_Remnant(YFS::YFS_Handler * yfs,const unsigned int & beam,
         const unsigned int & tag=0);

    bool FillBlob(Colour_Generator* colours, ATOOLS::ParticleMomMap* ktmap = nullptr,
                  const bool& copy = true) override;
    bool TestExtract(const ATOOLS::Flavour& flav, const ATOOLS::Vec4D& mom) override;
  };
}

#endif
