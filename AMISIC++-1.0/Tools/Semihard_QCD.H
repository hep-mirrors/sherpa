#ifndef Semihard_QCD_H
#define Semihard_QCD_H

#include "XS_Group.H"
#include "Channel_Interface.H"

namespace AMISIC {

  class Semihard_QCD: public EXTRAXS::XS_Group {
  private:

    int m_fsrmode;

    PHASIC::Channel_Interface *p_fsrinterface;

  public:

    // constructor
    Semihard_QCD(BEAM::Beam_Spectra_Handler *const beam,
		 PDF::ISR_Handler *const isr,
		 ATOOLS::Selector_Data *const seldata,
		 const ATOOLS::Flavour *flavours,
		 const int scalescheme,
		 const int kfactorscheme);

    // destructor 
    ~Semihard_QCD();

    // member functions
    void CreateFSRChannels();
    void CreateISRChannels();
    void InitIntegrators();

    // inline functions
    inline int FSRMode() const                { return m_fsrmode; }
    inline void SetFSRMode(const int fsrmode) { m_fsrmode=fsrmode; }

    inline PHASIC::Channel_Interface *const FSRInterface() const
    { return p_fsrinterface; }
    inline void SetFSRInterface(PHASIC::Channel_Interface *const fsr)
    { p_fsrinterface=fsr; }

  };// end of class Semihard_QCD

  //////////////////
  // Doxygen part //
  //////////////////

  /*!
    \file Semihard_QCD.H
    \brief Declares the class Semihard_QCD
  */

  /*!
    \class Semihard_QCD
    \brief Container class for processes used in multiple 
    interaction simulations

    This container is derived from EXTRAXS::XS_Group. It stores 
    the processes used during multiple interaction simulations.
    Semihard_QCD implements a special setup for the phasespace, 
    i.e. it generates only one isr channel, which is suitably 
    designed.
  */

}// end of namespace AMISIC

#endif
