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
		 const int scalescheme,
		 const int kfactorscheme,
		 const double scalefactor);

    // destructor 
    ~Semihard_QCD();

    // member functions
    void CreateFSRChannels();
    void InitIntegrators();

    // inline functions
    inline int FSRMode() const                { return m_fsrmode; }
    inline void SetFSRMode(const int fsrmode) { m_fsrmode=fsrmode; }

    inline PHASIC::Channel_Interface *const FSRInterface() const
    { return p_fsrinterface; }
    inline void SetFSRInterface(PHASIC::Channel_Interface *const fsr)
    { p_fsrinterface=fsr; }

  };// end of class Semihard_QCD

}// end of namespace AMISIC

#endif
