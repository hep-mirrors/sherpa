#ifndef SHRIMPS_Main_Shrimps_H
#define SHRIMPS_Main_Shrimps_H

#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Eikonals/Eikonal_Creator.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "SHRiMPS/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include <list>

namespace PDF    { class ISR_Handler; }
namespace BEAM   { class Beam_Spectra_Handler; }

namespace SHRIMPS {
  /*!
    \class Shrimps
    \brief Central class
    
    \var m_test
    Steers the testing (with SHRIMPS::m_test):
    -  m_test=-2:  Prints values of various Bessel functions from
                   SHRIMPS::Special_Functions to file.
		   Test method in SHRIMPS::Form_Factor
    -  m_test=-1:  Checks the SHRIMPS::Form_Factor, and in particular the
                   Fourier Transform, some results on screen, some on
		   file.  Test method in SHRIMPS::Form_Factor.
    -  m_test=1:   Checks the eikonal evaluation with a Gaussian form
                   factor and no damping/recombination in the DEqs 
		   (lambda=0).  Then the results are fully analytic.
		   Test function in SHRIMPS::Eikonal_Creator.
  */
  class Shrimps {
  private:
    run_mode::code    m_runmode; 
    weight_mode::code m_weightmode;
    deqmode::code     m_deq;
    int               m_NGWstates;
    int               m_test;
    
    std::vector<Form_Factor>   m_ffs;
    std::list<Omega_ik *>      m_eikonals;
    Cross_Sections             m_cross;
    std::vector<Continued_PDF> m_pdfs;
    Beam_Remnant_Handler     * p_beamremnants;
    Event_Generator          * p_generator;
    
    void GenerateXsecs();
    
    void InitialiseFormFactors();
    void InitialiseSingleChannelEikonals(const double & Ecms);
    void InitialiseCrossSections(const double & energy);
    void InitialiseBeamRemnants(BEAM::Beam_Spectra_Handler *const beam,
				PDF::ISR_Handler *const isr);
    void InitialiseEventGenerator();
    void PrintPDFandAlphaS();
  public:
    Shrimps(ATOOLS::Data_Reader * dr,
	  BEAM::Beam_Spectra_Handler *const beam,
	  PDF::ISR_Handler *const isr);
    ~Shrimps();

    int  GenerateEvent(ATOOLS::Blob_List * blobs);
    void CleanUp(const size_t & mode=0);

    ATOOLS::Return_Value::code FillBeamBlobs(ATOOLS::Blob_List * blobs);

    inline void   SetBeamBlob(ATOOLS::Blob *const beamblob,const int & beam) {
      p_beamremnants->SetBeamBlob(beamblob,beam);
    }
    inline int    ShowerMode()   const { return MBpars("shower_mode"); }
    inline double ShowerFac()    const { return MBpars("kt2_factor"); }
    inline double ShowerMinKT2() const { return p_generator->ShowerMinKT2(); }
    inline const double & XSec() const { return p_generator->XSec(); } 
    inline const double   Smin() const { return p_generator->Smin(); }
    inline const bool     IsLastRescatter() const { 
      return p_generator->IsLastRescatter(); 
    }
  };
}
#endif
