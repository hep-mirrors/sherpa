#ifndef SHRIMPS_Main_Shrimps_H
#define SHRIMPS_Main_Shrimps_H

#include "SHRiMPS/Main/Cluster_Algorithm.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Cross_Sections/Cross_Sections.H"
#include "SHRiMPS/Event_Generation/Event_Generator.H"
#include "SHRiMPS/Beam_Remnants/Remnant_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Flavour.H"
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
    weight_mode::code m_weightmode;
    deqmode::code     m_deq;
    int               m_test;
    Cross_Sections  * p_xsecs;
    Remnant_Handler * p_remnants;
    Event_Generator * p_generator;
    Cluster_Algorithm m_cluster;
    
    void GenerateXsecs();
    void ReadEnergiesFromFile(std::set<double> & energies,
			      std::string infile);
    void WriteOutElasticsYodaFile(const double & energy,std::string dirname);
    void WriteOutXSecsYodaFile(const std::set<double> & energies_tot,
			       const std::set<double> & energies_inel,
			       const std::set<double> & energies_el,
			       const std::set<double> & energies_sd,
			       const std::set<double> & energies_dd,
			       const std::vector<double> & xsectot,
			       const std::vector<double> & xsecinel,
			       const std::vector<double> & xsecelas,
			       const std::vector<double> & xsecsd,
			       const std::vector<double> & xsecdd,
			       std::string dirname);
    
    void InitialiseTheRun(BEAM::Beam_Spectra_Handler *const beam,
			  PDF::ISR_Handler *const isr);
    void InitialiseFormFactors();
    void InitialiseSingleChannelEikonals();
    void InitialiseRemnants(BEAM::Beam_Spectra_Handler *const beam,
			    PDF::ISR_Handler *const isr);
    void InitialiseTheEventGenerator();
    
    void TestShrimps(BEAM::Beam_Spectra_Handler *const beam,
		     PDF::ISR_Handler *const isr);
    void PrintAlphaS(const std::string & dirname);
    void PrintPDFs(const std::string & dirname);
    void TestEikonalGrids(const std::string & dirname);
    void TestCrossSections(const std::string & dirname);
    void TestEventGeneration(const std::string & dirname);
  public:
    Shrimps(BEAM::Beam_Spectra_Handler *const,
            PDF::ISR_Handler *const);
    ~Shrimps();

    int  GenerateEvent(ATOOLS::Blob_List * blobs);
    void CleanUp(const size_t & mode=0);

    ATOOLS::Return_Value::code MakeBeamBlobs(ATOOLS::Blob_List * blobs);
    
    inline int    ShowerMode()   const { return 1; }// MBpars("shower_mode"); }
    inline double ShowerFac()    const { return 1; }//MBpars("kt2_factor"); }
    inline double ShowerMinKT2() const { return p_generator->MinKT2(); }
    inline const double & XSec() const { return p_generator->XSec(); } 
    inline double Smin()         const { return p_generator->Smin(); }
    inline bool   IsLastRescatter() const { 
      return p_generator->IsLastRescatter(); 
    }
    inline double LadderTMax() const {return p_generator->TMax(); }  
    inline int    NLadders()   const {return p_generator->NLadders(); }  
    inline Cluster_Algorithm * GetClusterAlgorithm() { return &m_cluster; }
  };
}
#endif
