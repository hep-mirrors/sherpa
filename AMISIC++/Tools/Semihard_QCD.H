#ifndef AMISIC_Tools_Semihard_QCD_H
#define AMISIC_Tools_Semihard_QCD_H
/*!
  \file Semihard_QCD.H
  \brief Declares the class Semihard_QCD
*/

#include "EXTRA_XS/Main/Process_Group.H"
#include "PHASIC++/Channels/Channel_Interface.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

namespace AMISIC {

  class Semihard_QCD: public EXTRAXS::Process_Group,
		      public PHASIC::ME_Generator_Base {
  private:

    int m_fsrmode;

    PHASIC::Channel_Interface *p_fsrinterface;

  public:

    // constructor
    Semihard_QCD(ATOOLS::Data_Reader *const read);

    // destructor 
    ~Semihard_QCD();

    // member functions
    void CreateFSRChannels();
    void CreateISRChannels();
    void InitIntegrators();

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr) { return true; }
    Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }
    int PerformTests() { return 1; }
    bool NewLibraries() { return false; }

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}

    ATOOLS::Cluster_Amplitude *
    ClusterConfiguration(PHASIC::Process_Base *const proc,
			 const size_t &mode=0)
    { return NULL; }

    // inline functions
    inline int FSRMode() const                { return m_fsrmode; }
    inline void SetFSRMode(const int fsrmode) { m_fsrmode=fsrmode; }

    inline PHASIC::Channel_Interface *FSRInterface() const
    { return p_fsrinterface; }
    inline void SetFSRInterface(PHASIC::Channel_Interface *const fsr)
    { p_fsrinterface=fsr; }

  };// end of class Semihard_QCD
  /*!
    \class Semihard_QCD
    \brief Container class for processes used in multiple 
    interaction simulations

    This container is derived from EXTRAXS::Process_Group. It stores 
    the processes used during multiple interaction simulations.
    Semihard_QCD implements a special setup for the phasespace, 
    i.e. it generates only one isr channel, which is suitably 
    designed.
  */

}// end of namespace AMISIC

#endif
