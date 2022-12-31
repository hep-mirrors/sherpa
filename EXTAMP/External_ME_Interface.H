#ifndef EXTAMP_Main_ME_Provider_H
#define EXTAMP_Main_ME_Provider_H

#include "PHASIC++/Process/ME_Generator_Base.H"

namespace PHASIC {
  class Tree_ME2_Base;
  struct External_ME_Args;
}

namespace EXTAMP {

  class Dipole_Info;
  class CS_Dipole;

  class External_ME_Interface: public PHASIC::ME_Generator_Base {
    
  private :

    BEAM::Beam_Spectra_Handler* p_beam;
    PDF::ISR_Handler* p_isr;

  public :

    typedef std::map<size_t, ATOOLS::Flavour_Vector> Combinable_Map;

    External_ME_Interface();

    ~External_ME_Interface();

    virtual int PerformTests()  {return 1;};
    virtual bool NewLibraries() {return false;};

    /* Pure virtual from PHASIC::ME_Generator_Base. Instantiates and
       initializes a high-level group-type process (e.g. 93 93 -> 90 90) */
    virtual PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add);

    /* Determine if a PARTONIC Process as specified by args */
    static bool PartonicProcessExists(const PHASIC::External_ME_Args &args);

    /* Determine if a PARTONIC Process as specified by pi exists */
    static bool PartonicProcessExists(const PHASIC::Process_Info &pi);

    /* Instantiate a PARTONIC Process as specified by pi */
    static PHASIC::Process_Base *
    InstantiatePartonicProcess(const PHASIC::Process_Info &pi);
    
    /* Pure virtual from PHASIC::ME_Generator_Base. Initialize this class. */
    virtual bool Initialize(MODEL::Model_Base *const model,
			    BEAM::Beam_Spectra_Handler *const beam,
			    PDF::ISR_Handler *const isr);

    /* Use getter to load an external matrix element */
    static PHASIC::Tree_ME2_Base* 
    GetExternalBornME(const PHASIC::External_ME_Args& args);

    /* Given indices of external flavours i_1,i_2,i_3,... and the
       combined index k=(1<<i_1) & (1<<i_2) & (1<<i_3) & ... of some
       processs, this map contains all flavours that arise from
       diagrammatic clusterings of the external flavours */
    static Combinable_Map ConstructCombinableMap(const ATOOLS::Flavour_Vector&,
						 const PHASIC::Process_Info&,
						 const size_t& n_in);
  };

}

#endif
