#ifndef PHASIC_GGH_Process_Manager
#define PHASIC_GGH_Process_Manager

#include "PHASIC++/Process/Process_Base.H"

namespace PHASIC{

  class ME_Generators;

  /* 
     Class that reproduces the seriously fucked up process
     management of the Matrix_Element_Handler and
     the corresponding interplay mess. 
  */
  class GGH_Process_Manager {

    typedef std::vector<NLOTypeStringProcessMap_Map*> Map_Vector;

  private:

    /* 
       Mimics 'm_pmaps' member of Matrix_Element_Handler
       here only used for bookkeeping and memory management:
       all procs in these maps must be deleted in destructor 
    */
    Map_Vector m_maps;

    /* 
       Mimics 'm_procs' member of Matrix_Element_Handler, 
       contains processes of the form  93 93 -> 25 39
       that are actually groups of partonic subprocesses 
    */
    Process_Vector    m_procs;
    
    ME_Generators*    Generators();
    ME_Generators*    p_generators;
    void              InitializeProcess(const ATOOLS::Cluster_Amplitude& ampl, bool external);
    
  public:
    GGH_Process_Manager();
    ~GGH_Process_Manager();

    /* 
       High-level interface method for GGH_KFactor_Setter. 
       Returns a PARTONIC process corresponding to the
       Cluster_Amplitude. 'External' flag is for OpenLoops processes.
    */
    Process_Base* GetProcess(const ATOOLS::Cluster_Amplitude& ampl, bool external);
    Process_Base* GetProcess(const std::string& name              , bool external);
    inline void SetGenerators(ME_Generators* gens){ p_generators = gens;}

  };
  
}

#endif
