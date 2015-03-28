#ifndef OpenLoops_Interface_H
#define OpenLoops_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"


namespace OpenLoops {

  class OpenLoops_Interface: public PHASIC::ME_Generator_Base,
                             public ATOOLS::Terminator_Object {

    static std::string s_olprefix;
    static bool        s_ignore_model;

  public :
    OpenLoops_Interface() : ME_Generator_Base("OpenLoops") {}
    ~OpenLoops_Interface();

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    static int RegisterProcess(const PHASIC::Subprocess_Info& is,
                               const PHASIC::Subprocess_Info& fs,
                               int amptype);
    static void EvaluateTree(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res);
    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res, METOOLS::DivArrD& virt);
    static void EvaluateLoop2(int id, const ATOOLS::Vec4D_Vector& momenta,
                              double& res);

    static void SetParameter(const std::string & key, double value);
    static void SetParameter(const std::string & key, int value);
    static void SetParameter(const std::string & key, std::string value);
    static double GetDoubleParameter(const std::string & key);
    static int GetIntParameter(const std::string & key);

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add)
    { return NULL; }
    bool NewLibraries() { return false; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }

    int  PerformTests();
    void PrepareTerminate();
  };


}

#endif
