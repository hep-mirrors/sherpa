#ifndef GoSam_Interface_H
#define GoSam_Interface_H

#include "ATOOLS/Org/Terminator_Objects.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"


namespace GoSam {

  class GoSam_Interface: public PHASIC::ME_Generator_Base,
                         public ATOOLS::Terminator_Object {

    static std::string s_gosamprefix;
    static bool        s_ignore_model;
    static bool        s_exit_on_error;
    static bool        s_newlibs;

    void RegisterDefaults() const;

  public :
    GoSam_Interface();
    ~GoSam_Interface();

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    static int RegisterProcess(const PHASIC::Subprocess_Info& is,
                               const PHASIC::Subprocess_Info& fs);
    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& moms,
                             double& mu,
                             double& born, METOOLS::DivArrD& virt,
                             double& accu);

    static void SetParameter(const std::string & key, Complex value);
    static void SetParameter(const std::string & key, double value);
    static void SetParameter(const std::string & key, int value);
    static void SetParameter(const std::string & key, std::string value);
    static bool ExitOnError() { return s_exit_on_error; }

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add)
    { return NULL; }
    bool NewLibraries() { return s_newlibs; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }

    int  PerformTests();
    void PrepareTerminate();

    static std::map<int,std::string> s_procmap;
    static size_t      s_vmode;
  };


}

#endif
