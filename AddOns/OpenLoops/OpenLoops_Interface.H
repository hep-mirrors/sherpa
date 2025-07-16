#ifndef OpenLoops_Interface_H
#define OpenLoops_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "MODEL/Main/Model_Base.H"


namespace ATOOLS {

  extern KF_Table s_kftable;

}

namespace MODEL {

  class Model_Base;

}

namespace OpenLoops {

  enum AmplitudeType {
    Tree = 1,
    Loop2 = 12
  };

  class OpenLoops_Interface: public PHASIC::ME_Generator_Base,
                             public ATOOLS::Terminator_Object {

    static std::string s_olprefix;
    static bool        s_ignore_model;
    static bool        s_exit_on_error;
    static bool        s_ass_func;
    static int         s_ass_ew;
    static std::map<std::string, std::string> s_evgen_params;

    void RegisterDefaults() const;
    int TranslatedEWRenormalisationScheme() const;

  public :
    OpenLoops_Interface();
    ~OpenLoops_Interface();

    bool Initialize(MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr,
		    YFS::YFS_Handler *const yfs);

    static void SwitchMode(const int mode);

    static int RegisterProcess(const PHASIC::Subprocess_Info& is,
                               const PHASIC::Subprocess_Info& fs,
                               int amptype);
    static int RegisterProcess(const ATOOLS::Flavour_Vector& infl,
                               const ATOOLS::Flavour_Vector& oufl,
                               int amptype);
    static void EvaluateTree(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res);
    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res, METOOLS::DivArrD& virt);
    static void EvaluateLoop2(int id, const ATOOLS::Vec4D_Vector& momenta,
                              double& res);
    static void EvaluateAssociated(int id, const ATOOLS::Vec4D_Vector& momenta,
                                   int ass, double& res);
    static int ConvertAssociatedContributions(const asscontrib::type at);

    /* Returns <1,2,...|p_{mu_i}T_iT_j p_{nu_i}|1,2,...> / p^2 with
       i=emitter, j=spectator */
    static double EvaluateSpinCorrelator(int id, const Vec4D_Vector& momenta,
					 const Vec4D& p,
					 size_t emitter, size_t spectator,
					 AmplitudeType type);
    
    /* Evaluate all color correlators. Return the one corresponding to
       <1,2,...|TiTj|...,2,1>, with i=emitter, j=spectator */
    static double EvaluateColorCorrelator(int id, const Vec4D_Vector& momenta,
					  size_t emitter, size_t spectator,
					  AmplitudeType type);
    
    /* Evaluate all color correlators and fill matrix such that
       ccmatrix[i][j] = <1,2,...|TiTj|...,2,1> */
    static void
    PopulateColorCorrelatorMatrix(int id, const Vec4D_Vector& momenta,
				  double& born2, double* ccmatrix,
				  AmplitudeType type);

    static void SetParameter(const std::string & key, double value);
    static void SetParameter(const std::string & key, int value);
    static void SetParameter(const std::string & key, std::string value);
    static double GetDoubleParameter(const std::string & key);
    static int GetIntParameter(const std::string & key);

    /* Propagate model parameters to OpenLoops in SM */
    void SetParametersSM(const MODEL::Model_Base* model);
    
    /* Propagate model parameters to OpenLoops in UFO models */
    void SetParametersUFO(const MODEL::Model_Base* model);

    static bool ExitOnError() { return s_exit_on_error; }

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add)
    { return NULL; }
    bool NewLibraries() { return false; }
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
