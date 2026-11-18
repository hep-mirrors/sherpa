#ifndef GoSam_Interface_H
#define GoSam_Interface_H

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Terminator_Objects.H"

#include "PHASIC++/Process/External_ME_Args.H"

namespace MODEL {
  class Model_Base;
}

namespace GoSam {
  extern "C" void OLP_Start(const char * filename, int* success);
  extern "C" void OLP_EvalSubProcess(int,double*,double,double*,double*);
  extern "C" void OLP_EvalSubProcess2(int*,double*,double*,double*,double*);
  extern "C" void OLP_Option(const char * assignment, int* success);
  extern "C" void OLP_PrintParameter(const char * filename);
  extern "C" void OLP_SetParameter(const char* para, double* re, double* im, int* ierr);

  class GoSam_Interface: public PHASIC::ME_Generator_Base, public ATOOLS::Terminator_Object {
  private:
    static std::string s_contractName;
    static bool s_newLibraries;
    static bool s_contract;
    static std::string s_type;

    std::string m_contractBaseName;

    int m_drmode;
    std::string m_widthscheme;

  public:

    static int SetParamter(std::string para, std::complex<double> val, int errmode=1) {
      int ierr;
      double re = val.real();
      double im = val.imag();
      OLP_SetParameter(para.c_str(),&re,&im,&ierr);
      if(errmode > 0) {
        if(ierr == 0)
          THROW(fatal_error,"Failed to set parameter "+para+".");
        if(errmode > 1 && ierr > 1)
          THROW(fatal_error,"Failed to set parameter "+para+" with code "+ATOOLS::ToString(ierr));
      }
      return ierr;
    }

    static int SetParameter(std::string para, double val, int errmode=1) {
      int ierr;
      double im = 0;
      OLP_SetParameter(para.c_str(),&val,&im,&ierr);
      if(errmode > 0) {
        if(ierr == 0)
          THROW(fatal_error,"Failed to set parameter "+para+".");
        if(errmode > 1 && ierr > 1)
          THROW(fatal_error,"Failed to set parameter "+para+" with code "+ATOOLS::ToString(ierr));
      }
      return ierr;
    }

    static int SetMass(kf_code sherpaCode, std::string goSamCode="",int errmode=1) {
      if(goSamCode == "0") return errmode;
      std::string goSamName = "mass("+(goSamCode.empty()?ATOOLS::ToString(sherpaCode):goSamCode)+")";
      SetParameter(goSamName,ATOOLS::Flavour(sherpaCode).Mass(),errmode);
      return errmode;
    }

    static int SetWidth(kf_code sherpaCode, std::string goSamCode="",int errmode=1) {
      if(goSamCode == "0") return errmode;
      std::string goSamName = "width("+(goSamCode.empty()?ATOOLS::ToString(sherpaCode):goSamCode)+")";
      SetParameter(goSamName,ATOOLS::Flavour(sherpaCode).Width(),errmode);
      return errmode;
    }

    
    GoSam_Interface();
    void RegisterDefaults();
    bool Initialize(MODEL::Model_Base *const model,
                    BEAM::Beam_Spectra_Handler *const beam,
                    PDF::ISR_Handler *const isr,
                    YFS::YFS_Handler *const yfs);
    static int RegisterProcess(const PHASIC::Process_Info& pi,
                               const ATOOLS::Flavour_Vector& flavs,bool active);
    static int RegisterProcess(const PHASIC::External_ME_Args& args,bool active);
    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi,
                                            bool add)
      {
        return NULL;
      }
    bool NewLibraries() { return s_newLibraries; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
      (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }

    int  PerformTests();
    void PrepareTerminate();
  };
}

#endif
