#ifndef Recola_Interface_H
#define Recola_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"

#include "recola.h"

namespace Recola {

  class Recola_Interface: public PHASIC::ME_Generator_Base,
                             public ATOOLS::Terminator_Object {

    static std::string  s_recolaprefix;
    static bool         s_ignore_model;
    static bool         s_exit_on_error;
    static double       s_light_fermion_threshold;
    static unsigned int s_recolaProcIndex;
    static bool         s_processesGenerated;
    static double       s_default_alpha;
    static double       s_default_scale;
    static int          s_default_flav;
    static int          s_getPDF_default;
    static int          s_fixed;
    static std::vector<double> s_pdfmass;
    
  public :

    Recola_Interface() : ME_Generator_Base("Recola") {  }
    ~Recola_Interface();


    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    static std::string process2Recola(const PHASIC::Process_Info &pi);
    static std::string particle2Recola(const int p);
    static std::string particle2Recola(const std::string p);
    static void increaseProcIndex(){s_recolaProcIndex++;}
    static unsigned int getProcIndex(){return s_recolaProcIndex;}
    static bool checkProcGeneration(){return s_processesGenerated;}
    static void setProcGenerationTrue(){s_processesGenerated = true;}
   
    static int RegisterProcess(const PHASIC::Process_Info& pi,
			       int amptype);
    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res, METOOLS::DivArrD& virt);
    static bool ExitOnError() { return s_exit_on_error; }
    static double GetDefaultAlphaS() { return s_default_alpha; }
    static double GetFixed() { return s_fixed; }
    static double GetDefaultScale() { return s_default_scale; }
    static double GetDefaultFlav() { return s_default_flav; }
    static void SetDefaultFlav(int nlight) {  s_default_flav=nlight; }
    static int GetPDFDefault() { return s_getPDF_default; }
    static std::vector<double> GetPDFMasses() { return s_pdfmass; }


    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add)
    { return NULL; }
    bool NewLibraries() { return false; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }
    static int PDFnf(double scale, int maxn);
    int  PerformTests();
    void PrepareTerminate();

  };


}

#endif
