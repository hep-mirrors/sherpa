#ifndef Recola_Interface_H
#define Recola_Interface_H

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"

#include "recola.hpp"

namespace PHASIC {
  class External_ME_Args;
}

namespace Recola {

  class Recola_Interface: public PHASIC::ME_Generator_Base {

    static std::string  s_recolaprefix;
    static double       s_light_fermion_threshold;
    static unsigned int s_recolaProcIndex;
    static bool         s_processesGenerated;
    static int          s_default_flav;
    static double       s_default_alphaqcd;
    static double       s_default_scale;
    static int          s_fixed;
    static int          s_ewscheme;

    static std::vector<double> s_pdfmass;
    
  public :

    Recola_Interface() : ME_Generator_Base("Recola") {  }
    ~Recola_Interface();


    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    static std::string process2Recola(const ATOOLS::Flavour_Vector& fl);
    static std::string process2Recola(const PHASIC::Process_Info &pi)
    { return process2Recola(pi.ExtractFlavours()); }
    static std::string particle2Recola(const int p);
    static std::string particle2Recola(const std::string p);
    static void increaseProcIndex(){s_recolaProcIndex++;}
    static unsigned int getProcIndex(){return s_recolaProcIndex;}
    static bool checkProcGeneration(){return s_processesGenerated;}
    static void setProcGenerationTrue(){s_processesGenerated = true;}
   
    // TODO: Clean up
    static int RegisterProcess(const PHASIC::Process_Info& pi,
			       int amptype);

    static int RegisterBorn(const PHASIC::External_ME_Args& args,
			    const int& amptype);

    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& momenta,
                             double& res, METOOLS::DivArrD& virt);
    static int GetEWScheme() { return s_ewscheme; }
    static double GetFixed() { return s_fixed; }
    static int GetDefaultFlav() { return s_default_flav; }
    static double GetDefaultAlphaQCD() {  return s_default_alphaqcd; }
    static double GetDefaultScale() {  return s_default_scale; }
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
