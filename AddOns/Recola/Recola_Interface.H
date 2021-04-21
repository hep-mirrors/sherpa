#ifndef Recola_Interface_H
#define Recola_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"

#include "recola.h"

namespace Recola {

  struct amptype {
    enum code {
      none = 0,
      treetree = 1,
      treeloop = 2,
      looploop = 3
    };
  };

  std::ostream & operator<<(std::ostream & s,const amptype::code &at);

  class Recola_Interface: public PHASIC::ME_Generator_Base,
                          public ATOOLS::Terminator_Object {

    static std::string  s_recolaprefix;
    static bool         s_ignore_model;
    static bool         s_exit_on_error;
    static bool         s_use_iop_in_ewapprox;
    static double       s_light_fermion_threshold;
    static size_t       s_recolaProcIndex;
    static bool         s_processesGenerated;
    static double       s_default_alphaqcd;
    static double       s_default_scale;
    static int          s_default_flav;
    static int          s_getPDF_default;
    static int          s_fixed_flav;
    static double       s_ir_scale;
    static double       s_uv_scale;
    static int          s_collier_cache;
    static size_t       s_doint;
    static std::vector<double> s_pdfmass;

    static std::map<size_t,PHASIC::Process_Info> s_procmap;
    static std::map<size_t,PHASIC::asscontrib::type> s_asscontribs;

  public :

    Recola_Interface() : ME_Generator_Base("Recola") {  }
    ~Recola_Interface() {}


    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    static std::string process2Recola(const PHASIC::Process_Info &pi);
    static std::string particle2Recola(const int p);
    static std::string particle2Recola(const std::string p);

    static void increaseProcIndex()     { s_recolaProcIndex++; }
    static size_t getProcIndex()        { return s_recolaProcIndex; }
    static bool checkProcGeneration()   { return s_processesGenerated; }
    static void setProcGenerationTrue() { s_processesGenerated = true; }

    static int    RegisterProcess(const PHASIC::Process_Info& pi,
                                  const amptype::code& at);
    static void   EvaluateProcess(int id, const Vec4D_Vector &momenta,
                                  const size_t &voqcd, const size_t &boqcd,
                                  METOOLS::DivArrD &Vqcd, double &B,
                                  std::vector<double> &asscontribs,
                                  const bool& calcass);

    static bool   ExitOnError()        { return s_exit_on_error; }
    static bool   UseIopInEWapprox()   { return s_use_iop_in_ewapprox; }
    static double GetDefaultAlphaQCD() { return s_default_alphaqcd; }
    static double GetDefaultScale()    { return s_default_scale; }
    static double GetDefaultFlav()     { return s_default_flav; }
    static double GetFixedFlav()       { return s_fixed_flav; }
    static void   SetDefaultFlav(int nlight) {  s_default_flav=nlight; }
    static int    GetPDFDefault()      { return s_getPDF_default; }
    static std::vector<double> GetPDFMasses() { return s_pdfmass; }

    static double IRScale() { return s_ir_scale; }
    static double UVScale() { return s_uv_scale; }

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
