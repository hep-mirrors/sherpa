#ifndef Recola_Interface_H
#define Recola_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Phys/NLO_Types.H"

#include "recola.h"

namespace PHASIC {
  class External_ME_Args;
}

namespace Recola {

  class Recola_Interface: public PHASIC::ME_Generator_Base {

    static std::string  s_recolaprefix;
    static bool         s_ignore_model;
    static bool         s_exit_on_error;
    static double       s_light_fermion_threshold;
    static size_t       s_recolaProcIndex;
    static bool         s_processesGenerated;
    static size_t       s_default_flav;
    static double       s_default_alphaqcd;
    static double       s_default_scale;
    static size_t       s_fixed_flav;
    static int          s_ewscheme;
    static int          s_amptype;
    static bool         s_compute_poles;
    static bool         s_mass_reg;
    static double       s_photon_mass;
    static bool         s_check_mass;
    static bool         s_self_energy_off;
    static bool         s_use_decay;
    static size_t       s_getPDF_default;

    static std::vector<double> s_pdfmass;
    
    
    static std::map<size_t, bool> s_interference;
    static std::map<size_t,ATOOLS::asscontrib::type> s_asscontribs;
    static size_t       s_doint;
    
    void RegisterDefaults() const;

    
  public :

    Recola_Interface();
    ~Recola_Interface();


    bool Initialize(MODEL::Model_Base *const model,
            BEAM::Beam_Spectra_Handler *const beam,
            PDF::ISR_Handler *const isr,
            YFS::YFS_Handler *const yfs);

    static std::string process2Recola(const ATOOLS::Flavour_Vector& fl);
    static std::string process2Recola(const PHASIC::Process_Info &pi)
    { return process2Recola(pi.ExtractFlavours()); }
    static std::string particle2Recola(const int p);
    static std::string particle2Recola(const std::string p);

    static void increaseProcIndex(){s_recolaProcIndex++;}
    static size_t getProcIndex(){return s_recolaProcIndex;}
    static bool checkProcGeneration(){return s_processesGenerated;}
    static void setProcGenerationTrue(){s_processesGenerated = true;}
   
    // TODO: Clean up
    static size_t RegisterProcess(const PHASIC::Process_Info& pi,
                   int amptype);
    static void   GenerateProcesses(const double& alpha, const double& alphas,
                                    const double& muIR, const double& muUV,
                                    const double& muR2);

    static int RegisterProcess(const PHASIC::External_ME_Args& args,
                   const int& amptype);

    static void EvaluateLoop(int id, const ATOOLS::Vec4D_Vector& momenta, double& bornres, METOOLS::DivArrD& virt, 
                             std::vector<double> &asscontribs);
    static void EvaluateBorn(int id, const ATOOLS::Vec4D_Vector& momenta, double& bornres, int amptype);

    static bool   ExitOnError()              { return s_exit_on_error; }
    static size_t GetFixedFlav()             { return s_fixed_flav; }
    static size_t GetDefaultFlav()           { return s_default_flav; }
    static void   SetDefaultFlav(int nlight) { s_default_flav=nlight; }
    static double GetDefaultAlphaQCD()       { return s_default_alphaqcd; }
    static double GetDefaultScale()          { return s_default_scale; }
    static size_t GetPDFDefault()            { return s_getPDF_default; }
    static int    GetAmpType()               { return s_amptype; }

    static std::vector<double> GetPDFMasses() { return s_pdfmass; }


    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add)
    { return NULL; }
    bool NewLibraries() { return false; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }
    static size_t PDFnf(double scale, size_t maxn);
    int  PerformTests();
    void PrepareTerminate();

    static std::map<size_t,PHASIC::Process_Info> s_procmap;
    static size_t s_vmode;

  };


}

#endif
