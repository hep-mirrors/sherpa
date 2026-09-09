#ifndef Recola_Interface_H
#define Recola_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "ATOOLS/Phys/NLO_Types.H"

namespace PHASIC {
  struct External_ME_Args;
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
    static bool         s_use_decay;
    static size_t       s_getPDF_default;

    static std::vector<double> s_pdfmass;
    
    
    /*!
      Beam polarisation / parity-violating asymmetry, read from the same
      POLARIZATIONS block AMEGIC uses (Pol_Info.C). When PV_ASYMMETRY is on,
      EvaluatePolarised() sums Recola's helicity-resolved squared amplitudes
      into  sum_hel lambda_1 |M|^2 / sum_hel |M|^2  at LO and at NLO. The NLO
      one is the point: it is the genuine O(alpha) corrected asymmetry, which
      the Born-only AMEGIC column cannot produce.
    */
    static bool         s_polarised;
    static double       s_pol1, s_pol2;
    static double       s_pvasym_lo, s_pvasym_nlo;
    static bool         s_pvasym_valid;
    //! Polarised/unpolarised cross-section ratio for the declared beam
    //! polarisations, sum_hel (1+l1 P1)(1+l2 P2)|M|^2 / sum_hel |M|^2.
    static double       s_polwgt_lo, s_polwgt_nlo;
    //! Refuse to enumerate more helicity configurations than this.
    static size_t       s_pol_maxcfg;

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

    //! Is POLARIZATIONS: PV_ASYMMETRY switched on?
    static bool   PolarisedOn()   { return s_polarised; }
    //! Born and one-loop-corrected PV asymmetry ratios from the last evaluation.
    static bool   PVAsymValid()   { return s_pvasym_valid; }
    static double PVAsymLO()      { return s_pvasym_lo; }
    static double PVAsymNLO()     { return s_pvasym_nlo; }
    /*!
      Loop over all +-1 helicity assignments of the external legs, pull the
      polarised squared amplitudes out of Recola and build the asymmetry.
      Must be called after compute_process_rcl for process id.
    */
    static double PolWeightLO()   { return s_polwgt_lo; }
    static double PolWeightNLO()  { return s_polwgt_nlo; }
    /*!
      with_nlo=false when the process was computed at 'LO' only: asking Recola
      for an 'NLO' polarised amplitude it has not computed is invalid.
    */
    static void   EvaluatePolarised(int id, int boqcd, int voqcd,
                                    const ATOOLS::Flavour_Vector &flavs,
                                    bool with_nlo = true);

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

    static void SetAlphas(double alphas, double scale, int nflavour);

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
