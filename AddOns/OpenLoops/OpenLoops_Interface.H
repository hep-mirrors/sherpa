#ifndef OpenLoops_Interface_H
#define OpenLoops_Interface_H

#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include <dirent.h>


namespace OpenLoops {

  class OpenLoops_Particle {
  public:
    std::string m_name, m_faname;
    size_t m_order;
    OpenLoops_Particle() : m_name(""), m_faname(""), m_order(0) {}
    OpenLoops_Particle(std::string name, std::string faname, size_t order) :
      m_name(name), m_faname(faname), m_order(order) {}
  };


  class OpenLoops_Interface: public PHASIC::ME_Generator_Base {

    void InitializeParticles();

    static std::string MatchOptions(std::vector<std::string> options,
                                    int oew, int oqcd, int nloop);
    static int SelectInfo(const dirent *entry);
    static bool Order(const std::pair<size_t, ATOOLS::Flavour> &a,
                      const std::pair<size_t, ATOOLS::Flavour> &b);

  public :
    OpenLoops_Interface() : ME_Generator_Base("OpenLoops") {}
    ~OpenLoops_Interface() {}

    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr);

    void OpenLoopsInit(double mur2, double alpha_QED, double alpha_S);

    static std::string GetProcessPermutation(const ATOOLS::Flavour_Vector& flavs_orig,
                                             std::vector<int>& permutation);
    static std::pair<std::string, std::string> ScanFiles(std::string& process,
                                                         int oew, int oqcd, int nloop);
    static ATOOLS::Flavour_Vector MapFlavours(const ATOOLS::Flavour_Vector& orig);

    PHASIC::Process_Base *InitializeProcess(const PHASIC::Process_Info &pi, bool add)
    { return NULL; }
    int  PerformTests() { return 1; }
    bool NewLibraries() { return false; }
    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}
    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode,const double &kt2) { return NULL; }

    static std::string s_olprefix;
    static std::vector<std::string> s_allowed_libs;
    static int s_pole_mode, s_amp_switch, s_amp_switch_rescue, s_nf;
    static bool s_generate_list;
    static std::set<std::string> s_shoppinglist;
    static std::map<int, OpenLoops_Particle> s_particles;
  };

}



typedef void (*Amp2Func)(double* moms,
                         double* B,
                         double* V_finite, double* V_eps, double* V_eps2,
                         double* I_finite, double* I_eps, double* I_eps2);

typedef void (*PermutationFunc)(int* permutation);

extern "C" {

  void parameters_init_(double* Mass_E, double* Mass_M, double* Mass_L,
                        double* Mass_U, double* Mass_D, double* Mass_S,
                        double* Mass_C, double* Width_C, double* Mass_B, double* Width_B, double* Mass_T, double* Width_T,
                        double* Mass_W, double* Width_W, double* Mass_Z, double* Width_Z, double* Mass_H, double* Width_H,
                        double* Coupl_Alpha_QED, double* Coupl_Alpha_QCD,
                        int* last_switch, int* amp_switch, int* amp_switch_rescue, int* use_coli_cache, int* check_Ward_tree, int* check_Ward_loop, int* out_symmetry);


  void loop_parameters_init_(double* renscale, double* fact_UV, double* fact_IR, double* pole1_UV, double* pole1_IR, double* pole2_IR,
                             int* polenorm_swi, int* N_quarks, int* nq_nondecoupled,
                             double* opp_rootsvalue, double* opp_limitvalue, double* opp_thrs, int* opp_idig, int* opp_scaloop,
                             int* sam_isca, int* sam_verbosity, int* sam_itest, int* fermion_loops, int* nonfermion_loops,
                             int* CT_on, int* R2_on, int* IR_on, int* polecheck, double* set_C_PV_threshold, double* set_D_PV_threshold, int* set_DD_red_mode);
  void openloops_welcome_(char* str, int* len_str);
  
}

#endif
