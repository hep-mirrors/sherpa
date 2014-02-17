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

    static void SetParameter(const std::string & key, double value);
    static void SetParameter(const std::string & key, int value);
    static void SetParameter(const std::string & key, std::string value);
    static double GetDoubleParameter(const std::string & key);
    static int GetIntParameter(const std::string & key);
    static void FlushParameters();

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
    (PHASIC::Process_Base *const proc,const size_t &mode) { return NULL; }

    static std::string s_olprefix;
    static std::vector<std::string> s_allowed_libs;
    static bool s_generate_list;
    static std::set<std::string> s_shoppinglist;
    static std::map<int, OpenLoops_Particle> s_particles;
  };


}



typedef void (*Amp2Func)(double* moms, double* M2L0, double* M2L1, double* IRL1,
                         double* M2L2, double* IRL2);

typedef void (*PermutationFunc)(int* permutation);

#endif
