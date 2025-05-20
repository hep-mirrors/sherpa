#ifndef ALPACA_Tools_HeavyIon_Parameters_H
#define ALPACA_Tools_HeavyIon_Parameters_H

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <map>
#include <string>




namespace ALPACA {

  struct xsec_form {
    enum code {
      Gauss     = 2,
      BlackDisc = 1
    };
  };
  std::ostream & operator<<(std::ostream & s, 
			    const xsec_form::code & xsecform);


  typedef std::map<std::string,double> StringDoubleMap;
  typedef StringDoubleMap::iterator    SDiter;

  class HeavyIon_Parameters {
  private:

    StringDoubleMap               m_params;
    xsec_form::code               m_xsecform;


    //Thermal equilibrium parameters (see Event_Generator.H for more detailed description)
    double                m_taumax;
    bool                  m_fixed_sigma;
    double                m_fixed_sigma_val;
    bool                  m_fixed_gamma;
    double                m_fixed_gamma_val;
    bool                  m_only_gluons;
    double                m_p_scale;
    bool                  m_splitting_merging;
    bool                  m_show_event_information;
    bool                  m_show_progress_bar;
    double                m_tau_restart;
    double                m_tsample_min, m_tsample_max;
    bool                  m_include_bose_factors;
    double                m_p_min;
    int                   m_N_include;
    bool                  m_elastic_scattering;
    double                m_tau_max_scaling_limit;
    double                m_alpha_s;
    bool                  m_gaussian_kT2;
    double                m_kT2_reg;
    double                m_f_r;
    double                m_f_delta_p;
    bool                  m_f_shell;
    int                   m_f_N_max;
    double                m_formation_time;
    int                   m_timekeeper;
    double                m_lambda;
    double                m_m2_min_scale;

    double                m_OE_mult_scatter;
    double                m_OE_mult_merge;
    double                m_OE_mult_split;

    double                m_test_double;
    bool                  m_test_bool;


  public:
    HeavyIon_Parameters();
    ~HeavyIon_Parameters();

    void   Init();
    void   RegisterDefaults() const;
    void   ReadParameters();
    double operator()(std::string keyword);
    void   Set(std::string keyword, double value) {
      if (m_params.find(keyword)!=m_params.end()) m_params[keyword] = value;
      else {
	      std::cerr<<"Error in HeavyIon_Parameters::Set: '" << keyword << "' not found.\n";
      }
    }
    xsec_form::code             XSec_Form()             const { return m_xsecform; }
    double                      TauMax()                const { return m_taumax; }
    std::pair<bool, double>     FixedSigma()            const { return std::make_pair(m_fixed_sigma, m_fixed_sigma_val); }
    std::pair<bool, double>     FixedGamma()            const { return std::make_pair(m_fixed_gamma, m_fixed_gamma_val); }
    bool                        OnlyGluons()            const { return m_only_gluons; }
    bool                        SplittingMerging()      const { return m_splitting_merging; }
    bool                        ShowEventInformation()  const { return m_show_event_information; }
    bool                        ShowProgressBar()       const { return m_show_progress_bar; }
    double                      TauRestart()            const { return m_tau_restart; }
    double                      tsampleMin()            const { return m_tsample_min; }
    double                      tsampleMax()            const { return m_tsample_max; }
    bool                        IncludeBoseFactors()    const { return m_include_bose_factors; }
    double                      pMin()                  const { return m_p_min; }
    int                         NInclude()              const { return m_N_include; }
    bool                        ElasticScattering()     const { return m_elastic_scattering; }
    double                      TauMaxScalingLimit()    const { return m_tau_max_scaling_limit; }
    double                      AlphaS()                const { return m_alpha_s; }
    bool                        GaussiankT2()           const { return m_gaussian_kT2; }
    double                      kT2Reg()                const { return m_kT2_reg; }
    double                      fr()                    const { return m_f_r; }
    double                      fDeltap()               const { return m_f_delta_p; }
    bool                        fShell()                const { return m_f_shell; }
    int                         fNMax()                 const { return m_f_N_max; }
    double                      FormationTime()         const { return m_formation_time; }
    int                         Timekeeper()            const { return m_timekeeper; }
    double                      Lambda()                const { return m_lambda; }
    double                      M2MinScale()            const { return m_m2_min_scale; }

    double                      OEMultScatter()         const { return m_OE_mult_scatter; }
    double                      OEMultMerge()           const { return m_OE_mult_merge; }
    double                      OEMultSplit()           const { return m_OE_mult_split; }

    double                      TestDouble()           const { return m_test_double; }
    bool                        TestBool()             const { return m_test_bool; }
    
  };
  
  extern HeavyIon_Parameters HIPars;
}

#endif
