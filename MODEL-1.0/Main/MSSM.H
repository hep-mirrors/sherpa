#ifndef MSSM_h
#define MSSM_h

#include "Model_Base.H"

namespace MODEL {
  class MSSM : public Model_Base {
  private :
    std::string m_generator;
    std::string m_benchmark,m_scenario;

    bool        m_scales,m_unification;

    void InitializeBenchmarkPoint() {};
    void ReadInFile();
  public :
    MSSM(std::string,std::string);
    ~MSSM()                     {}
    bool RunSpectrumGenerator();
    void FillMasses()                     {}
    void FillWidths()                     {}
    bool FillDecay(ATOOLS::Decay_Table *) { return 0; }
  };
  /*!
    \class MSSM
    \brief contains the MSSM (no surprise)
  */
  /*!
    \var std::string MSSM::m_generator
    The spectrum generator to be used. At the moment, only Isajet/Isasusy is available.
    \todo Interface more SUSY generators.
  */
  /*!
    \var std::string MSSM::m_benchmark
    This parameters is read in from the data file, default is string(""). with it, it should
    be possible to select one of the SUSY benchmark points and have them initialized without
    further inputs by the spectrum generator of choice. At the moment this option has not been
    implemented.
    \todo Enable direct access to the SUSY benchmark points.
  */
  /*!
    \var std::string MSSM::m_scenario
    In case, no benchmark point has been selected for the SUSY spectrum, a SUSY breaking 
    scenario is read in, the default is mSUGRA.
  */
  /*!
    \var bool MSSM::m_scales
    This flag is set to true, if the SUSY scale is set by the user. 
  */
  /*!
    \var bool MSSM::m_unification
    This flag is set to true, if unification is enforced.
  */
  /*!
    \fn MSSM::MSSM(std::string,std::string)
    The constructor first initalizes the subset of the Standard Model, MODEL::SM. 
    Its parameters (i.e. constants, functions, matrices) are then copied into the newly 
    initialized pointers to the searchable maps. 
  */
  /*!
    \fn bool MSSM::RunSpectrumGenerator()
    In this method the spectrum generator is specified and initialised to run.
    At the moment ,the only option is to run Hdecay.
  */
  /*!
    \fn void MSSM::ReadInFile()
    This method is called by the constructor and selects the model framework. Either a specific
    SUSY benchmark point is chosen (an option that is foreseen only at the moment) or a SUSY
    breaking scenario has to be specified. If the latter option is selected the corresponding
    input parameters are read in. Their default values are:
    - AMSB: \f[M_0 = M_{3/2} = \tan\beta = A_0 = 0,\;\mbox{\rm sign}\mu = 1\f]
    - mGMSB : \f[M_{\rm mes}=2000,\; \Lambda_M = \tan\beta = 0,\;
                 {\cal R}_{\rm grav} = \mbox{\rm sign}\mu = 1,\;
	  	 N_{\rm mes} = 1\f]
    - GMSB : \f[M_{\rm mes}=2000,\; \Lambda_M = \tan\beta = 0,\;
                 {\cal R}_{\rm grav} = \mbox{\rm sign}\mu = 1,\;
	 	 N_{\rm mes} = 1\f]
	     \f[{\cal R}_{\rm gauge} = 1,\; 
	        \Delta_M^{H_u} = \Delta_M^{H_d} = \Delta_M^{Y} = 1,\;
	        N_5(U1) = N_5(SU2) = N_5(SU3) = 1\f]
    - mSUGRA : \f[M_0 = M_{1/2} = \tan\beta = A_0 = 0,\;\mbox{\rm sign}\mu = 1\f]
      (This set of default avalues applies also for the SUGRA breaking
       scenario with enforced unification)
    - SUGRA : \f[M_0 = M_{1/2} = \tan\beta = A_0 = 0,\;\mbox{\rm sign}\mu = 1\f]
              all other parameters are set to \f$10^20$, the SUSY scale is set to 0.

  */
  /*!
    \fn void MSSM::FillMasses()
    \todo This method should be implemented.
  */
  /*!
    \fn void MSSM::FillWidths()
    This method should be implemented. A good way would be to calculate couplings in vertices
    from decay widths to guarantee internal consistency.
  */
}

#endif


