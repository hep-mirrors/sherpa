#ifndef Run_Parameter_h
#define Run_Parameter_h

#include <string>
#include "Data_Read.H"
#include "Flavour.H"
#include "Matrix.H"

#define DEVELOP
#ifndef DEVELOP
// don't check, do ordinary return
#define SAVE_RETURN(x,name) return x;
#else
// check if sucessfully filled everytime it is used.
// "IsDefined" must be declared for all types!
#include "Data_Return.H"
#define SAVE_RETURN(x,name) return IsDefined(x,name);
#endif


namespace AORGTOOLS {  
  //! An all purpose class for reading in parameters.
  class Run_Parameter {
    std::string         path;

  public:
    Run_Parameter()    {gen.output = 0;}
    ~Run_Parameter();
    void Init(std::string);
    
    // CONSTs!

    class {  // general parameters
      int                output;
      std::string        runtime;
      long               nevents;
      double             cms_energy;
      double             accu;
      Switch::code       massive;      
      APHYTOOLS::Flavour beam1, beam2;
      std::string        rpa_id;
      std::string        flavour_id;
    public:
      inline long NumberOfEvents()             { SAVE_RETURN(nevents,"nevents");}
      inline long SetNumberOfEvents(long N)    { return nevents=N;}
      inline double Ecms()                     { SAVE_RETURN(cms_energy,"cms_energy");}
      inline double SetEcms(double v)          { return cms_energy=v;}
      // output switches
      inline int Output()                      { SAVE_RETURN(output,"output"); }
      inline int Error()                       { return (output>=0); }
      inline int Events()                      { return (output>=1); }
      inline int Tracking()                    { return (output>=2);}
      inline int Debugging()                   { return (output>=3);}
      inline int SetOutput(int o)              { return output=o;}
      inline void Output_Off()                 { output = -1;}
      // misc
      inline double Accu()                     { SAVE_RETURN(accu,"accu");}
      // beam particles
      inline const APHYTOOLS::Flavour Beam1()  { SAVE_RETURN(beam1,"beam1");}
      inline const APHYTOOLS::Flavour Beam2()  { SAVE_RETURN(beam2,"beam2");}

      inline Switch::code Masses()             { SAVE_RETURN(massive,"massive");}
      friend class Run_Parameter;
      // for the tests.
      inline void SetBeam1(const APHYTOOLS::Flavour _b) { beam1 = _b; }
      inline void SetBeam2(const APHYTOOLS::Flavour _b) { beam2 = _b; }

      inline std::string GetRpaID()        { return rpa_id; }
      inline std::string GetFlavourID()    { return flavour_id; }

      inline void  SetRpaID(const std::string & _id)        { rpa_id = _id; }
      inline void  SetFlavourID(const std::string & _id)    { flavour_id = _id; }
      
    } gen;

    class { // matrix element
      std::string       model_file;
      Model_Type::code  the_model;
      String_Type::code m_string; 
      Switch::code      model_mass;
      Switch::code      coulomb_corr;
      int               cut_scheme;
      int               kfactorscheme;
      int               scalescheme;
      AMATOOLS::Matrix<4>        Theta; // NC QED
    public:
      inline std::string ModelFile()         { SAVE_RETURN(model_file,"model_file");}
      inline Model_Type::code Model()        { SAVE_RETURN(the_model,"the_model");}
      inline Switch::code UsingModelMass()   { SAVE_RETURN(model_mass,"model_mass");}
      inline Switch::code UsingCoulombCorr() { SAVE_RETURN(coulomb_corr,"coulomb_corr");}
      inline String_Type::code String()      { SAVE_RETURN(m_string,"m_string");}
      inline void SetString(const String_Type::code& stype) { m_string = stype;}
      inline int CutScheme()                { return cut_scheme; }
      inline void  SetCutScheme(int _cs )   { cut_scheme=_cs; }
      inline int    KFactorScheme()         { SAVE_RETURN(kfactorscheme,"kfactorscheme"); }
      inline int    ScaleScheme()           { SAVE_RETURN(scalescheme,"scalescheme"); }
      AMATOOLS::Matrix<4>* GetTheta()       { return &Theta;}
      friend class Run_Parameter;
    } me;

    class { // integration
      int    integrator;
      int    jetfinder;
      double ycut;
      double accu;
    public:
      inline double Ycut()    { return ycut; }
      inline double Accu()   { return accu; }
      inline int Integrator() { return integrator; }
      inline int Jetfinder()  { return jetfinder; } 
      inline void Set_Ycut(double _ycut) { ycut = _ycut; }
      // inline double Ycut()   { SAVE_RETURN(ycut,"ycut");}
      // inline int JetFinder() { SAVE_RETURN(jet_scheme,"jet_scheme");}
      friend class Run_Parameter;
    } integ;


    class { // parton shower
      double           final_q02;
      double           initial_q02;   // t0 for initial state shower
      double           lambda_qcd;
    public:
      inline double FinalQ02()          { SAVE_RETURN(final_q02,"final_q02"); }
      inline double InitialQ02()        { SAVE_RETURN(initial_q02,"initial_q02");}
      inline double Lambda()            { SAVE_RETURN(lambda_qcd,"lambda_qcd");}
      inline double SetLambda(double l) { return lambda_qcd=l;}
      // min transversal momentum squared
      friend class Run_Parameter;
    } pshower ;

    class { // dipole shower
      int              k_perp_scheme;
      int              k_ord_scheme;
      double           kt_min;

    public:
      inline int KtScheme()            { SAVE_RETURN(k_perp_scheme,"k_perp_scheme");}
      inline int SetKtScheme(int k)    { return k_perp_scheme=k;}
      inline int KOrd()                { SAVE_RETURN(k_ord_scheme,"k_ord_scheme");}
      inline int SetKOrd(int k)        { return k_ord_scheme=k;}
      inline double KtMin()            { SAVE_RETURN(kt_min,"kt_min");}
      inline double SetKtMin(double k) { return kt_min=k;}
      friend class Run_Parameter;
    } dshower ;

    class cl_consts { // general coupling consts (possibly running)
      double           alpha_s_fix;

      Switch::code     running_masses;
      int              running_widths;
      Switch::code     running_aqed;

      Switch::code     running_alphas;
      double           aqed_mz;
      double           sin2_tw;
    public:
      inline Switch::code IsASRunning()      { SAVE_RETURN(running_alphas,"running_alphas");}
      inline Switch::code IsAQEDRunning()    { SAVE_RETURN(running_aqed,"running_aqed");}
      inline double AqedMZ()                 { SAVE_RETURN(aqed_mz,"aqed_mz");}
      inline double Sin2TW()                 { SAVE_RETURN(sin2_tw,"sin2_tw");}
      inline double FixedAlphaS()            { SAVE_RETURN(alpha_s_fix,"alpha_s_fix"); }
      inline double SetFixedAlphaS(double a) { return alpha_s_fix=a;}     
      inline Switch::code IsMassRunning()    { SAVE_RETURN(running_masses,"running_masses");}
      inline int IsWidthRunning()            { SAVE_RETURN(running_widths,"running_widths");}

      double Mass(APHYTOOLS::Flavour,double);
      double RunningMass(APHYTOOLS::Flavour,double);
      double Width(APHYTOOLS::Flavour,double);
      friend class Run_Parameter;
    } consts;

    class { // test == only for the experts
      int    analysis;
    public:
      inline int Analysis()          { SAVE_RETURN(analysis,"analysis"); }
      inline int SetAnalysis(int k)  { return analysis=k;}
      friend class Run_Parameter;
    } test;


    double Picobarn()                       { return 3.89379656e8;} // pb/GeV^2
    double c()                              { return 299792458.e3;} // mm/s
    double hbar()                           { return 6.58211889e-25;} // GeV s
    std::string GetPath()                   { return path;} 
    std::string SetPath(std::string _path)  { path = _path;} 
  };

  extern Run_Parameter rpa;
}


#endif  // Run_Parameter_h
