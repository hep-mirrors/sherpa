#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"
#include "Tools.H"

namespace HADRONS {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Lepton : public HD_ME_Base {
    int m_nutau, m_nulep, m_lep;
    double m_GF2;                           // Fermi constant
    Complex m_cR1, m_cL1, m_cR2, m_cL2;     // couplings for Z function
    public:
      Tau_Lepton( int nout, ATOOLS::Flavour *fl );
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p ); 
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pseudo : public HD_ME_Base {
    int     m_nutau, m_pion;
    bool    m_pionmode;             // either tau->pi or tau->kaon
    double  m_fxx2;                 // decay constants (f_pi, f_K)
    double  m_Vxx2;                 // CKM ME (Vud, Vus)
    double  m_GF2;                  // Fermi constant
    Complex m_cR, m_cL;             // coupling for X function
    public:
      Tau_Pseudo( int nout, ATOOLS::Flavour *fl );
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p );
      inline void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Two_Pion : public HD_ME_Base {
    int    m_nutau, m_pion_ch, m_pion0;
    bool   m_pionmode;              // either tau->2pi or tau->2kaon
    int    m_running, m_ff;     // running width, ff model
    double m_Vud2;                  // CKM ME
    double m_fxx;                   // decay constant
    double m_GF2;                   // Fermi constant
    double m_frho, m_grpp;          // rho-pipi coupling
    double m_CG;                    // Clebsch-Gordon-coefficient
    Complex m_cR, m_cL;             // coupling for X function
    double m_MR, m_GR, m_MR2;
    double m_beta, m_lambda, m_delta, m_gamma;
    double m_MRR, m_GRR, m_MRRR, m_GRRR, m_MRR2, m_MRRR2;
    double m_gammaR, m_gammaRR, m_gammaRRR;
    double m_m, m_m2;               // mass of pion or kaon
    double m_m2_pi, m_m2_K;         // mass pf pion, kaon
    public:
      Tau_Two_Pion( int nout, ATOOLS::Flavour *fl );
      double operator()(const ATOOLS::Vec4D *_p );
      Complex A( double x, double y );
      Complex FormFactor( double s ); 
      double Using_Hels( const ATOOLS::Vec4D *_p );
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pion_Kaon : public HD_ME_Base {
     
    class FF_Base {
      public :
      FF_Base() { };
      virtual Complex VectorFormFactor( double ) = 0;
      virtual Complex ScalarFormFactor( double ) = 0;
      virtual void SetModelParameters( GeneralModel ) = 0;
      virtual void SetMasses2( double, double, double ) {};
    };
      

    class RChT : public FF_Base {       // Guerrero, Tesis Doctoral
       
      double m_MK2, m_GK, m_fpi2;       // K* resonance par's, pi decay constant
      double m_MK02, m_GK0;             // K0* resonance par's
      double m_renorm2;                 // renormalisation scale
      double m_mPi2, m_mK2, m_mEta2, m_mPi, m_mK, m_mEta;  // mass of pi, K, eta
      double m_Sigma_KP, m_Delta_KP;    // sum, difference of masses^2
      
      Complex JBar( double s, double MP2, double MQ2, double Sigma, double Delta );
      Complex JBarBar( double s, double MP2, double MQ2, double Sigma, double Delta );
      Complex Mr( double s, double MP2, double MQ2 );
      Complex L( double s, double MP2, double MQ2 );
      double MassWidthVector( double s );
      double MassWidthScalar( double s );
       
      public:
      RChT() {};
      ~RChT() {};
      Complex VectorFormFactor( double );
      Complex ScalarFormFactor( double );
      void SetModelParameters( GeneralModel _md );
      void SetMasses2( double, double, double );
    }; 

    class KS : public FF_Base {     // Guerrero, Tesis Doctoral
       
      int *  m_part;                // internal numbers for particles
      double m_MR, m_MRR;
      double m_MR2, m_MRR2;
      double m_GR, m_GRR;
      double m_m2, m_mK2;           // pion, kaon mass
      double m_beta;
      int    m_running;
      public:
      KS() {};
      ~KS() {};
      void SetModelParameters( GeneralModel _md );
      Complex VectorFormFactor( double );
      Complex ScalarFormFactor( double );
    }; 

    int      m_nutau, m_pion, m_kaon;
    bool     m_chpionmode;      // pion = pi+ or pi0 
    double   m_Vus2;            // CKM ME
    double   m_fpi2;            // decay constant
    double   m_GF2;             // Fermi constant
    Complex  m_cR, m_cL;        // coupling for X function
    double   m_ms[4];           // masses^2 
    double   m_Delta_KP;        // difference of masses^2

    FF_Base * p_ff;
     
    public:
      Tau_Pion_Kaon( int nout, ATOOLS::Flavour *fl );
      ~Tau_Pion_Kaon() { if (p_ff!=NULL) delete p_ff; }
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p );
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Three_Pseudo : public HD_ME_Base {
     
    class FF_Base {
      double m_ms[3];                       // masses of particles^2
      int  * m_part;
      public :
      FF_Base() { };
      virtual void SetMode( int ) {};
      virtual void SetPath( std::string ) = 0;
      virtual Complex FormFactor( int, double, double, double ) = 0;
      virtual void SetModelParameters( GeneralModel ) = 0;
      inline void SetInternalNumbers( int * part ) { m_part = part; }
      inline void SetOutgoingMasses2( double ms1, double ms2, double ms3 ) {
        m_ms[0] = ms1;                      // mass of first pseudo
        m_ms[1] = ms2;                      // mass of second pseudo
        m_ms[2] = ms3;                      // mass of third pseudo
      }
      inline double Mass2(int i)    { return m_ms[i]; }
      inline int ParticleNo(int i)  { return m_part[i]; }
    };
      

    class RChT : public FF_Base {       // Dumm, Pich, Portoles, hep-ph/0312183 
      inline double Sqrt_Lambda( double _a, double _b, double _c ) {
        return sqrt(Tools::Lambda(_a,_b,_c)); }
       
      double  MassWidthVector( double s ); 
      double  MassWidthVector( int a, double s ); 
      double  MassWidthAxial( double Q2 ); 
       
      ATOOLS::Histogram *    CreatePhiHistogram();     // pre-calculation of phi(q2)
      double  Phi( double Q2 ); 
      double  IntegralPhi( double );
      double  FFunc( double, double, double );
       
      double m_l0, m_exp_alpha;         // fitting parameter
      double m_MA, m_MA2;               // axial resonance mass
      double m_msV[2];                  // masses of vector resonances
      double m_widthV[2];                 
      double m_GA_at_MA2;               // on-shell axial width
      double m_Phi_at_MA2;              // Phi(MA2)
      double m_m, m_m2, m_mK2;          // masses of pion and kaon
      double m_fpi2;                    // decay constant
      double m_gammaR;                  // global constant for GV
      ATOOLS::Histogram * p_phi;        // pointer on phi histogram
      std::string m_path;               // path of Decaydata files

      int    m_mode;                    // internal mode representation
      bool   m_twoident;                // true if two vector resonances are identical
      int    m_running;                 // running width

      double m_GV, m_FV, m_FV2;         // couplings
      double m_FA, m_FA2;
      double m_lsum;                    // fit parameters
      double m_l1, m_l2;
       
      public:
      RChT() : p_phi (NULL) { };
      ~RChT() { if (p_phi!=NULL) delete p_phi; }
      Complex FormFactor( int, double, double, double ); 
      void SetModelParameters( GeneralModel _md );
      void SetPath( std::string s ) { m_path = s; }
      void SetMode( int );
    }; 

    class KS : public FF_Base {     // Decker, Finkemeier, Mirkes, hep-ph/9310270
      // order of Pseudos: see table 1 of the paper
       
      inline Complex BW_V( int a, int b, double s ) {
        if (m_running & 2 ) {
          return Tools::BreitWigner( s, m_msV[a],  m_widthV[a], Mass2(a), Mass2(b), 1. );
        }
        else {
          return Tools::BreitWigner( s, m_msV[a], sqrt(m_msV[a])*m_widthV[a] );
        }
      }
      
      inline Complex BW_v( int a, int b, double s ) {
        if (m_running & 1 )
          return Tools::BreitWigner( s, m_msv[a],  m_widthv[a], Mass2(a), Mass2(b), 1. );
        else
          return Tools::BreitWigner( s, m_msv[a], sqrt(m_msv[a])*m_widthv[a] );
      }
       
      inline double Sqrt_Lambda( double _a, double _b, double _c ) { 
        return sqrt(Tools::Lambda(_a,_b,_c)); }
         
      Complex Trho( double s );
      Complex TKstar( double s );
      Complex BW_A( double s );
      Complex Tvector1( int a, int b, double x );
      Complex Tvector2( int a, int b, double x );
      Complex TSvector( int a, int b, int c, double Q2, double s, double t );
      Complex Tgen( int a, int b, int c, double s, double t);
       
      ATOOLS::Histogram * CreateGHistogram();       // pre-calculation of G(q2)
      double  G( double s );
      double  IntegralG( double );
       
      double m_MA2, m_MAA2;                 // axial resonances
      double m_MA, m_MAA;                   // axial resonances
      double m_GA, m_GAA;
      double m_alpha, m_Beta[2];            // relative strength of resonances
      double m_msV[2], m_msv[2];            // masses of vector resonances
      double m_widthV[2], m_widthv[2];
      int    m_running; 
      double m_exp_alpha;                   // exponent in GA
      double m_G_at_MA2, m_G_at_MAA2;
      ATOOLS::Histogram * p_G;      // pointer on G histogram
      std::string m_path;           // path of Decaydata files

      int    m_mode;                // internal mode representation
      bool   m_twoident;            // true if two vector resonances are identical
      bool   m_G123;                // false if there is only one amplitude
      double m_X123, m_ms123;
       
      public:
      KS() : p_G (NULL) { } 
      ~KS() { if (p_G!=NULL) delete p_G; }
      Complex FormFactor( int, double, double, double );
      void SetMode( int );
      void SetModelParameters( GeneralModel _md );
      void SetPath( std::string s ) { m_path = s; }
    }; 

    int    m_nutau, m_pseudo_3, m_pseudo_1, m_pseudo_2;
    int    m_part[4];
    double m_Vud, m_Vus;        // CKM ME
    double m_fpi2;              // decay constant
    double m_GF2;               // Fermi constant
    Complex m_cR, m_cL;         // X function couplings
    double m_ms[5];             // masses^2
    int    m_mode;              // code for decay mode
    double m_A123;              // coupling constant
    double m_global;            // global constant
    void   SetA123();
    FF_Base * p_ff;
    public:
       
      Tau_Three_Pseudo( int nout, ATOOLS::Flavour *fl );
      ~Tau_Three_Pseudo() {
        if (p_ff!=NULL) delete p_ff;
      }
     
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p );
      Complex FormFactor( int j, double Q2, double s, double t ); 
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Four_Pion_3 : public HD_ME_Base {

    class LorenzBase {
      public :
        Complex         * p_X;
       
        ATOOLS::Vec4D   * p_p;
        ATOOLS::Vec4D     m_r[5];
        double            m_s[5];
        double            m_q2;
       
      LorenzBase( ) {};
      virtual Complex   operator()( int ) {};
      virtual void      SetModelParameters( GeneralModel ) {};
      void              SetPrivates( Complex *, ATOOLS::Vec4D * );
    };

    class KS : public LorenzBase {  // CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
            // this class does not only provide form factors
            // but the complete Lorentz structure
       

      double  m_fpi2;                   // pion decay constant
      double  m_Go3p, m_Frho, m_grop;   // coupling constants
      double  m_mpi2, m_mpi02;          // mass^2 of pion, pion0
      double  m_MR, m_MRR, m_MRRR;      // mass of rho resonances
      double  m_MR2, m_MRR2, m_MRRR2;   // mass^2 of rho resonances
      double  m_GR, m_GRR, m_GRRR;      // width of rho resonances
      double  m_MO, m_MO2, m_GO;        // omega resonance
      double  m_MS, m_MS2, m_GS;        // sigm0 resonance
      double  m_MF, m_MF2, m_GF;        // f0 resonance
      double  m_MA, m_MA2, m_GA;        // axial resonance
      double  m_beta, m_gamma;          // relative strength of vector resonances
      double  m_sigma;                  // relative strength of vector resonance
      double  m_R[5];                   // important factor
      Complex OmegaPi(),  m_Beta_opi[4];
      Complex AonePi(),   m_Beta_api[4];
      Complex SigmaRho(), m_Beta_srh[4];
      Complex FzeroRho(), m_Beta_frh[4];
       
      double  Dots( int, int );
      Complex Fk( double, Complex * );
      Complex Trho( double );
      Complex TTrho( double );
       
      public : 
       
      KS() : LorenzBase() {};
      Complex operator()( int );
      void    SetModelParameters( GeneralModel _md );
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
    double       m_ms[6];
    int          m_inter[6];            // internal numbering
     
    double       m_Vud2,m_Vus2;         // CKM ME
    double       m_GF2;                 // Fermi constant
    Complex      m_cR, m_cL;                // couplings in X functions
    LorenzBase * p_lorenz;              // pointer to current contributions
    Complex      m_Alpha[m_ncontrib];   // weight of each contribution
    Complex      m_SumAlpha;            // sum of all weights
     
    ATOOLS::Vec4D   m_p[6];             // momenta
    Complex         m_X[5];             // value of an X function
     
    public:
     
    Tau_Four_Pion_3( int nout, ATOOLS::Flavour *fl );
    ~Tau_Four_Pion_3() {
      if (p_lorenz!=NULL) delete p_lorenz;
    }
   
    double operator()(const ATOOLS::Vec4D *_p );
    double Using_Hels( const ATOOLS::Vec4D *_p );
    void   SetModelParameters( GeneralModel _md );

  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Four_Pion_1 : public HD_ME_Base {

    class LorenzBase {
      public :
      Complex * p_X;
       
        ATOOLS::Vec4D   * p_p;
        ATOOLS::Vec4D     m_r[5];
        double            m_s[5];
        double            m_q2;
       
      LorenzBase( ) {};
      virtual Complex operator()() {};
      virtual void SetModelParameters( GeneralModel ) {};
      void SetPrivates( Complex *, ATOOLS::Vec4D * );
    };

    class KS : public LorenzBase {  // CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
       

      double  m_fpi2;                   // pion decay constant
      double  m_mpi2, m_mpi02;          // mass^2 of pion, pion0
      double  m_MR, m_MRR, m_MRRR;      // mass of rho resonances
      double  m_MR2, m_MRR2, m_MRRR2;       // mass^2 of rho resonances
      double  m_GR, m_GRR, m_GRRR;      // width of rho resonances
      double  m_beta, m_gamma;          // relative strength of vector resonances
       
      double  Dots( int, int );
      Complex Trho( double );
       
      public : 
       
      KS() : LorenzBase() {};
      Complex operator()();
      void SetModelParameters( GeneralModel _md );
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
    double       m_ms[6];
    int          m_inter[6];            // internal numbering
     
    double       m_Vud2,m_Vus2;             // CKM ME
    double       m_GF2;                 // Fermi constant
    Complex      m_cR, m_cL;                // couplings in X functions
    LorenzBase * p_lorenz;              // pointer to current contributions
     
    ATOOLS::Vec4D   m_p[6];             // momenta
    Complex         m_X[5];             // value of an X function
     
    public:
     
    Tau_Four_Pion_1( int nout, ATOOLS::Flavour *fl );
    ~Tau_Four_Pion_1() {
      if (p_lorenz!=NULL) delete p_lorenz;
    }
   
    double operator()(const ATOOLS::Vec4D *_p );
    double Using_Hels( const ATOOLS::Vec4D *_p );
    void SetModelParameters( GeneralModel _md );

  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Eta_Two_Pion : public HD_ME_Base {
     
    int    m_nutau, m_eta, m_pion, m_pion0;
    double m_Vud2;                  // CKM ME
    double m_fpi2;                  // decay constant
    double m_GF2;                   // Fermi constant
    double m_global;
    Complex m_cR, m_cL;             // X function couplings
    double m_ms[5];                 // masses^2
     
    public:
       
      Tau_Eta_Two_Pion( int nout, ATOOLS::Flavour *fl );
     
      double operator()(const ATOOLS::Vec4D *_p );
      double Using_Hels( const ATOOLS::Vec4D *_p );
      void   SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Doxygen part                                                     %%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*!
    \file Tau_Decays_MEs.H
    \brief Declares the classes HADRONS::Tau_Lepton, HADRONS::Tau_Pion, and HADRONS::Tau_Two_Pion

    This file can be found in directory \c ME_Library
   */

  /*!
    \class Tau_Lepton
    \brief Tools to calculate the ME for leptonic decay channel
   */   
  /*!
    \class Tau_Pion
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi/K\f$
   */
  /*!
    \class Tau_Two_Pion
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi \pi^0\f$
   */   
  /*!
    \class Tau_Three_Pseudo
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \Pi \Pi \Pi\f$
   */   
  /*!
    \class Tau_Four_Pion
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau 4 \pi \f$
   */   
}
