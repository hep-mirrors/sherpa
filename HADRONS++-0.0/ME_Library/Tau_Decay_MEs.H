#ifndef Tau_Decays_h
#define Tau_Decays_h
 
#include "HD_ME_Base.H"
#include "Flavour.H"
#include "Vector.H"
#include "MathTools.H"
#include "Histogram.H"
#include "Message.H"
#include "Tools.H"
#include <utility>
#include "ResonanceFlavour.H"

namespace HADRONS {

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Lepton : public HD_ME_Base {
    int m_nutau, m_nulep, m_lep;
    double m_global;                        // GF/sqrt(2)
    Complex m_cR1, m_cL1, m_cR2, m_cL2;     // couplings for Z function
    public:
      Tau_Lepton( int nout, ATOOLS::Flavour *fl );
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void   SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pseudo : public HD_ME_Base {
    int     m_nutau, m_pion;
    bool    m_pionmode;             // either tau->pi or tau->kaon
    double  m_global;               // GF * f_had * V_CKM
    Complex m_cR, m_cL;             // coupling for X function
    public:
      Tau_Pseudo( int nout, ATOOLS::Flavour *fl );
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      inline void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Two_Pion : public HD_ME_Base {
    int    m_nutau, m_pion_ch, m_pion0;
    bool   m_pionmode;              // either tau->2pi or tau->2kaon
    int    m_ff;                    // ff model
    double m_global;                // GF * V_CKM * Clebsch Gordon
    double m_fpi;                   // pion decay constant
    double m_frho, m_grpp;          // rho-pipi coupling
    Complex m_cR, m_cL;             // coupling for X function
    ResonanceFlavour m_R, m_RR, m_RRR;  // resonances
    double m_beta, m_gamma;         // rel. strength
    double m_gammaR, m_gammaRR, m_gammaRRR;
    double m_m, m_m2;               // mass of outgoing particle
    double m_m2_pi, m_m2_K;         // mass of pion, kaon
    Complex A( double x, double y );
    Complex FormFactor( double s );
    public:
      Tau_Two_Pion( int nout, ATOOLS::Flavour *fl );
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Pion_Kaon : public HD_ME_Base {
     
    class FF_Base {
      protected:
        ResonanceFlavour m_R, m_R0;
        double m_fpi2;
        double m_mPi2, m_mK2, m_mEta2, m_mPi, m_mK, m_mEta;  // mass of pi, K, eta
        double m_Sigma_KP, m_Delta_KP;    // sum, difference of masses^2
      public :
        FF_Base(GeneralModel _md);
        virtual Complex VectorFormFactor( double ) = 0;
        virtual Complex ScalarFormFactor( double ) = 0;
        void SetMasses2( double, double, double );
    };
      

    class RChT : public FF_Base {       // Guerrero, Tesis Doctoral
       
      double m_MK2, m_GK;               // K* resonance par's, pi decay constant
      double m_MK02, m_GK0;             // K0* resonance par's
      double m_renorm2;                 // renormalisation scale
      double m_cd, m_cm;                // LEC
      
      Complex JBar( double s, double MP2, double MQ2, double Sigma, double Delta );
      Complex JBarBar( double s, double MP2, double MQ2, double Sigma, double Delta );
      Complex Mr( double s, double MP2, double MQ2 );
      Complex L( double s, double MP2, double MQ2 );
      double MassWidthVector( double s );
      double MassWidthScalar( double s );
       
      public:
      RChT(GeneralModel _md);
      ~RChT() {};
      Complex VectorFormFactor( double );
      Complex ScalarFormFactor( double );
    }; 

    class KS : public FF_Base {     
       
      public:
      KS(GeneralModel _md);
      ~KS() {};
      Complex VectorFormFactor( double );
      Complex ScalarFormFactor( double );
    }; 

    int      m_nutau, m_pion, m_kaon;
    bool     m_chpionmode;      // pion = pi+ or pi0 
    double   m_Vus2;            // CKM ME
    double   m_global;          // GF * Vus / 2.;
    Complex  m_cR, m_cL;        // coupling for X function
    double   m_ms[4];           // masses^2 
    double   m_Delta_KP;        // difference of masses^2

    FF_Base * p_ff;
     
    public:
      Tau_Pion_Kaon( int nout, ATOOLS::Flavour *fl );
      ~Tau_Pion_Kaon() { if (p_ff!=NULL) delete p_ff; }
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Three_Pseudo : public HD_ME_Base {
     
    class FF_Base {
      double        m_ms[3];            // masses of particles^2
       
      protected :
      // resonances
      ResonanceFlavour m_A, m_AA, m_V[2], m_VV[2], m_AnoV, m_AnoVV, m_AnoVVV;
      double m_alpha, m_Beta[2];            // relative strength of resonances
      double m_BetaV[2];                    // relative strength for V, V', V''
      double m_AlphaV;                      // relative strength for K* as V_ij (anomalous)

      // parameters
      double m_fpi2;                        // pion decay constant
      int    m_mode;                // internal mode representation
      bool   m_deltas;              // Delta S 
      bool   m_G123;                // false if there is only one amplitude
      double m_X123, m_ms123;
       
      public :
      FF_Base(int mode, std::string path, GeneralModel _md, double * ms );
      ~FF_Base() {}
      virtual Complex FormFactor( int, double, double, double ) = 0;
      inline double Mass2(int i)    { return m_ms[i]; }
    };
      

    class RChT : public FF_Base {       // Dumm, Pich, Portoles, hep-ph/0312183 
      double  MassWidthVector( double s ); 
      double  MassWidthVector( int a, double s ); 
      double  MassWidthAxial( double Q2 ); 
       
      double  FFunc( double, double, double );
       
      // resonances
      double m_MO, m_MO2, m_GO;
      double m_gammaR;                  // global constant for GV
       
      // parameters
      double m_l0, m_exp_alpha;         // fitting parameter
      double m_m, m_m2, m_mK2;          // masses of pion and kaon
      double m_GV, m_FV, m_FV2;         // couplings
      double m_FA, m_FA2;
      double m_lsum;                    // fit parameters
      double m_l1, m_l2;
       
      public:
      RChT(int mode, std::string path, GeneralModel _md, double * _ms);
      ~RChT() {}
      Complex FormFactor( int, double, double, double ); 
    }; 

    class KS : public FF_Base {     // Decker, Finkemeier, Mirkes, hep-ph/9310270
      // order of Pseudos: see table 1 of the paper
         
      // methods for axial and scalar FF
      Complex BW_A( double s );                     // resonance A
      Complex BW_V( int a, double s );              // resonance V_ij
      Complex BW_VV( int a, double s );             // resonance V_ij'
      Complex Tvector1( int a, int b, double x );
      Complex Tvector2( int a, int b, double x );
      Complex TSvector( int a, int b, int c, double Q2, double s, double t );
      Complex Tgen( int a, int b, int c, double s, double t);
      // methods for vector FF
      Complex BW_Vano( int a, double q2 );
      Complex T_V( double q2 );
      Complex T_V13V23( double s, double t );

      public:
      KS(int mode, std::string path, GeneralModel _md, double * _ms);
      ~KS() {}
      Complex FormFactor( int, double, double, double );
    }; 

    int    m_nutau, m_pseudo_3, m_pseudo_1, m_pseudo_2;
    double m_Vud, m_Vus;        // CKM ME
    Complex m_cR, m_cL;         // X function couplings
    double m_ms[5];             // masses^2
    int    m_mode;              // code for decay mode
    double m_global, m_B123;    // global constant
    FF_Base * p_ff;
    public:
       
      Tau_Three_Pseudo( int nout, ATOOLS::Flavour *fl );
      ~Tau_Three_Pseudo() {
        if (p_ff!=NULL) delete p_ff;
      }
     
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
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
      void              SetPrivates( Complex *, ATOOLS::Vec4D * );
    };

    class KS : public LorenzBase {  // CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
            // this class does not only provide form factors
            // but the complete Lorentz structure
       

      double  m_fpi2;                   // pion decay constant
      double  m_Go3p, m_Frho, m_grop;   // coupling constants
      double  m_mpi2, m_mpi02;          // mass^2 of pion, pion0
      ResonanceFlavour m_Rho, m_RR, m_RRR;  // rho resonances
      ResonanceFlavour m_O, m_S, m_F, m_A;  // omega, sigma, f0, a1 resonance
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
       
      KS(std::string path, GeneralModel _md );
      Complex operator()( int );
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
    double       m_ms[6];
    int          m_inter[6];            // internal numbering
     
    double       m_global;              // GF/sqrt(2)*Vud
    Complex      m_cR, m_cL;            // couplings in X functions
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
   
    void   operator()( 
        const ATOOLS::Vec4D  * _p, 
        std::vector<Complex> * _ampls_tensor, 
        std::vector<std::pair<int,int> > * _indices,
        int k0_n );
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
      void SetPrivates( Complex *, ATOOLS::Vec4D * );
    };

    class KS : public LorenzBase {  // CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
       

      double  m_fpi2;                   // pion decay constant
      double  m_mpi2, m_mpi02;          // mass^2 of pion, pion0
      ResonanceFlavour m_Rho, m_RR, m_RRR;  // rho resonances
      double  m_beta, m_gamma;          // relative strength of vector resonances
       
      double  Dots( int, int );
      Complex Trho( double );
       
      public : 
       
      KS(GeneralModel _md);
      Complex operator()();
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
    double       m_ms[6];
    int          m_inter[6];            // internal numbering
     
    double       m_global;              // GF/sqrt(2)*Vud
    Complex      m_cR, m_cL;                // couplings in X functions
    LorenzBase * p_lorenz;              // pointer to current contributions
     
    ATOOLS::Vec4D   m_p[6];             // momenta
    Complex         m_X[5];             // value of an X function
     
    public:
     
    Tau_Four_Pion_1( int nout, ATOOLS::Flavour *fl );
    ~Tau_Four_Pion_1() {
      if (p_lorenz!=NULL) delete p_lorenz;
    }
   
    void   operator()( 
        const ATOOLS::Vec4D  * _p, 
        std::vector<Complex> * _ampls_tensor, 
        std::vector<std::pair<int,int> > * _indices,
        int k0_n );
    void SetModelParameters( GeneralModel _md );

  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  class Tau_Eta_Two_Pion : public HD_ME_Base {
     
    int    m_nutau, m_eta, m_pion, m_pion0;
    double m_global;                // 2/3*GF*Vud/fpi
    Complex m_cR, m_cL;             // X function couplings
    double m_ms[5];                 // masses^2
     
    public:
       
      Tau_Eta_Two_Pion( int nout, ATOOLS::Flavour *fl );
     
      void   operator()( 
          const ATOOLS::Vec4D  * _p, 
          std::vector<Complex> * _ampls_tensor, 
          std::vector<std::pair<int,int> > * _indices,
          int k0_n );
      void   SetModelParameters( GeneralModel _md );
  };

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// Doxygen part                                                     %%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  /*!
    \file Tau_Decay_MEs.H
    \brief Declares the ME classes for tau decays

    This file can be found in the directory \c ME_Library

    The implemted classes are
    - HADRONS::Tau_Lepton
    - HADRONS::Tau_Pseudo
    - HADRONS::Tau_Two_Pion
    - HADRONS::Tau_Pion_Kaon
    - HADRONS::Tau_Three_Pseudo
    - HADRONS::Tau_Four_Pion_1
    - HADRONS::Tau_Four_Pion_3
    - HADRONS::Tau_Eta_Two_Pion
    .
   */

  /*!
    \class Tau_Eta_Two_Pion
    \brief Tools to calculate the ME for eta-two-pion decay channel
   */   
  /*!
    \class Tau_Lepton
    \brief Tools to calculate the ME for leptonic decay channel
   */   
  /*!
    \class Tau_Pseudo
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi/K\f$
   */
  /*!
    \class Tau_Two_Pion
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \pi \pi^0\f$, \f$\tau \to \nu_\tau K K^0\f$
   */   
  /*!
    \class Tau_Pion_Kaon
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau K\pi\f$
   */   
  /*!
    \class Tau_Three_Pseudo
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau \Pi \Pi \Pi\f$
   */   
  /*!
    \class Tau_Three_Pseudo::FF_Base
    \brief Base class for form factor parameterisations
   */   
  /*!
    \class Tau_Three_Pseudo::RChT
    \brief RChT parameteristion 
   */   
  /*!
    \class Tau_Three_Pseudo::KS
    \brief KS parameteristion 
   */   
  /*!
    \class Tau_Four_Pion_1
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau 4 \pi \f$ (one-prong)
   */   
  /*!
    \class Tau_Four_Pion_3
    \brief Tools to calculate the ME for \f$\tau \to \nu_\tau 4 \pi \f$ (three-prong)
   */   
}

#endif
