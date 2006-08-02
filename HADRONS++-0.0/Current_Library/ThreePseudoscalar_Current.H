#ifndef ThreePseudoscalar_Current_H
#define ThreePseudoscalar_Current_H

#include "Current_Base.H"
#include "ResonanceFlavour.H"
#include "Tools.H"

namespace HADRONS {
  class ThreePseudoscalar_Current : public Current_Base {
     
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

    int    m_pseudo_3, m_pseudo_1, m_pseudo_2;
    double m_Vud, m_Vus;        // CKM ME
    double m_ms[3];             // masses^2
    int    m_mode;              // code for decay mode
    double m_global, m_B123;    // global constant
    FF_Base * p_ff;
    
  public:
    ThreePseudoscalar_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name);
    ~ThreePseudoscalar_Current() {
      if (p_ff!=NULL) delete p_ff;
    }
    Complex FormFactor( int j, double Q2, double s, double t ); 
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}
#endif
