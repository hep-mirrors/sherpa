#ifndef HADRONS_Current_Library_VA_0_PPP_H
#define HADRONS_Current_Library_VA_0_PPP_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {
  class VA_0_PPP : public Current_Base {
     
    class FF_Base {
      double        m_ms[3];            // masses of particles^2
       
      protected :
      // resonances
        ResonanceFlavour m_A, m_AA, m_V[2], m_VV[2];
        double m_alpha, m_Beta[2];            // relative strength of resonances

      // parameters
        double m_fpi2;                        // pion decay constant
        int    m_mode, m_kaon_mode;                // internal mode representation
        bool   m_deltas;              // Delta S 
       
      public :
        FF_Base(int mode, int kaon_mode, std::string path, GeneralModel _md, double * ms );
        virtual ~FF_Base();
        virtual Complex FormFactor( int, double, double, double ) = 0;
        inline double Mass2(int i)    { return m_ms[i]; }
    };
      

    class RChT : public FF_Base {       // Dumm, Pich, Portoles, hep-ph/0312183 
	   
	  private: 
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
        RChT(int mode, int kaon_mode, std::string path, GeneralModel _md, double * _ms);
        ~RChT() {}
        Complex FormFactor( int, double, double, double ); 
    }; 

    class KS : public FF_Base {     // Decker, Finkemeier, Mirkes, hep-ph/9310270
      // order of Pseudos: see table 1 of the paper
		 
	  private: 
	  bool   m_G123;                // false if there is only one amplitude
	  double m_X123, m_ms123; 		// constants for this particular parameterisation

	  // resonances for vector FF (anomaly)
	  ResonanceFlavour m_AnoV, m_AnoVV, m_AnoVVV;
	  double m_BetaV[2];                    // relative strength for V, V', V''
	  double m_AlphaV;                      // relative strength for K* as V_ij (anomalous)
         
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
        KS(int mode, int kaon_mode, std::string path, GeneralModel _md, double * _ms);
        ~KS() {}
        Complex FormFactor( int, double, double, double );
    }; 

	class KS95 : public FF_Base {     // Finkemeier, Mirkes, hep-ph/9503474

	  private: 
	  // resonances in the scalar, vector FF
	  ResonanceFlavour  m_V_VectorFF, m_VV_VectorFF, m_VVV_VectorFF;
	  ResonanceFlavour  m_Omega_VectorFF, m_Phi_VectorFF;
	  double 			m_beta_VectorFF, m_gamma_VectorFF, m_eps_VectorFF;
	   
	  // methods for the FF
	  Complex Axial(double q2);
	  Complex TVector(int i, double s);
	  Complex TOmega(double s);
	  Complex Vector(int i, double s, double u);
	  Complex G3(double q2, double s, double t, double u );

	  public:
		KS95(int mode, int kaon_mode, std::string path, GeneralModel _md, double * _ms);
		~KS95() {}
		Complex FormFactor( int, double, double, double );
	}; 

    int    m_pseudo_3, m_pseudo_1, m_pseudo_2;
    double m_Vud, m_Vus;        // CKM ME
    double m_ms[3];             // masses^2
    int    m_mode, m_kaon_mode; // code for decay mode
    double m_global, m_B123;    // global constant
    FF_Base * p_ff;
    
  public:
    VA_0_PPP(const ATOOLS::Flavour_Vector& flavs,
             const std::vector<int>& indices, const std::string& name);
    ~VA_0_PPP() {
      if (p_ff!=NULL) delete p_ff;
    }
    Complex FormFactor( int j, double Q2, double s, double t ); 
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}
#endif
