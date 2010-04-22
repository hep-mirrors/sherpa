#ifndef HADRONS_Current_Library_VA_0_PiPiPiPi3Charged_H
#define HADRONS_Current_Library_VA_0_PiPiPiPi3Charged_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {
  class VA_0_PiPiPiPi3Charged : public Current_Base {
  // indices: pi+ = 0, pi0 = 1, pi- = 2, pi- = 3
    class LorenzBase {
      public :
//         Complex         * p_X;
       
        const ATOOLS::Vec4D * p_p;
        int             * p_i;
		ATOOLS::Vec4D     m_q;
        double            m_q2;
       
      LorenzBase( ) {};
      virtual ~LorenzBase() {};
      virtual ATOOLS::Vec4C operator()() { return ATOOLS::Vec4C(); };
      void                 SetPrivates(const ATOOLS::Vec4D *, int * i);
    };

	class Novo : public LorenzBase { // Novosibirsk parameterisation hep-ph/0201149, hep-ph/0312240

	  double m_mpi2;
	  SimpleResonanceFlavour m_rho, m_sigma, m_omega;
	  SimpleResonanceFlavour m_a1;
	  ATOOLS::Histogram * p_G, * p_Go, * p_a1width;
	  Complex m_z;
	  double m_Lambda2;
	   
	  double G( double q2 );
	  double Go( double q2 );
	  double F2_a1( double s );
	  Complex Da1( double s );
	  double hrho( double s );
	  Complex Drho( double s );
	  Complex Dsigma( double s );
	  double F2_o( double s );
	  Complex Do( double s );
	  ATOOLS::Vec4C t1( int a, int b, int c, int d );
	  ATOOLS::Vec4C t2( int a, int b, int c, int d );
	  ATOOLS::Vec4C t3( int a, int b, int c, int d );
	   
	  public:

		Novo(std::string path, GeneralModel _md );
          ~Novo() {
            delete p_G;
            delete p_Go;
            delete p_a1width;
          }
		ATOOLS::Vec4C operator()();

	};


    class KS : public LorenzBase {  // CLEO parameterisation hep-ex/9908024 and CERN-TH.6793/93
            // this class does not only provide form factors
            // but the complete Lorentz structure
       

	  ATOOLS::Vec4D     m_r[4];
	  double            m_s[4];
	   
	  static const int  m_ncontrib = 4;   // # contributions to total hadr current
	  Complex      m_Alpha[m_ncontrib];   // weight of each contribution
	  Complex      m_SumAlpha;            // sum of all weights
	   
      double  m_fpi2;                   // pion decay constant
      double  m_Go3p, m_Frho, m_grop;   // coupling constants
      double  m_mpi2, m_mpi02;          // mass^2 of pion, pion0
      ResonanceFlavour m_Rho, m_RR, m_RRR;  // rho resonances
      ResonanceFlavour m_O, m_S, m_F, m_A;  // omega, sigma, f0, a1 resonance
      double  m_beta, m_gamma;          // relative strength of vector resonances
      double  m_sigma;                  // relative strength of vector resonance
      double  m_R[4];                   // important factor
	  ATOOLS::Vec4C OmegaPi();  Complex m_Beta_opi[4];
      ATOOLS::Vec4C AonePi();   Complex m_Beta_api[4];
      ATOOLS::Vec4C SigmaRho(); Complex m_Beta_srh[4];
      ATOOLS::Vec4C FzeroRho(); Complex m_Beta_frh[4];
       
      double  Dots( int, int );
      Complex Fk( double, Complex * );
      Complex Trho( double );
      Complex TTrho( double );
       
      public : 
       
      KS(std::string path, GeneralModel _md );
      ATOOLS::Vec4C operator()();
    }; 

//     int          m_inter[6];            // internal numbering
     
    double       m_global;              // GF/sqrt(2)*Vud
    LorenzBase * p_lorenz;              // pointer to current contributions
     
//     ATOOLS::Vec4D   m_p[6];             // momenta
//     Complex         m_X[5];             // value of an X function
     
  public:
    VA_0_PiPiPiPi3Charged(ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    ~VA_0_PiPiPiPi3Charged() {
      if (p_lorenz!=NULL) delete p_lorenz;
    }
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D * moms);
  };
}
#endif

// to look at: m_r, m_R, m_s, Dots()
