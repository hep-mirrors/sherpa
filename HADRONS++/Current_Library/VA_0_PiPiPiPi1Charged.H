#ifndef HADRONS_Current_Library_VA_0_PiPiPiPi1Charged_H
#define HADRONS_Current_Library_VA_0_PiPiPiPi1Charged_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Math/Histogram.H"

namespace HADRONS {
  class VA_0_PiPiPiPi1Charged : public Current_Base {
  // indices: pi- = 0, pi0 = 1, pi0 = 2, pi0 = 3

    class LorenzBase {
      public :
        const ATOOLS::Vec4D   * p_p;
        int * p_i;
	ATOOLS::Vec4D     m_q;
        ATOOLS::Vec4D     m_r[4];
        double            m_s[4];
        double            m_q2;
       
      LorenzBase( ) {};
      virtual ~LorenzBase() {};
      virtual ATOOLS::Vec4C operator()() { return ATOOLS::Vec4C(); };
      void SetPrivates(const ATOOLS::Vec4D *, int* i);
    };

    class Novo : public LorenzBase {  // Novosibirsk parameterisation hep-ex/0201149, hep-ph/0312240

	  double m_mpi2;
	  SimpleResonanceFlavour m_rho, m_sigma;
	  SimpleResonanceFlavour m_a1;
	  ATOOLS::Histogram * p_G, * p_a1width, * p_zforma1;
	  Complex m_z;
	  double m_Lambda2;
	   
	  double G( double q2 );
	  double F2_a1( double s );
	  Complex Da1( double s );
	  double hrho( double s );
	  Complex Drho( double s );
	  Complex Dsigma( double s );
	  ATOOLS::Vec4C t1( int a, int b, int c, int d );
	  ATOOLS::Vec4C t2( int a, int b, int c, int d );

      public : 
      Novo(GeneralModel _md, std::string path);
      ~Novo() {
        delete p_G;
        delete p_zforma1;
        delete p_a1width;
      }
      ATOOLS::Vec4C operator()();
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
      ATOOLS::Vec4C operator()();
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    int          m_nutau, m_pion1, m_pion2, m_pion3, m_pion0;
     
    double       m_global;
    LorenzBase * p_lorenz;              // pointer to current contributions
     
  public:
    VA_0_PiPiPiPi1Charged(ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    ~VA_0_PiPiPiPi1Charged() {
      if (p_lorenz!=NULL) delete p_lorenz;
    }
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D * moms);
  };
}
#endif

// to look at: m_r, m_R, m_s, Dots()
