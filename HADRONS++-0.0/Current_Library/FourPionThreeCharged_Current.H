#ifndef FourPionThreeCharged_Current_H
#define FourPionThreeCharged_Current_H

#include "Current_Base.H"
#include "ResonanceFlavour.H"
#include "Tools.H"

namespace HADRONS {
  class FourPionThreeCharged_Current : public Current_Base {

    class LorenzBase {
      public :
//         Complex         * p_X;
       
        ATOOLS::Vec4D   * p_p;
        ATOOLS::Vec4D     m_q;
        ATOOLS::Vec4D     m_r[4];
        double            m_s[4];
        double            m_q2;
       
        LorenzBase( ) {};
        virtual ComplexVec4D   operator()( int ) {};
        void              SetPrivates( ATOOLS::Vec4D * );
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
      ComplexVec4D OmegaPi();  Complex m_Beta_opi[4];
      ComplexVec4D AonePi();   Complex m_Beta_api[4];
      ComplexVec4D SigmaRho(); Complex m_Beta_srh[4];
      ComplexVec4D FzeroRho(); Complex m_Beta_frh[4];
       
      double  Dots( int, int );
      Complex Fk( double, Complex * );
      Complex Trho( double );
      Complex TTrho( double );
       
      public : 
       
        KS(std::string path, GeneralModel _md );
        ComplexVec4D operator()( int );
    }; 

    static const int  m_ncontrib = 4;   // # contributions to total hadr current
    double       m_ms[4];
//     int          m_inter[6];            // internal numbering
    
    double       m_global;
    LorenzBase * p_lorenz;              // pointer to current contributions
    Complex      m_Alpha[m_ncontrib];   // weight of each contribution
    Complex      m_SumAlpha;            // sum of all weights
     
//     ATOOLS::Vec4D   m_p[4];             // momenta
//     Complex         m_X[5];             // value of an X function
     
    public:
     
      FourPionThreeCharged_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name);
      ~FourPionThreeCharged_Current() {
        if (p_lorenz!=NULL) delete p_lorenz;
      }
      void SetModelParameters( struct GeneralModel _md );
      void Calc();
  };

}
#endif
