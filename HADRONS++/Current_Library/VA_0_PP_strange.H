#ifndef HADRONS_Current_Library_VA_0_PP_strange_H
#define HADRONS_Current_Library_VA_0_PP_strange_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"


namespace HADRONS {
  class VA_0_PP_strange : public Current_Base {
    class FF_Base {
      protected:
        ResonanceFlavour m_R, m_R0;
        double m_fpi2;
        double m_mPi2, m_mK2, m_mEta2, m_mPi, m_mK, m_mEta;  // mass of pi, K, eta
        double m_Sigma_KP, m_Delta_KP;    // sum, difference of masses^2
      public :
        FF_Base(GeneralModel _md);
        virtual ~FF_Base();
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
	  double MuOf( double m2 );
      double MassWidthVector( double s );
      double MassWidthScalar( double s );
       
      public:
        RChT(GeneralModel _md);
        ~RChT() {};
        Complex VectorFormFactor( double );
        Complex ScalarFormFactor( double );
    }; 

    class KS : public FF_Base {     
       
	  ResonanceFlavour  m_RR;
	  double 			m_beta;
      public:
        KS(GeneralModel _md);
        ~KS() {};
        Complex VectorFormFactor( double );
        Complex ScalarFormFactor( double );
    }; 

    bool     m_chpionmode;      // pion = pi+ or pi0 
//     double   m_Vus2;            // CKM ME
    double   m_global;          //
    double   m_ms[4];           // masses^2 
    double   m_Delta_KP;        // difference of masses^2

    FF_Base * p_ff;
    
  public:
    VA_0_PP_strange(const ATOOLS::Flavour_Vector& flavs,
                    const std::vector<int>& indices, const std::string& name) :
      Current_Base(flavs, indices, name) {};
    ~VA_0_PP_strange() {
      if (p_ff!=NULL) delete p_ff;
    }
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}
#endif
