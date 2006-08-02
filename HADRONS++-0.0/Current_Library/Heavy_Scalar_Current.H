#ifndef Heavy_Scalar_Current_H
#define Heavy_Scalar_Current_H

#include "Current_Base.H"


namespace HADRONS {

  class Heavy_Scalar_Current : public Current_Base {

    class FF_Base {
      protected:
        double m_mB, m_mS;
        double m_Fplus, m_F0;
        bool   m_calced;
      public:
        FF_Base( struct GeneralModel _md, double* _masses ) {
          m_mB = _masses[0];
          m_mS = _masses[1];
        }
        virtual void CalcFFs( Vec4D pB, Vec4D pS ) = 0;
        double Fplus() { if(!m_calced) abort(); return m_Fplus; }
        double F0()    { if(!m_calced) abort(); return m_F0; }
    };

    class HQET : public FF_Base {
      protected:
        double m_rho2;
        double m_c;
      public:
        HQET( struct GeneralModel _md, double* _masses );
        void CalcFFs( Vec4D pB, Vec4D pS );
    };

    class HQET2 : public FF_Base {
      protected:
        double m_rho2;
        double m_v1_1;
      public:
        HQET2( struct GeneralModel _md, double* _masses );
        void CalcFFs( Vec4D pB, Vec4D pS );
    };

    class ISGW : public FF_Base {
      protected:
        double m_msb;
        double m_msd;
        double m_bb2;
        double m_msq;
        double m_bx2;
        double m_kap2;
      public:
        ISGW( struct GeneralModel _md, double* _masses );
        void CalcFFs( Vec4D pB, Vec4D pS );
    };


    class ISGW2 : public FF_Base {
      protected:
        double m_msb;
        double m_msd;
        double m_bb2;
        double m_mbb;
        double m_nf;
        double m_mbx;
        double m_nfp;
        double m_msq;
        double m_bx2;
      public:
        ISGW2( struct GeneralModel _md, double* _masses );
        void CalcFFs( Vec4D pB, Vec4D pS );
        double Getas( double massq, double massx );
        double GetGammaji ( double z );
    };

    double m_Vxx;
    FF_Base* p_ff;
  public:
    Heavy_Scalar_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}
// indices for B_HeavyVector_Current: 0 must be B, 1 HeavyVector.
#endif
