#ifndef Heavy_Vector_Current_H
#define Heavy_Vector_Current_H

#include "Current_Base.H"

namespace HADRONS {

  class Heavy_Vector_Current : public Current_Base {

    class FF_Base {
      protected:
        double m_mB, m_mV;
        double m_A0, m_A1, m_A2, m_A3, m_V;
        bool   m_calced;
      public:
        FF_Base( struct GeneralModel _md, double* _masses ) {
          m_mB   = _masses[0];
          m_mV   = _masses[1];
          m_calced = false;
        }
        virtual void CalcFFs( double q2 ) = 0;
        double A0() { if(!m_calced) abort(); return m_A0; }
        double A1() { if(!m_calced) abort(); return m_A1; }
        double A2() { if(!m_calced) abort(); return m_A2; }
        double A3() { if(!m_calced) abort(); return m_A3; }
        double V()  { if(!m_calced) abort(); return m_V; }
    };

    class HQET : public FF_Base {
      protected:
        double m_R1;
        double m_R2;
        double m_RStar;
        double m_rho2;
      public:
        HQET( struct GeneralModel _md, double* _masses );
        void CalcFFs( double q2 );
    };

    class HQET2 : public FF_Base {
      protected:
        double m_r1_1;
        double m_r2_1;
        double m_ha1_1;
        double m_rho2;
      public:
        HQET2( struct GeneralModel _md, double* _masses );
        void CalcFFs( double q2 );
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
        void CalcFFs( double q2 );
    };

    class ISGW2 : public FF_Base {
      protected:
        double m_msb;
        double m_msd;
        double m_bb2;
        double m_msq;
        double m_bx2;
        double m_mbb;
        double m_nf;
        double m_cf;
        double m_mbx;
        double m_nfp;
      public:
        ISGW2( struct GeneralModel _md, double* _masses );
        void CalcFFs( double q2 );
        double Getas( double massq, double massx );
        double GetGammaji ( double z );
    };

    double m_Vxx;
    FF_Base* p_ff;
  public:
    Heavy_Vector_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}
// indices for Heavy_Vector_Current: 0 must be B, 1 HeavyVector.
#endif
