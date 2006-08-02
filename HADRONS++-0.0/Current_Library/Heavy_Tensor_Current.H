#ifndef Heavy_Tensor_Current_H
#define Heavy_Tensor_Current_H

#include "Current_Base.H"

namespace HADRONS {

  class Heavy_Tensor_Current : public Current_Base {

    class FF_Base {
      protected:
        double m_mB, m_mT;
        double m_h, m_k, m_bplus, m_bminus;
        bool   m_calced;
      public:
        FF_Base( struct GeneralModel _md, double* _masses ) {
          m_mB   = _masses[0];
          m_mT   = _masses[1];
          m_calced = false;
        }
        virtual void CalcFFs( double q2 ) = 0;
        double h() { if(!m_calced) abort(); return m_h; }
        double k() { if(!m_calced) abort(); return m_k; }
        double bplus() { if(!m_calced) abort(); return m_bplus; }
        double bminus() { if(!m_calced) abort(); return m_bminus; }
    };

    class ISGW : public FF_Base {
      protected:
         double m_msb;
         double m_msd;
         double m_bb2;
         double m_msq;
         double m_bx2;
      public:
        ISGW( struct GeneralModel _md, double* _masses );
        void CalcFFs( double q2 );
    };

    class ISGW2 : public FF_Base {
      protected:
         double m_msb;
         double m_msd;
         double m_bb2;
         double m_mbb;
         double m_nf;
         double m_msq;
         double m_bx2;
         double m_mbx;
         double m_nfp;
      public:
        ISGW2( struct GeneralModel _md, double* _masses );
        void CalcFFs( double q2 );
      double Getas( double mass );   
    };

    double m_Vxx;
    FF_Base* p_ff;
  public:
    Heavy_Tensor_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}

#endif
