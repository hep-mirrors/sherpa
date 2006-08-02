#ifndef V_A_Current_H
#define V_A_Current_H

#include "Current_Base.H"


namespace HADRONS {
  class V_A_Current : public Current_Base {
    class FF_Base {
      public:
        FF_Base( struct GeneralModel _md ) {}
        virtual double ff( double q2 ) = 0;
    };

    class Simple_FF : public FF_Base { /* dummy */
      public:
      Simple_FF( struct GeneralModel _md );
      double ff( double q2 ) { return 1.0; }
    };

    Complex m_cR, m_cL;
    FF_Base* p_ff;
  public:
    V_A_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {
      m_cL=Complex( 1.0, 0.0);
      m_cR=Complex( 0.0, 0.0);
    };
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}
#endif
// indices for V_A_Current: 0 is "the barred spinor" in the current,
// 1 is the not-barred one
