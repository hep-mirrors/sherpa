#ifndef HADRONS_Current_Library_VA_F_F_H
#define HADRONS_Current_Library_VA_F_F_H

#include "HADRONS++/Current_Library/Current_Base.H"


namespace HADRONS {
  class VA_F_F : public Current_Base {
    class FF_Base {
      public:
        FF_Base( struct GeneralModel _md ) {}
        virtual ~FF_Base();
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
    VA_F_F(const ATOOLS::Flavour_Vector& flavs,const std::vector<int>& indices,
           const std::string& name) :
      Current_Base(flavs, indices, name) {
      m_cL=Complex( 1.0, 0.0);
      m_cR=Complex( 0.0, 0.0);
    };
    ~VA_F_F() {
      if (p_ff!=NULL) delete p_ff;
    }
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}
#endif
// indices for VA_F_F: 0 is "the barred spinor" in the current,
// 1 is the not-barred one
