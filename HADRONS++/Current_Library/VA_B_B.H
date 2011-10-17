#ifndef HADRONS_Current_Library_VA_B_B_H
#define HADRONS_Current_Library_VA_B_B_H

#include "HADRONS++/Current_Library/Current_Base.H"


namespace HADRONS {

namespace VA_B_B_FFs
{
  class FormFactor_Base {
    protected:
      double m0, m1;
      double m_V1, m_V2, m_V3, m_A1, m_A2, m_A3;
      bool   m_calced;
    public:
      FormFactor_Base(struct GeneralModel model, double* masses) :
        m0(masses[0]), m1(masses[1]), m_calced(false) {}
      virtual void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) = 0;
      double V1() { if(!m_calced) abort(); return m_V1; }
      double V2() { if(!m_calced) abort(); return m_V2; }
      double V3() { if(!m_calced) abort(); return m_V3; }
      double A1() { if(!m_calced) abort(); return m_A1; }
      double A2() { if(!m_calced) abort(); return m_A2; }
      double A3() { if(!m_calced) abort(); return m_A3; }
  };
}

class VA_B_B : public Current_Base
{
  double m_V_CKM;
  VA_B_B_FFs::FormFactor_Base* p_ff;
  bool m_unnatural; // unnatural parity (1/2^+ --> 1/2^-)
  double m_v, m_a;

public:
  VA_B_B(const ATOOLS::Flavour_Vector& flavs,const std::vector<int>& indices,
         const std::string& name) :
    Current_Base(flavs, indices, name) {};
  ~VA_B_B();
  void SetModelParameters( struct GeneralModel _md );
  void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
};
}

#endif
