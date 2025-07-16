#ifndef HADRONS_Current_Library_VA_P_V_H
#define HADRONS_Current_Library_VA_P_V_H

#include "HADRONS++/Current_Library/Current_Base.H"

namespace HADRONS
{

namespace VA_P_A_FFs
{
  class FormFactor_Base {
  protected:
    const ATOOLS::Flavour_Vector& m_flavs;
    const std::vector<int>& p_i;
    double m_m0, m_m1;
    double m_V0, m_V1, m_V2, m_A;
    bool   m_calced;
  public:
    FormFactor_Base(GeneralModel model, double* masses,
                    const ATOOLS::Flavour_Vector& flavs,
                    const std::vector<int>& indices) :
      m_flavs(flavs), p_i(indices) {
      m_m0=masses[0];
      m_m1=masses[1];
      m_V0=m_V1=m_V2=m_A=0.0;
      m_calced=false;
    }
    virtual ~FormFactor_Base();
    virtual void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) = 0;
    double V0()  { if(!m_calced) abort(); return m_V0; }
    double V1()  { if(!m_calced) abort(); return m_V1; }
    double V2()  { if(!m_calced) abort(); return m_V2; }
    double A()   { if(!m_calced) abort(); return m_A; }
  };
}

class VA_P_A : public Current_Base {
  double m_Vxx, m_cV;
  VA_P_A_FFs::FormFactor_Base* p_ff;
public:
  VA_P_A(const ATOOLS::Flavour_Vector& flavs,const std::vector<int>& indices,
         const std::string& name) :
    Current_Base(flavs, indices, name) {};
  ~VA_P_A() {
    if (p_ff!=NULL) delete p_ff;
  }
  void SetModelParameters( struct GeneralModel _md );
  void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
};

}
#endif
