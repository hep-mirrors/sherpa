#ifndef HADRONS_Current_Library_VA_P_V_H
#define HADRONS_Current_Library_VA_P_V_H

#include "HADRONS++/Current_Library/Current_Base.H"

namespace HADRONS
{

namespace VA_P_V_FFs
{
  class FormFactor_Base {
  protected:
    const ATOOLS::Flavour_Vector& m_flavs;
    const std::vector<int>& p_i;
    double m_m0, m_m1;
    double m_A0, m_A1, m_A2, m_A3, m_V;
    bool   m_calced;
  public:
    FormFactor_Base(GeneralModel model, double* masses,
                    const ATOOLS::Flavour_Vector& flavs,
                    const std::vector<int>& indices) :
      m_flavs(flavs), p_i(indices) {
      m_m0=masses[0];
      m_m1=masses[1];
      m_A0=m_A1=m_A2=m_A3=m_V=0.0;
      m_calced=false;
    }
    virtual void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) = 0;
    double A0()  { if(!m_calced) abort(); return m_A0; }
    double A1()  { if(!m_calced) abort(); return m_A1; }
    double A2()  { if(!m_calced) abort(); return m_A2; }
    double A3()  { if(!m_calced) abort(); return m_A3; }
    double V()   { if(!m_calced) abort(); return m_V; }
  };
}

class VA_P_V : public Current_Base {
  double m_Vxx, m_cV;
  VA_P_V_FFs::FormFactor_Base* p_ff;
public:
  VA_P_V(const ATOOLS::Flavour_Vector& flavs,const std::vector<int>& indices,
         const std::string& name) :
    Current_Base(flavs, indices, name) {};
  ~VA_P_V() {
    if (p_ff!=NULL) delete p_ff;
  }
  void SetModelParameters( struct GeneralModel _md );
  void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
};

}
#endif
