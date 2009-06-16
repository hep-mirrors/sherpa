#ifndef HADRONS_Current_Library_VA_P_T_H
#define HADRONS_Current_Library_VA_P_T_H

#include "HADRONS++/Current_Library/Current_Base.H"

namespace HADRONS
{

namespace VA_P_T_FFs
{
  class FormFactor_Base {
  protected:
    ATOOLS::Flavour* m_flavs;
    int* p_i;
    double m_m0, m_m1;
    double m_h, m_k, m_bplus, m_bminus;
    bool   m_calced;
  public:
    FormFactor_Base(struct GeneralModel model, double* masses,
                    ATOOLS::Flavour* flavs, int* indices) {
      m_flavs=flavs;
      p_i=indices;
      m_m0=masses[0];
      m_m1=masses[1];
      m_calced=false;
    }
    virtual void CalcFFs( ATOOLS::Vec4D p0, ATOOLS::Vec4D p1 ) = 0;
        double h()      { if(!m_calced) abort(); return m_h; }
        double k()      { if(!m_calced) abort(); return m_k; }
        double bplus()  { if(!m_calced) abort(); return m_bplus; }
        double bminus() { if(!m_calced) abort(); return m_bminus; }
  };
}


class VA_P_T : public Current_Base {
  double m_Vxx;
  VA_P_T_FFs::FormFactor_Base* p_ff;
public:
  VA_P_T(ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
    Current_Base(flavs, n, indices, name) {};
  ~VA_P_T() {
    if (p_ff!=NULL) delete p_ff;
  }
  void SetModelParameters( struct GeneralModel _md );
  void Calc(const ATOOLS::Vec4D * moms);
};

}

#endif
