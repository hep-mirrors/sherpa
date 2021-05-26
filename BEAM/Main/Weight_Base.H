#ifndef BEAM_Main_Weight_Base_H
#define BEAM_Main_Weight_Base_H

#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Kinematics_Base.H"
#include "ATOOLS/Math/Vector.H"

namespace BEAM {
  class Weight_Base {
  protected:
    Kinematics_Base * p_kinematics;
    Beam_Base      ** p_beams;
    double            m_weight;
    std::string       m_keyid;
  public:
    Weight_Base(Kinematics_Base * kinematics);
    ~Weight_Base();

    virtual void AssignKeys(ATOOLS::Integration_Info *const info)=0;
    virtual bool Calculate(const double & scale)=0;
    virtual const double operator()(ATOOLS::Flavour * flin=NULL) const { return m_weight; }
  };
}
#endif