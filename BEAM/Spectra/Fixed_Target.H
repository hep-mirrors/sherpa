#ifndef BEAM_Main_Fixed_Target_H
#define BEAM_Main_Fixed_Target_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Function_Base.H"

namespace BEAM {

  
  class Fixed_Target : public Beam_Base {
  private:
    void RegisterDefaults();
    void selfTest(std::string filename);

  public:
    Fixed_Target(const ATOOLS::Flavour beam, const double &ene,
    const double &pol,  const int &dir);
    ~Fixed_Target();

    inline Beam_Base *Copy() { return new Fixed_Target(*this); }

    bool CalculateWeight(double, double);

    ATOOLS::Vec4D m_beam_lab, m_fixed_target;
  };
}

#endif
