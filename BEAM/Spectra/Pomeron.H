#ifndef BEAM_Main_Pomeron_H
#define BEAM_Main_Pomeron_H

#include "BEAM/Main/Beam_Base.H"

namespace BEAM {

class Pomeron final : public Beam_Base {
private:
  double m_tMax, m_xMax, m_proton_mass;
  double m_A, m_B, m_alpha_intercept, m_alpha_slope;
public:
  Pomeron(const ATOOLS::Flavour beam,const double energy,const double pol,
      const int dir);
  ~Pomeron() override = default;

  Beam_Base *Copy() override;
  bool CalculateWeight(double, double) override;
  void SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i=0) override;
  void FixPosition() override;
};// end of class Pomeron
}// end of namespace BEAM
#endif
