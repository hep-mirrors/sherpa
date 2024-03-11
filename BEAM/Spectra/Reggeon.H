#ifndef BEAM_Main_Reggeon_H
#define BEAM_Main_Reggeon_H

#include "BEAM/Main/Beam_Base.H"

namespace BEAM {

  class Reggeon final : public Beam_Base {
  private:
    double m_tMax, m_tMin, m_xMax, m_xMin, m_proton_mass;
    double m_A, m_B, m_alpha_intercept, m_alpha_slope, m_n;

  public:
    Reggeon(const ATOOLS::Flavour beam, const double energy, const double pol,
            const int dir);
    ~Reggeon() override = default;

    Beam_Base* Copy() override;
    bool CalculateWeight(double, double) override;
    void SetOutMomentum(const ATOOLS::Vec4D& out, const size_t& i = 0) override;
    void FixPosition() override;
    void FixNormalisation();
    inline double Xmax() const override { return m_xMax; }
    inline double Xmin() const override { return m_xMin; }
  };// end of class Reggeon
}// end of namespace BEAM
#endif
