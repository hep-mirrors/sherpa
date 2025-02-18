#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Spectra/EPA_FF.H"

namespace ATOOLS {
  class Settings;
}

namespace BEAM {
  class EPA final : public Beam_Base {
  private:
    EPA_ff_type  m_type;
    EPA_FF_Base* p_ff;
    double       m_mass;
    double       m_aqed, m_pref;
    double       m_q2, m_pt2max;
    double       m_xmin, m_xmax;
    double       m_b;

    size_t m_plotting;

    void Initialise();
    void RegisterDefaults() const;

    void Tests();

  public:
    EPA(const ATOOLS::Flavour& beam, double energy, double pol, int dir);
    ~EPA() override = default;
    Beam_Base* Copy() override { return new EPA(*this); }

    bool CalculateWeight(double x, double q2) override;
    void FixPosition() override;
    void SetOutMomentum(const ATOOLS::Vec4D& out) override;

    inline double Xmax() const override { return m_xmax; }
    inline double Xmin() const override { return m_xmin; }
  };// end of class EPA

}// end of namespace BEAM

#endif
