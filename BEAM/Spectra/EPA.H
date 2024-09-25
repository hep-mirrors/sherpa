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
    double       m_mass, m_charge, m_gamma;
    double       m_aqed, m_pref;
    double       m_x, m_q2, m_pt2max;
    double       m_Rmin, m_Rmax, m_Rcut, m_xmin, m_xmax, m_bmin, m_bmax;

    double m_Q02;

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
    void SetOutMomentum(const ATOOLS::Vec4D& out, const size_t& i = 0) override;

    inline double Xmax() const override { return m_xmax; }
    inline double Xmin() const override { return m_xmin; }
  };// end of class EPA

}// end of namespace BEAM

#endif
