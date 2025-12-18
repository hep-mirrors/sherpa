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
    EPA_ff_type  m_fftype;
    EPA_FF_Base* p_ff;
    double       m_pref;
    double       m_pt2max;

    bool m_plotting;

    void Initialise();
    void RegisterDefaults() const;

    static void Tests();

  public:
    EPA(const ATOOLS::Flavour& beam, double energy, double pol, int dir);
    ~EPA() override = default;
    Beam_Base* Copy() override { return new EPA(*this); }

    bool CalculateWeight(double x, double q2) override;
    void FixPosition() override;
    void SetOutMomentum(const ATOOLS::Vec4D& out) override;

    double Xmax() const override { return p_ff->Xmax(); }
    double Xmin() const override { return p_ff->Xmin(); }

    double ImpactParameter() const override { return p_ff->ImpactParameter(); }
  };// end of class EPA

}// end of namespace BEAM

#endif
