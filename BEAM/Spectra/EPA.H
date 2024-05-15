#ifndef BEAM_Main_EPA_H
#define BEAM_Main_EPA_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Function_Base.H"

namespace BEAM {

  class EPA final : public Beam_Base {
  private:

    class CosInt final : public ATOOLS::Function_Base {
    public:
      double GetValue(double T) override { return -cos(T)/T; }
      double operator()(double T) override { return GetValue(T); }
      double GetCosInt(double X);
    };

    double m_mass, m_charge, m_gamma, m_minR, m_maxR;
    double m_q2Max, m_q2Min, m_pt_min, m_aqed, m_theta_max;
    double m_xmin, m_xmax;

    bool m_lo_epa;

    void RegisterDefaults() const;

    double phi(double x, double qq) const;

    void selfTest(const std::string& filename);

    int m_formfactor;

  public:

    EPA(const ATOOLS::Flavour beam,const double energy,const double pol,
	const int dir);

    ~EPA() override = default;

    Beam_Base *Copy() override;
    bool CalculateWeight(double, double) override;
    void FixPosition() override;
    void SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i=0) override;
    inline double Xmax() const override { return m_xmax; }
    inline double Xmin() const override { return m_xmin; }
  };// end of class EPA

}// end of namespace BEAM

#endif
