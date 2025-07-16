#ifndef BEAM_Main_DM_beam_H
#define BEAM_Main_DM_beam_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Function_Base.H"

namespace BEAM {
  struct DM_type {
    enum code {
      Boltzmann = 1,
      none      = 0
    };
  };
  std::ostream &operator<<(std::ostream &s,const DM_type::code & type);
  
  class DM_beam : public Beam_Base {
  private:
    double m_mass, m_temperature, m_BesselK2;
    DM_type::code  m_formfactor;
    bool   m_relativistic, m_debugging; // Whether or not to multiply by v
    // long double m_bessel_2; // result of K(2,m/T)

    void RegisterDefaults();
    void selfTest(std::string filename);
  public:
    DM_beam(const ATOOLS::Flavour beam,const double & temp,
	    const DM_type::code & form, const bool & relativistic,
	    const int dir);
    ~DM_beam() {}

    inline Beam_Base *Copy() { return new DM_beam(*this); }

    bool CalculateWeight(double, double);
    void SetTemperature(const double & temp) { m_temperature = temp; }
    
    inline const double & Temperature() { return m_temperature; }
  };
}

#endif
