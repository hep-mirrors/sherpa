#ifndef BEAM_Main_Target_Matter_H
#define BEAM_Main_Target_Matter_H

#include "BEAM/Main/Beam_Base.H"

namespace BEAM {
  class Target_Matter : public Beam_Base {
  private:
    double m_Z, m_A;

    static ATOOLS::Flavour SampleFlavour(const double& Z, const double& A);

  public:
    Target_Matter(const double& Z, const double& A, const double& energy,
                  const double& polarisation, const int& dir);

    Beam_Base * Copy();
    ATOOLS::Flavour Remnant() { return ATOOLS::Flavour(kf_none); }
    bool CalculateWeight(const double, const double) { return true; }

    bool IsMultiSpecies() const { return true; }
    ATOOLS::Flavour DrawFlavour() const { return SampleFlavour(m_Z, m_A); }
  };
}
#endif
