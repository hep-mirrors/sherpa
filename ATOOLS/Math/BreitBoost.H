#ifndef ATOOLS_Phys_BreitBoost_H
#define ATOOLS_Phys_BreitBoost_H
#include "ATOOLS/Math/Poincare.H"

namespace ATOOLS {

  class Cluster_Amplitude;

  class BreitBoost : public Poincare_Sequence {
  private:
    /// helper doing the actual work in the constructors
    void _init(const Vec4D& lepin, const Vec4D& lepout, const Vec4D& hadin);

    /// save dis variables
    double m_Q2, m_x;

  public:
    /*
      Construct a boost into the Breit frame.

      The breit frame is defined by the incoming \p lepin and outgoing \p lepout
      lepton momenta, and the incoming hadron momentum \p hadin.

      Incoming momenta are assumed incoming, i.e. with positive energies. The
      hadron momentum should be projected to zero mass if necessary beforehand.
    */
    BreitBoost(const Vec4D& lepin, const Vec4D& lepout, const Vec4D& hadin);

    /*
      Constructor using information from the amplitude.

      Convenience funtion, takes lepin and lepout from \p ampl, and
      hadin from the global beam information, projected to zero mass.
    */
    BreitBoost(Cluster_Amplitude *const ampl);

    /// Convenience function to apply boost to Cluster_Amplitude
    void Apply(Cluster_Amplitude *const ampl) const;

    /// access dis variables
    inline double Q2() const { return m_Q2; }
    inline double x()  const { return m_x; }
  };

}
#endif
