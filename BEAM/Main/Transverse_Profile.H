#ifndef BEAM_Main_Transverse_Profile_H
#define BEAM_Main_Transverse_Profile_H

#include <iosfwd>
#include <utility>
#include <vector>

namespace BEAM {

  enum class Transverse_Profile_Type { Gaussian, UniformSphere, WoodsSaxon };
  std::ostream& operator<<(std::ostream&, Transverse_Profile_Type);
  std::istream& operator>>(std::istream&, Transverse_Profile_Type&);

  /**
   * Normalised, azimuthally symmetric probability distribution of a hard
   * interaction point relative to its beam centre. Distances are in 1/GeV.
   *
   * SampleRadius maps a uniform coordinate to the radial CDF, including the
   * transverse area measure. No density, Jacobian or nuclear A factor is to
   * be multiplied into the event weight. A zero RMS radius is pointlike.
   * This object has no RNG, settings access or event-dependent state.
   */
  class Transverse_Profile {
  private:
    Transverse_Profile_Type m_type;
    double m_rms_radius, m_scale, m_skin_depth;
    // (CDF, s^2 / RMS^2) nodes; empty for analytic or pointlike profiles.
    std::vector<std::pair<double, double>> m_inverse_cdf;

    void BuildWoodsSaxon();

  public:
    // skin_depth is used only by WoodsSaxon; zero is its hard-sphere limit.
    Transverse_Profile(Transverse_Profile_Type type, double rms_radius,
                       double skin_depth = 0.);

    Transverse_Profile_Type Type() const { return m_type; }
    double RMSRadius() const { return m_rms_radius; }
    double SkinDepth() const { return m_skin_depth; }

    /// Inverse radial CDF; requires a finite coordinate in [0,1).
    double SampleRadius(double u) const;
  };

} // namespace BEAM

#endif
