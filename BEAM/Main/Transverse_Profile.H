#ifndef BEAM_Main_Transverse_Profile_H
#define BEAM_Main_Transverse_Profile_H

#include "ATOOLS/Math/Lookup_Tables.H"

#include <iosfwd>
#include <memory>

namespace BEAM {

  class Woods_Saxon;

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
    double m_rms_radius, m_scale;
    // Woods-Saxon only: the radial CDF as a function of s^2
    std::shared_ptr<const ATOOLS::OneDim_Table> p_cdf;

  public:

    /// Gaussian or projected uniform sphere of the given transverse RMS radius
    Transverse_Profile(Transverse_Profile_Type type, double rms_radius);
    /// Transverse projection of a Woods-Saxon density
    explicit Transverse_Profile(const Woods_Saxon& density);

    Transverse_Profile_Type Type() const { return m_type; }
    double RMSRadius() const { return m_rms_radius; }

    /// Inverse radial CDF; requires a finite coordinate in [0,1).
    double SampleRadius(double u) const;
  };

} // namespace BEAM

#endif
