#ifndef YFS_Main_Leg_H
#define YFS_Main_Leg_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"

namespace YFS {

  // One leg of a dipole as the form factors see it: an energy and the mass
  // that goes with it.
  // Build from a four-vector where there is one, which cannot violate it. The
  // explicit constructor exists for rest-frame energies, which are computed
  // rather than carried by any vector, and it checks.
  struct Leg {
    double E, m;

    explicit Leg(const ATOOLS::Vec4D &p) : E(p[0]), m(p.Mass()) {}

    Leg(double energy, double mass) : E(energy), m(mass) {
      if (E < m)
        THROW(fatal_error, "YFS::Leg built with E < m: the energy and the mass "
                           "do not describe the same particle.");
    }

    double P2() const { return E*E - m*m; }
  };

}// end of namespace YFS

#endif
