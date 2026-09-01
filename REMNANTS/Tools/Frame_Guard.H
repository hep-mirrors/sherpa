#ifndef REMNANTS_Tools_Frame_Guard_H
#define REMNANTS_Tools_Frame_Guard_H

#include "ATOOLS/Math/Poincare.H"
#include <array>
#include <memory>

namespace ATOOLS { class Blob_List; }

namespace REMNANTS {
  class Remnant_Base;

  // RAII bracket around the remnant-breakup construction: boosts the full
  // blob list into the c.m. frame of the two colliding (bunch) particles on
  // construction and back to the lab on destruction, covering all exits
  // including exceptions.  While active, it hands the frame to the remnants
  // through their window pointer, so that their incoming momenta are read in
  // the same frame as the blob contents.  The frame never leaks outside the
  // scope of this guard.
  //
  // The bracket only activates when it can gain phase space: at least one
  // beam with a spectrum (variable energy), a remnant strategy with kinematic
  // content, and a non-trivial boost.  The pair momentum must be
  // longitudinal: the breakup kinematics reconstruct the scatters
  // back-to-back along the z axis, which only holds in frames reached by
  // longitudinal boosts.
  class Frame_Guard {
  private:
    ATOOLS::Blob_List* p_bloblist;
    std::array<std::shared_ptr<Remnant_Base>, 2> p_remnants;
    ATOOLS::Poincare   m_cms;
    bool               m_active;

  public:
    Frame_Guard(ATOOLS::Blob_List* bloblist,
                const std::array<std::shared_ptr<Remnant_Base>, 2>& remnants,
                const bool& kinematics);
    ~Frame_Guard();

    inline bool Active() const { return m_active; }
  };
}

#endif
