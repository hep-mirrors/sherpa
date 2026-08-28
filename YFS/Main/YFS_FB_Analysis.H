#ifndef YFS_Main_YFS_FB_Analysis_H
#define YFS_Main_YFS_FB_Analysis_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Weights.H"

#include <string>
#include <vector>

namespace YFS {

  // Conventions for what "forward" and "backward" mean. Each measures
  // cos(theta) of the reference (negatively charged) outgoing particle
  // relative to a reference axis; forward is cos(theta) > 0.
  namespace fbdef {
    enum code {
      lab           = 0,  // polar angle w.r.t. the incoming e- direction
      pair_rest     = 1,  // same, boosted into the outgoing pair rest frame
      collins_soper = 2   // beam-bisector axis in the pair rest frame
    };
  }

  // Short tag used to label the split weights, e.g. "NNLO_CS_F".
  std::string FBTag(fbdef::code def);

  // Splits each named perturbative contribution in a Weights object into
  // forward/backward variants, for one or more angle definitions. The variants
  // are added back into the same Weights object as new named entries; downstream
  // (Single_Process) they pick up the universal ME x PDF x flux factors like any
  // other multiweight, so summing a column gives the forward/backward cross
  // section for that contribution.
  class YFS_FB_Analysis {
  public:
    // enabled: which angle definitions to emit. Empty => all three.
    // target_kf: if non-zero, the reference particle must have this |kf code|
    //            (e.g. 211 pions, 13 muons); otherwise the first negatively
    //            charged outgoing particle is used.
    YFS_FB_Analysis(const std::vector<fbdef::code>& enabled = {},
                    long int target_kf = 0);

    // cos(theta*) of the reference particle under the given definition.
    // Returns a value in [-1,1], or -2 if it cannot be evaluated.
    double CosTheta(fbdef::code def,
                    const ATOOLS::Vec4D_Vector& p,
                    const ATOOLS::Flavour_Vector& fl) const;

    // For every named entry currently in w (the "Nominal" entry is left alone)
    // and every enabled definition, add "<name>_<tag>_F" and "<name>_<tag>_B":
    // one holds the contribution's weight for forward events (0 otherwise) and
    // the other for backward events.
    void SplitWeights(ATOOLS::Weights& w,
                      const ATOOLS::Vec4D_Vector& p,
                      const ATOOLS::Flavour_Vector& fl) const;

    const std::vector<fbdef::code>& Definitions() const { return m_defs; }

  private:
    // Indices of the reference (negative) particle and its positive partner.
    // Either may be -1 if not found.
    void ReferenceIndices(const ATOOLS::Flavour_Vector& fl,
                          int& neg, int& pos) const;

    // Direction of the incoming e- (negatively charged incoming particle),
    // falling back to the +z beam.
    ATOOLS::Vec4D IncomingMinus(const ATOOLS::Vec4D_Vector& p,
                                const ATOOLS::Flavour_Vector& fl) const;

    std::vector<fbdef::code> m_defs;
    long int m_target_kf;
  };

}  // namespace YFS

#endif
