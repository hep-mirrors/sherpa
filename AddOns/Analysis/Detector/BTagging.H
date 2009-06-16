#ifndef Analysis_Detector_BTagging_H
#define Analysis_Detector_BTagging_H

#include "AddOns/Analysis/Detector/Reconstructed_Object.H"

namespace ANALYSIS {
  class BTagging {
    double m_tagprob,m_mistaglight,m_mistagcharm;
    bool m_foundcharm;

    bool FindB(Reconstructed_Object *);
    bool TagIdentified();
    bool Mistag();
  public:
    BTagging(const double tagprob,const double mistag1,const double mistag2);
    ~BTagging();

    bool Tag(Reconstructed_Object *);
  };
}

#endif
