#ifndef SHERPA_SoftPhysics_Hadron_Init_H
#define SHERPA_SoftPhysics_Hadron_Init_H

#include "ATOOLS/Phys/KF_Table.H"

#include <unordered_set>

namespace SHERPA {

  class Hadron_Init {
  public:
    Hadron_Init() {}
    void Init();

  private:
    void OverrideProperties();

    std::unordered_set<kf_code> m_addedhadrons;
  };

}

#endif
