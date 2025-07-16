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
    void InitHadron(const kf_code&,
                    const bool checkinitialised,
                    const double &mass,
                    const double &radius,
                    const double &width,
                    const int icharge,
                    const int spin,
                    const bool majorana,
                    const bool on,
                    const int stable,
                    const std::string& idname,
                    const std::string& texname);
    void InitHadron(const kf_code&,
                    const double &mass,
                    const double &radius,
                    const double &width,
                    const int icharge,
                    const int strong,
                    const int spin,
                    const bool majorana,
                    const bool on,
                    const int stable,
                    const bool massive,
                    const std::string& idname,
                    const std::string& antiname,
                    const std::string& texname,
                    const std::string& antitexname);
    void OverrideProperties();

    std::unordered_set<kf_code> m_addedhadrons;
  };

}

#endif
