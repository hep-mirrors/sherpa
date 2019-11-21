#ifndef ATOOLS_Phys_KF_Table_H
#define ATOOLS_Phys_KF_Table_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  void OutputHadrons(std::ostream&);
  void OutputParticles(std::ostream&);
  void OutputContainers(std::ostream&);

  class Settings;

  typedef std::map<kf_code,Particle_Info*> KFCode_ParticleInfo_Map;

  class KF_Table: public KFCode_ParticleInfo_Map {
  public:

    ~KF_Table();

    kf_code KFFromIDName(const std::string &idname) const;

  };// end of class KF_Table

  extern KF_Table s_kftable;

}

#endif
