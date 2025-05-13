#ifndef ATOOLS_Phys_KF_Table_H
#define ATOOLS_Phys_KF_Table_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  void OutputHadrons(std::ostream&);
  void OutputParticles(std::ostream&);
  void OutputContainers(std::ostream&);

  void AddParticle(kf_code kfc, double mass, double radius, double width,
                   int icharge, int strong, int spin, bool majorana, bool on,
                   int stable, bool massive, const std::string &idname,
                   const std::string &antiname, const std::string &texname,
                   const std::string &antitexname);

  // version with no strong charge, and no massive flag, and no anti names
  void AddParticle(kf_code kfc, double mass, double radius, double width,
                   int icharge, int spin, bool majorana, bool on, int stable,
                   const std::string &idname, const std::string &texname);

  typedef std::map<kf_code, Particle_Info *> KFCode_ParticleInfo_Map;

  class KF_Table: public KFCode_ParticleInfo_Map {
  public:

    ~KF_Table();

    kf_code KFFromIDName(const std::string &idname) const;

  };// end of class KF_Table

  extern KF_Table s_kftable;

}

#endif
