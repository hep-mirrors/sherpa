#ifndef ATOOLS_Phys_KF_Table_H
#define ATOOLS_Phys_KF_Table_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  void OutputHadrons(std::ostream&);
  void OutputParticles(std::ostream&);
  void OutputContainers(std::ostream&);

  void AddParticle(kf_code kfc, double mass, double radius, double width,
                   int icharge, int strong, int spin, int majorana, bool on,
                   int stable, bool massive, const std::string &idname,
                   const std::string &antiname, const std::string &texname,
                   const std::string &antitexname, bool dummy = false,
                   bool group = false);

  void AddOrUpdateParticle(kf_code kfc, double mass, double radius,
                           double width, int icharge, int strong, int spin,
                           int majorana, bool on, int stable, bool massive,
                           const std::string &idname,
                           const std::string &antiname,
                           const std::string &texname,
                           const std::string &antitexname, bool dummy = false,
                           bool group = false);

  void AddParticle(kf_code kfc, double mass, double radius, double width,
                   int icharge, int spin, int majorana, bool on, int stable,
                   const std::string &idname, const std::string &texname);

  void AddParticle(kf_code kfc, double mass, double radius, double width,
                   int icharge, int spin, bool on, int stable,
                   const std::string &idname, const std::string &texname);

  void AddParticle(kf_code kfc, double mass, double radius, int icharge,
                   int spin, int formfactor, const std::string &idname,
                   const std::string &texname);

  void ClearParticles();

  typedef std::map<kf_code, Particle_Info *> KFCode_ParticleInfo_Map;

  class KF_Table: public KFCode_ParticleInfo_Map {
  public:

    ~KF_Table();

    kf_code KFFromIDName(const std::string &idname) const;

  };// end of class KF_Table

  extern KF_Table s_kftable;

}

#endif
