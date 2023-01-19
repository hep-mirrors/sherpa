#ifndef PHOTONS_PhotonSplitting_Kinematics_Base_H
#define PHOTONS_PhotonSplitting_Kinematics_Base_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"

namespace PHASIC { class Jet_Finder; }

using namespace ATOOLS;

namespace PHOTONS {
  class Kinematics_FF {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    double GetY(const double Q2,const double kt2,const double z,
		const double mi2,const double mj2,const double mk2) const;
    double GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double mk2) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double mi2,const double mj2,const double mk2) const;
    double GetVirt(const double &Q2,const double &y,const double &z,
      const double mi2,const double mj2,const double mk2) const;
    bool MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double mk2, const double mij2);
  };

  class Kinematics_FI {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    double GetY(const double Q2,const double kt2,const double z,
		const double mi2,const double mj2,const double ma2) const;
    double GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double ma2) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double mi2,const double mj2,const double ma2) const;
    double GetVirt(const double &Q2,const double &y,const double &z,
      const double mi2,const double mj2,const double ma2) const;
    bool MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double ma2, const double mij2);
  };
}

#endif