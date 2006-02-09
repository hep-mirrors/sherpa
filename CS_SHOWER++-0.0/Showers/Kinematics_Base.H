#ifndef Kinematics_Base_h
#define Kinematics_Base_h

#include "Parton.H"
#include "Vector.H"

namespace CS_SHOWER {
  class Kinematics_Base {
    const static ATOOLS::Vec3D s_ex, s_ey, s_ez;
  public:
    Kinematics_Base() {}
    virtual ~Kinematics_Base() {}
    virtual Parton * MakeKinematics(Parton *, ATOOLS::Flavour &) = 0;
  };

  class Kinematics_FF : public Kinematics_Base {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    Parton * MakeKinematics(Parton *, ATOOLS::Flavour &);
  };
  
  class Kinematics_FI : public Kinematics_Base {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    Parton * MakeKinematics(Parton *, ATOOLS::Flavour &);
  };

  class Kinematics_IF : public Kinematics_Base {
  public:
    Kinematics_IF() {}
    virtual ~Kinematics_IF() {}
    Parton * MakeKinematics(Parton *, ATOOLS::Flavour &);
  };
}

#endif
