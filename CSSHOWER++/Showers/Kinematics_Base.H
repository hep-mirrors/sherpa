#ifndef CSSHOWER_Showers_Kinematics_Base_H
#define CSSHOWER_Showers_Kinematics_Base_H

#include "CSSHOWER++/Tools/Parton.H"
#include "ATOOLS/Math/Vector.H"

namespace PHASIC { class Jet_Finder; }

namespace CSSHOWER {
  class Kinematics_Base {
  public:
    const static ATOOLS::Vec3D s_ex, s_ey, s_ez;
    PHASIC::Jet_Finder *p_jf;
    ATOOLS::Mass_Selector *p_ms;
  public:
    Kinematics_Base(): p_jf(NULL), p_ms(NULL) {}
    virtual ~Kinematics_Base() {}
    virtual int MakeKinematics(Parton *const,const ATOOLS::Flavour &,
			       const ATOOLS::Flavour &,Parton *&) = 0;
    inline PHASIC::Jet_Finder *SetJF(PHASIC::Jet_Finder *const jf) 
    { PHASIC::Jet_Finder *ojf(p_jf); p_jf=jf; return ojf; }
    inline PHASIC::Jet_Finder *JF() const { return p_jf; }
    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }
  };

  class Kinematics_FF : public Kinematics_Base {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
  };
  
  class Kinematics_FI : public Kinematics_Base {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
  };

  class Kinematics_IF : public Kinematics_Base {
  public:
    Kinematics_IF() {}
    virtual ~Kinematics_IF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
  };

  class Kinematics_II : public Kinematics_Base {
  public:
    Kinematics_II() {}
    virtual ~Kinematics_II() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
  };
}

#endif
