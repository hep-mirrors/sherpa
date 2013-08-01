#ifndef CSSHOWER_Showers_Kinematics_Base_H
#define CSSHOWER_Showers_Kinematics_Base_H

#include "CSSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"

namespace PHASIC { class Jet_Finder; }

namespace CSSHOWER {
  class Sudakov;
  class Kinematics_Base {
  public:
    ATOOLS::Mass_Selector *p_ms;
  public:
    Kinematics_Base() {}
    virtual ~Kinematics_Base() {}
    virtual int MakeKinematics(Parton *const,const double &,
			       const double &,const ATOOLS::Flavour &,
			       Parton *&,const int mode=0) = 0;
    void SetFixVec(Parton *const p,ATOOLS::Vec4D mom,
		   const PHASIC::Kin_Args &lt,const int mode) const;
    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }
    inline ATOOLS::Mass_Selector *MS() const { return p_ms; }
  };

  class Kinematics_FF : public Kinematics_Base {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flb,
		const bool force=false) const;
  };
  
  class Kinematics_FI : public Kinematics_Base {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flb,
		const bool force=false) const;
  };

  class Kinematics_IF : public Kinematics_Base {
  public:
    Kinematics_IF() {}
    virtual ~Kinematics_IF() {}
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const bool force=false) const;
  };

  class Kinematics_II : public Kinematics_Base {
  public:
    Kinematics_II() {}
    virtual ~Kinematics_II() {}
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const bool force=false) const;
  };
}

#endif
