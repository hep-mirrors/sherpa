#ifndef POWHEG_Showers_Kinematics_Base_H
#define POWHEG_Showers_Kinematics_Base_H

#include "POWHEG/Tools/Parton.H"
#include "ATOOLS/Math/Vector.H"

namespace PHASIC { class Jet_Finder; }

namespace POWHEG {
  class Sudakov;
  class Kinematics_Base {
  public:
    const static ATOOLS::Vec3D s_ex, s_ey, s_ez;
    PHASIC::Jet_Finder *p_jf;
    ATOOLS::Mass_Selector *p_ms;
    Sudakov *p_sud;
  public:
    Kinematics_Base(): p_jf(NULL), p_ms(NULL) {}
    virtual ~Kinematics_Base() {}
    virtual int MakeKinematics(Parton *const,const ATOOLS::Flavour &,
			       const ATOOLS::Flavour &,Parton *&) = 0;
    virtual double GetY(const double &Q2,const double &kt2,const double &z,
			const double &s1,const double &s2,const double &s3) const = 0;
    double GetS(const double &Q2,const double &y,
		const double &mi2,const double &mj2,const double &mk2) const;
    double GetZ(const double &Q2,const double &sij,const double &y,const double &zt,
		const double &mi2,const double &mk2) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &mi2,const double &mj2,const double &mk2) const;
    double ConstructLN(const double &Q2,const double &sij,
		       const double &mij2,const double &mk2,
		       const ATOOLS::Vec4D &Q,ATOOLS::Vec4D &pk,
		       ATOOLS::Vec4D &l,ATOOLS::Vec4D &n) const;
    inline PHASIC::Jet_Finder *SetJF(PHASIC::Jet_Finder *const jf) 
    { PHASIC::Jet_Finder *ojf(p_jf); p_jf=jf; return ojf; }
    inline PHASIC::Jet_Finder *JF() const { return p_jf; }
    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }
    inline void SetSudakov(Sudakov *const sud) { p_sud=sud; }
  };

  class Kinematics_FF : public Kinematics_Base {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3) const;
  };
  
  class Kinematics_FI : public Kinematics_Base {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3) const;
  };

  class Kinematics_IF : public Kinematics_Base {
  public:
    Kinematics_IF() {}
    virtual ~Kinematics_IF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3) const;
  };

  class Kinematics_II : public Kinematics_Base {
  public:
    Kinematics_II() {}
    virtual ~Kinematics_II() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3) const;
  };
}

#endif
