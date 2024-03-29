#ifndef MCATNLO_Showers_Kinematics_Base_H
#define MCATNLO_Showers_Kinematics_Base_H

#include "MCATNLO/Tools/Parton.H"
#include "ATOOLS/Math/Vector.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"

namespace PHASIC { class Jet_Finder; }

namespace MCATNLO {
  class Sudakov;
  class Kinematics_Base {
  public:
    const static ATOOLS::Vec3D s_ex, s_ey, s_ez;
    const ATOOLS::Mass_Selector *p_ms;
    Sudakov *p_sud;
    int m_evolscheme;
  public:
    Kinematics_Base(): m_evolscheme(0), p_ms(NULL) {}
    virtual ~Kinematics_Base() {}
    virtual int MakeKinematics(Parton *const,const ATOOLS::Flavour &,
			       const ATOOLS::Flavour &,Parton *&) = 0;
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
    inline void SetMS(const ATOOLS::Mass_Selector *const ms) { p_ms=ms; }
    inline void SetSudakov(Sudakov *const sud) { p_sud=sud; }
    inline void SetEvolScheme(int evol) { m_evolscheme=evol; }
    inline int EvolScheme() const { return m_evolscheme; }
    inline bool ValidateDipoleKinematics(const double & mi2, const double & mj2,
					 const double & mk2, const PHASIC::Kin_Args& args)
    {
      if (ATOOLS::sqr(args.m_pi[0]) < mi2
	  || ATOOLS::sqr(args.m_pj[0]) < mj2
	  || ATOOLS::sqr(args.m_pk[0]) < mk2) {
	msg_IODebugging() << METHOD << "(): Less energy than mass found\n";
	return false;
      }
      return true;
    }
  };

  class Kinematics_FF : public Kinematics_Base {
  public:
    Kinematics_FF() {}
    virtual ~Kinematics_FF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
		const bool force=false) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &mi2,const double &mj2,const double &mk2,
		  const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const;
  };
  
  class Kinematics_FI : public Kinematics_Base {
  public:
    Kinematics_FI() {}
    virtual ~Kinematics_FI() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc,
		const bool force=false) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &mi2,const double &mj2,const double &ma2,
		  const ATOOLS::Flavour &fla,const ATOOLS::Flavour &flc) const;
  };

  class Kinematics_IF : public Kinematics_Base {
  public:
    Kinematics_IF() {}
    virtual ~Kinematics_IF() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		const bool force=false) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &ma2,const double &mi2,const double &mk2,
		  const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const;
  };

  class Kinematics_II : public Kinematics_Base {
  public:
    Kinematics_II() {}
    virtual ~Kinematics_II() {}
    int MakeKinematics(Parton *const,const ATOOLS::Flavour &,const ATOOLS::Flavour &,Parton *&);
    double GetY(const double &Q2,const double &kt2,const double &z,
		const double &s1,const double &s2,const double &s3,
		const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc,
		const bool force=false) const;
    double GetKT2(const double &Q2,const double &y,const double &z,
		  const double &ma2,const double &mi2,const double &mb2,
		  const ATOOLS::Flavour &flb,const ATOOLS::Flavour &flc) const;
  };
}

#endif
