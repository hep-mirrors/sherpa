#ifndef CSSHOWER_Showers_Kinematics_Base_H
#define CSSHOWER_Showers_Kinematics_Base_H

#include "CSSHOWER++/Tools/Parton.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"

namespace PHASIC { class Jet_Finder; }

namespace CSSHOWER {
  class Sudakov;
  class Kinematics_Base {
  public:
    const ATOOLS::Mass_Selector *p_ms;
    int m_evolscheme;
    int MakeForcedSplittingKinematics(Parton *const split,Parton *& pnew,
				      PHASIC::Kin_Args & args,const int mode=0) const;
  public:
    Kinematics_Base(): m_evolscheme(0), p_ms(NULL) {}
    virtual ~Kinematics_Base() {}
    virtual int MakeKinematics(Parton *const,const double &,
			       const double &,const ATOOLS::Flavour &,
			       Parton *&,const int mode=0) = 0;
    inline void SetMS(const ATOOLS::Mass_Selector *const ms) { p_ms=ms; }
    inline const ATOOLS::Mass_Selector *MS() const { return p_ms; }
    inline void SetEvolScheme(int evol) { m_evolscheme=evol; }
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
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
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
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
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
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
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
    int MakeKinematics(Parton *const,const double &,
		       const double &,const ATOOLS::Flavour &,
		       Parton *&,const int mode=0);
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
