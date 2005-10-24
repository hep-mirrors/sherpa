#ifndef Timelike_Kinematics_H
#define Timelike_Kinematics_H

#include "Jet_Finder.H"

namespace APACIC {

  class Knot;

  class Timelike_Kinematics {
  private:

    ATOOLS::Jet_Finder *p_jf;

    int m_zscheme, m_anglescheme;

    int ShuffleZ(Knot * const mo) const;
    int ShuffleMomenta(Knot * const mo) const; 

    void ConstructVectors(Knot *const mo,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2) const;
    void ConstructDreiBein(Knot *const mo,
			   ATOOLS::Vec3D &n1,ATOOLS::Vec3D &n2) const;

    bool CheckKinematics(Knot *const mo,const int &first) const;
    bool CheckVector(const ATOOLS::Vec4D &mom) const;

  public:

    // constructor
    Timelike_Kinematics(ATOOLS::Jet_Finder *const jf);

    // member functions
    bool CheckZRange(Knot *const mo, 
		     const ATOOLS::Flavour *const d1flav, 
		     const ATOOLS::Flavour *const d2flav) const;
    bool CheckZRange(const double &z,const double &E2,const double &ta, 
		     const double &tb,const double &tc) const;

    double GetZ(const double &zp, const double &ta, 
		const double &tb, const double &tc) const;
    double GetZ(const double &zp, const double &ta, 
		const double &tb, const double &tc, 
		const double &t01, const double &t02) const;

    double GetRelativeKT2(const double &z,const double &E2,const double &ta, 
			  const double &t1,const double &t2) const;

    double GetOpeningAngle(Knot *const mo) const;
    double GetOpeningAngle(const double &z,const double &E2,
			   const double &ta,
			   const double &tb, const double &tc) const;
    double GetDeflectionAngle(Knot *const mo) const;
    double GetDeflectionAngle(const double &z,const double &E2,
			      const double &ta,
			      const double &tb, const double &tc) const;

    int  Shuffle(Knot * const mo,const int first) const;

    void DoSingleKinematics(Knot * const mo,const bool force=false) const;

    bool DoKinematics(Knot * const mo) const;
    void BoostDaughters(Knot * const mo) const; 

    int UpdateDaughters(Knot *const mo,const bool force=false) const;


    bool ArrangeColourPartners(ATOOLS::Particle *const aup, 
			       Knot *const d1,Knot *const d2) const;

    double LightConeZ(Knot *const knot) const;
    double LightConeZ(const double &z,const double &E2,
		      const double &ta,const double &tb,
		      const double &tc) const;

    // inline functions
    inline void SetZScheme(const int &zscheme)    { m_zscheme=zscheme;    }
    inline void SetAngleScheme(const int &scheme) { m_anglescheme=scheme; }

    inline int ZScheme() const     { return m_zscheme;     }
    inline int AngleScheme() const { return m_anglescheme; }

  };// end of class Timelike_Kinematics

}// end of namespace APACIC

#endif




