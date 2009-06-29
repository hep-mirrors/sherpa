#ifndef APACIC_Showers_Timelike_Kinematics_H
#define APACIC_Showers_Timelike_Kinematics_H
/*!
  \file Timelike_Kinematics.H
  \brief Declares the class Timelike_Kinematics
*/

#include "ATOOLS/Phys/Blob.H"

namespace APACIC {

  class Knot;

  class Timelike_Kinematics {
  private:

    ATOOLS::Mass_Selector *p_ms;

    int m_zscheme;

    int ShuffleZ(Knot * const mo,const bool update=true) const;
    int ShuffleMomenta(Knot * const mo,const bool update=true,
		       const bool gpm=false) const; 

    void ConstructVectors(Knot *const mo,
			  ATOOLS::Vec4D &p1,ATOOLS::Vec4D &p2) const;
    void ConstructDreiBein(Knot *const mo,
			   ATOOLS::Vec3D &n1,ATOOLS::Vec3D &n2) const;

    bool CheckKinematics(Knot *const mo,const int &first,
			 const bool gpm=false) const;

    inline bool CheckVector(const ATOOLS::Vec4D &mom) const
    { return (mom.Nan() || mom[0]<0.0)?false:true; }

    bool ReconstructDaughter(Knot * const d) const; 

  public:

    // constructor
    Timelike_Kinematics();

    // member functions
    double GetZ(const double &zp,const double &ta, 
		const double &tb,const double &tc, 
		const double &t01=0.0,const double &t02=0.0) const;

    double GetRelativeKT2(const double &z,const double &E2,const double &ta, 
			  const double &t1,const double &t2) const;

    double GetOpeningAngle(Knot *const mo) const;
    double GetOpeningAngle(const double &z,const double &E2,
			   const double &ta,
			   const double &tb,const double &tc) const;
    double GetDeflectionAngle(Knot *const mo) const;
    double GetDeflectionAngle(const double &z,const double &E2,
			      const double &ta,
			      const double &tb,const double &tc) const;

    int  Shuffle(Knot * const mo,const int first) const;

    bool DoSingleKinematics(Knot * const mo,const bool force=false) const;

    bool DoKinematics(Knot * const mo) const;
    bool BoostDaughters(Knot * const mo) const; 
    bool BoostDaughter(Knot * const d) const; 

    int UpdateDaughters(Knot *const mo,const bool force=false) const;

    double LightConeZ(Knot *const knot) const;
    double LightConeZ(const double &z,const double &E2,
		      const double &ta,const double &tb,
		      const double &tc) const;
    double EnergyZ(const double &zlc,const double &E2,
		   const double &ta,const double &tb,
		   const double &tc) const;

    // inline functions
    inline void SetZScheme(const int &scheme) { m_zscheme=scheme; }

    inline int ZScheme() const { return m_zscheme; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

    inline ATOOLS::Mass_Selector *MS() const { return p_ms; }

  };// end of class Timelike_Kinematics
  /*!
    \class Timelike_Kinematics
    \brief Implements procedures to construct final state parton shower 
    kinematics.
    
    This class is employed by Timelike_Sudakov and Final_State_Shower, 
    to construct final state parton shower kinematics and to calculate 
    kinematics variables from the information given in the corresponding 
    Knot.
  */

}// end of namespace APACIC

#endif




