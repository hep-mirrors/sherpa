#ifndef Timelike_Kinematics_H
#define Timelike_Kinematics_H

#include "Knot.H"
#include "Jet_Finder.H"

namespace APACIC {

  class Timelike_Kinematics {
  private:

    double m_pt2min, m_pt_factor;
    int    m_pt_scheme,m_zrange_scheme,m_angle_scheme;
    int    m_type;

    ATOOLS::Jet_Finder *p_jf;

    int ShuffleZ(Knot * const mo) const;
    int ShuffleMoms(Knot * const mo) const; 

    void ConstructTwoVectors(Knot *const mo,ATOOLS::Vec4D &p1vec,
			     ATOOLS::Vec4D &p2vec) const;
    void ConstructDreiBein(Knot *const mo,ATOOLS::Vec3D &n1,
			   ATOOLS::Vec3D &n2) const;
    bool FailedKinCheck(const int first, Knot *const mo) const;
    bool CheckVector(const ATOOLS::Vec4D mom) const;
    void BoostDaughters(ATOOLS::Vec4D pold, ATOOLS::Vec4D pnew, 
			const ATOOLS::Vec4D &pmom, Knot *const mo) const;

  public:

    // constructor
    Timelike_Kinematics(ATOOLS::Jet_Finder *const jf,double pt2min);

    // member functions
    void UpdateDaughters(Knot *const mo,const bool force=false) const;
    bool CheckZRange(Knot const *const mo, 
		     ATOOLS::Flavour const * const d1_flav, 
		     ATOOLS::Flavour const * const d2_flav) const;
    bool CheckZRange(const double z, const double E2, 
		     const double t, const double t1, 
		     const double t2) const;
    double CalcZShift(const double z, const double t, 
		      const double t1, const double t2) const;
    double CalcZShift(const double z, const double t, 
		      const double t1, const double t2, 
		      const double t01, const double t02) const;
    double CalcKt2(const double z, const double E2, 
		   const double t, const double t1, const double t2) const;
    double CalculateAngle(const Knot *);
    double CalculateAngle(const double,const double,const double);

    bool ExtraJetCheck(Knot const * const mo, 
		       Knot const * const d1,Knot const * const d2) const;

    int  Shuffle(Knot * const mo,const int first) const;
    void DoSingleKinematics(Knot * const mo,const bool force=false) const;
    bool DoKinematics(Knot * const mo) const;
    void BoostDaughters(Knot * const mo) const; 

    bool ArrangeColourPartners(ATOOLS::Particle const * const aup, 
			       Knot const * const d1, 
			       Knot const * const d2) const;

    // inline functions
    inline void SetZRangeScheme(const int &zscheme) 
    { m_zrange_scheme=zscheme; }
    inline void SetPTCalcScheme(const int scheme,const double pref) 
    { m_pt_scheme=scheme; m_pt_factor=pref; }
    inline void SetAngleScheme(const int scheme)                    
    { m_angle_scheme=scheme; }

  };// end of class Timelike_Kinematics

  /*!
    \file
    \brief contains the class APACIC::Timelike_Kinematics
  */

  /*!
    \class Timelike_Kinematics
    \brief  Provides all neccessary routines to determine the kinematics of the
            final state shower

    Provides all neccessary routines to determine the kinematics of the
    final state shower, esp.:
     - DoKinematics() and
     - Shuffle()
  */

}

#endif




