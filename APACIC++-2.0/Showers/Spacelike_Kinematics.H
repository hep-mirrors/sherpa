#ifndef Spacelike_Kinematics_H
#define Spacelike_Kinematics_H

#include "Timelike_Kinematics.H"
#include "Tree.H"

namespace APACIC {

  class Spacelike_Kinematics {
  private:

    ATOOLS::Jet_Finder  *p_jf;
    Timelike_Kinematics *p_kin;

    ATOOLS::Poincare m_boost, m_rot;

    bool ResetEnergies(Knot *const knot);

    void RoBoIni(Knot *const k,ATOOLS::Poincare &rot,ATOOLS::Poincare &boost);
    void RoBoFin(Knot *const k,ATOOLS::Poincare &rot,ATOOLS::Poincare &boost);

    void BoostPartial(const int mode,Knot *const si,const ATOOLS::Vec4D &v_si);
    void BoostPartial(const int mode,Knot *const mo,Knot *const si, 
		      const ATOOLS::Vec4D &v_mo,const ATOOLS::Vec4D &v_si);

  public:

    // constructor
    Spacelike_Kinematics(ATOOLS::Jet_Finder *const jf);
    
    // destructor
    ~Spacelike_Kinematics();

    // member functions
    void InitKinematics(Tree **const trees,Knot *const k1, Knot *const k2,
			const int &first);
    bool DoKinematics(Tree **const trees,Knot *const act,Knot *const part,
		      const int &leg,const bool test);

    double BoostInCMS(Tree **const trees,Knot *const k1, Knot *const k2);
    double BoostFromCMS(Tree **const trees);

    ATOOLS::Vec4D BoostInLab(Tree **const trees);

    void ResetMomenta(Knot *const knot,Tree *const tree,const int &mode=0);

    static double CalculateMaxT(Knot *const k1,Knot *const k2);

    inline Timelike_Kinematics *Kinematics() const { return p_kin; }

  };// end of class Spacelike_Kinematics

}// end of namespace APACIC

#endif
