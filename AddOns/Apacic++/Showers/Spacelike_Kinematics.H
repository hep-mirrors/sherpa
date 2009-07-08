#ifndef APACIC_Showers_Spacelike_Kinematics_H
#define APACIC_Showers_Spacelike_Kinematics_H

#include "AddOns/Apacic++/Showers/Timelike_Kinematics.H"
#include "AddOns/Apacic++/Main/Tree.H"

namespace APACIC {

  class Spacelike_Kinematics {
  private:

    Timelike_Kinematics *p_kin;

    ATOOLS::Mass_Selector *p_ms;

    ATOOLS::Poincare m_boost, m_rot;

    bool ResetEnergies(Knot *const knot);

    void RoBoIni(Knot *const k,ATOOLS::Poincare &rot,ATOOLS::Poincare &boost);
    void RoBoFin(Knot *const k,ATOOLS::Poincare &rot,ATOOLS::Poincare &boost);

    void BoostPartial(const int mode,Knot *const si,const ATOOLS::Vec4D &v_si);
    void BoostPartial(const int mode,Knot *const mo,Knot *const si, 
		      const ATOOLS::Vec4D &v_mo,const ATOOLS::Vec4D &v_si);

  public:

    // constructor
    Spacelike_Kinematics();
    
    // destructor
    ~Spacelike_Kinematics();

    // member functions
    int InitKinematics(Tree **const trees,Tree *const tree,const int tree1,
		       Knot *const k1, Knot *const k2);
    bool DoKinematics(Tree **const trees,Tree *const tree,const int &leg,
		      Knot *const act,Knot *const part,const bool test);

    double BoostInCMS(Tree **const trees,Tree *const tree,const int tree1,
		      Knot *const k1, Knot *const k2);
    double BoostFromCMS(Tree **const trees,Tree *const tree);

    ATOOLS::Vec4D BoostInLab(Tree **const trees,Tree *const tree);

    void ResetMomenta(Knot *const knot,Tree *const tree,const int &mode=0);

    static double CalculateMaxT(Knot *const k1,Knot *const k2,bool diced=true);

    inline Timelike_Kinematics *Kinematics() const { return p_kin; }

    inline void SetMS(ATOOLS::Mass_Selector *const ms) { p_ms=ms; }

    inline ATOOLS::Mass_Selector *MS() const { return p_ms; }

  };// end of class Spacelike_Kinematics

}// end of namespace APACIC

#endif
