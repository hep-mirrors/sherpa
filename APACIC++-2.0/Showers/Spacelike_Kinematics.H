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

    bool   ResetEnergies(Knot *);
    inline bool CheckVector(const ATOOLS::Vec4D &vec) 
    { return (vec.Abs2()>0 && vec.Abs2()<0)?true:false; }

    void RoBoIni(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void RoBoFin(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void BoostPartial(const int mode, Knot * si, const ATOOLS::Vec4D & v_si);
    void BoostPartial(const int mode, Knot * mo, Knot * si, 
		      const ATOOLS::Vec4D & v_mo, const ATOOLS::Vec4D & v_si);

  public:

    // constructor
    Spacelike_Kinematics(ATOOLS::Jet_Finder *const jf);
    
    // destructor
    ~Spacelike_Kinematics();

    // member functions
    void InitKinematics(Tree ** trees, Knot * k1, Knot * k2,int first);
    bool DoKinematics(Tree **,Knot *, Knot *,int,int first,bool test);
    bool JetVeto(Knot *, Knot *);
    static double CalculateMaxT(Knot * k1,Knot * k2);
    double BoostInCMS(Tree **,Knot *, Knot *);
    double BoostFromCMS(Tree **);
    ATOOLS::Vec4D BoostInLab(Tree **);

    void ResetMomenta(Knot*,Tree*,int=0);

  };// end of class Spacelike_Kinematics

}// end of namespace APACIC

#endif
