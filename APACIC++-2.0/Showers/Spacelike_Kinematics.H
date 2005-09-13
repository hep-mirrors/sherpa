#ifndef Spacelike_Kinematics_H
#define Spacelike_Kinematics_H

#include "Timelike_Kinematics.H"
#include "Tree.H"
#include "Knot.H"
#include "Vector.H"
#include "Selector.H"
#include "Poincare.H"
#include "Run_Parameter.H"

namespace APACIC {

  class Spacelike_Kinematics {
  private:
    ATOOLS::Jet_Finder  *p_jf;
    Timelike_Kinematics *p_kin;
    ATOOLS::Poincare     m_boost, m_rot;

    int m_type;

    bool   ResetEnergies(Knot *);
    inline bool CheckVector(const ATOOLS::Vec4D &vec) { return ((vec.Abs2()>0)&&(vec.Abs2()<0)) ? true : false; }

    void RoBoIni(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void RoBoFin(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void BoostPartial(const int mode, Knot * si, const ATOOLS::Vec4D & v_si);
    void BoostPartial(const int mode, Knot * mo, Knot * si, 
		      const ATOOLS::Vec4D & v_mo, const ATOOLS::Vec4D & v_si);

  public:
    inline Spacelike_Kinematics(const double pt2minFS, ATOOLS::Jet_Finder * const jf) : 
      p_jf(jf), p_kin(new Timelike_Kinematics(p_jf,pt2minFS)) {}
    inline ~Spacelike_Kinematics() { if (p_kin) { delete p_kin; p_kin = NULL; }}
    void InitKinematics(Tree ** trees, Knot * k1, Knot * k2,int first);
    bool DoKinematics(Tree **,Knot *, Knot *,int,int first,bool test);
    inline bool JetVeto(const ATOOLS::Vec4D &v) { return (p_jf->TwoJets(v)) ? true : false; } 
    bool JetVeto(Knot *, Knot *);
    static double CalculateMaxT(Knot * k1,Knot * k2);
    double BoostInCMS(Tree **,Knot *, Knot *);
    double BoostFromCMS(Tree **);
    ATOOLS::Vec4D BoostInLab(Tree **);

    void ResetMomenta(Knot*,Tree*,int=0);
    inline void SetJetvetoPt2(const double q2i,const double q2f) { p_jf->SetShowerPt2(q2i); }
  };

  /*!
    \file
    \brief contains the class APACIC::Spacelike_Kinematics
  */

  /*!
    \class Spacelike_Kinematics
    \brief  Provides all neccessary routines to determine the kinematics of the
            initial state shower

    Provides all neccessary routines to determine the kinematics of the
    initial state shower, esp.:
     - InitKinematics() and
     - DoKinematics() 
  */

}

#endif
