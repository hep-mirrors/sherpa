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
    Timelike_Kinematics *  kink;

    //! For vetos on extra jet from Shower.
    ATOOLS::Jet_Finder* jf;
    ATOOLS::Poincare     boost,rot;
    int m_type;
    /*!
      Starting from a Knot with given energy E all energies of the subsequent
      knots are reset recursively as z*E and (1-z)*E, respectively. 
      In case any of the knots has then t > E^2 a zero will be returned, else
      the returned value is 1.
    */ 
    bool   ResetEnergies(Knot *);
    /*!
      Check individual momenta :
      nan values, negative energies ?
    */
    bool   CheckVector(ATOOLS::Vec4D);



    void RoBoIni(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void RoBoFin(Knot * k, ATOOLS::Poincare & rot, ATOOLS::Poincare & boost);
    void BoostPartial(const int mode, Knot * si, const ATOOLS::Vec4D & v_si);
    void BoostPartial(const int mode, Knot * mo, Knot * si, const ATOOLS::Vec4D & v_mo, const ATOOLS::Vec4D & v_si);
  public:
    //! Constructor to initialise the Jet_Finder
    Spacelike_Kinematics(double pt2minFS, ATOOLS::Data_Read * const );
    //! Destructor to delete the Jet_Finder and the Timelike_Kinematics
    ~Spacelike_Kinematics() {
      if (kink) delete kink;
    }
    /*!
      Initializes kinematics of two incoming spacelike vectors and sets them on the
      z-axis of their cm system.
    */
    void InitKinematics(Tree ** trees, Knot * k1, Knot * k2,int first);
    /*!
      Does the spacelike kinematics in each step. If need be, the timelike
      kinematics of a final state branch is done as well recursively.
      The vectors of the mother and the sister of k1 will be constructed in
      the cm frame of k1 and k2. Both trees are boosted in that system.
    */
    bool DoKinematics(Tree **,Knot *, Knot *,int,int first,bool test);
    /*!
      Checking for extra jets in the initial shower.
      If an emission gives rise to an extra jet, true is returned.
    */
    bool JetCheck(Knot *);
    /*!
      Checking for extra jets in the initial shower.
      If an emission gives rise to an extra jet, true is returned.
    */
    bool JetVeto(const ATOOLS::Vec4D &);
    /*
      Checking for extra jets in the initial shower.
      If an emission gives rise to an extra jet, true is returned.
    */
    bool JetVeto(Knot *, Knot *);
    /*
      Checking for extra jets in the initial shower.
      If an emission gives rise to an extra jet, true is returned.
    */
    bool KinCheck(Knot *,bool);
    /*!
      Knowing the virtual masses of k1 and k2 the maximal allowed timelike t
      of the sister of k1 is calculated in dependence on k1->t, k1->prev->t,
      and k2->t.
    */
    double CalculateMaxT(Knot * k1,Knot * k2);
    /*!
      Boosts two trees related to the two incoming Knots into the
      cm system of these two konts and returns the cm energy squared
      of the two knots.
      Note: Both vectors will be oriented along the z-axis !
    */
    double BoostInCMS(Tree **,Knot *, Knot *);
    double BoostFromCMS(Tree **);

    /*!
      Boosts the initiators of the two trees and thus the trees themselves
      into the lab system and returns the sum of the fourvector of their roots in the 
      lab system.

      Note: Both initiators will be oriented along the z-axis !
    */
    ATOOLS::Vec4D BoostInLab(Tree **);


    void ResetMomenta(Knot*,Tree*,int=0);

    void SetJetvetoPt2(const double q2i,const double q2f) { 
      jf->SetShowerPt2(q2i); 
      kink->SetJetvetoPt2(q2f); 
    }
  };

  /*!
  \fn void Spacelike_Kinematics::ResetMomenta(Knot*,Tree*,int=0);
  \brief ResetsMomenta after jetveto
  */

}
#endif
