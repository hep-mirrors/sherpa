#ifndef Timelike_Kinematics_H
#define Timelike_Kinematics_H

#include "Knot.H"
#include "Jet_Finder.H"
#include "Run_Parameter.H"

namespace ATOOLS {
  class Data_Read;
}

namespace APACIC {
  class Timelike_Kinematics {
    //! For vetos on extra jet from Shower.
    double                  pt2min,t0;
    bool                    jetveto;
    ATOOLS::Jet_Finder * jf;
    int                     pt_scheme; 
    int                     mass_scheme;
    int      m_type;
    int      m_losejet_veto;
  public:
    //-----------------------------------------------------------------------
    //----------------------- Constructors ----------------------------------
    //----------------------------------------------------------------------- 
    //! Constructor to initialise the Jet_Finder
     Timelike_Kinematics(double _pt2min, ATOOLS::Data_Read * const );
    //! Destructor to delete the Jet_Finder
    ~Timelike_Kinematics() {
      if (jf)   delete jf;
    };
    //-----------------------------------------------------------------------
    //------------------- Checks for kinematics and its evaluation ----------
    //----------------------------------------------------------------------- 
    // ! check if daughts fit and choose one daughter to be diced again if neccessary
    // 1 - nothing changed
    // 0 - selected a daughter to be diced again 
    bool CheckZRange(Knot * mo);
    // check for extra jets and throw veto if neccessary
    bool JetVeto(double,double,double,double,double);

    /*!
      Wrapper for both ShuffleZ and ShuffleMoms. ShuffleZ is called for
      first==0 and ShuffleMoms is called for first==1 (e.g. when FillBranch is
      in used from InitializeJets).
    */
    bool Shuffle(Knot * mo,int first);
    /*! 
      Determines recursively the momenta of partons, using the linked knots
      and the kinematical information within to reconstruct momenta.
      We'll have to redo this for an AO Shower.
    */
    bool DoKinematics(Knot *);

    void SetJetVeto(bool _jetveto) { jetveto = _jetveto; }

    bool ExtraJetCheck(Knot * mo, Knot * d1, Knot * d2);

    bool ArrangeColourPartners(ATOOLS::Particle *,Knot *,Knot *);
    void SetJetvetoPt2(const double pt2) { jf->SetShowerPt2(pt2); }
  private:
    //-----------------------------------------------------------------------
    //------------------- Service for the Checks ----------------------------
    //----------------------------------------------------------------------- 
    /*!
      Given the off-shell masses of the - before on-shell -  daughters,
      the z of the mother, i.e. the energy fractions of the daughters are altered
      in here and the kinematics is checked then.
    */
    bool ShuffleZ(Knot *);
    /*!
      This is like ShuffleZ, only with daughters with known momenta, i.e. from a ME.
      The same checks etc. applies anyhow
    */
    bool ShuffleMoms(Knot *); 
    /*! 
      Check basic kinematic of daughters (1==OK, 0==NO):
      1. Timelike daughters at all ?
      2. Physical opening angle (which is then set) ?
      3. Physical deflection angle ?
      4. New jet from shower (test only if first==1)?
      5. Daughters to be checked ? 
    */
    bool KinCheck(int first,Knot *);
    /*!
      Check individual momenta :
      nan values, negative energies ?
    */
    bool CheckVector(ATOOLS::Vec4D);
    /*!
      Boost/Rotate daughters of given knot, according to difference in pold and pnew.
     */
    void BoostDaughters(Knot * ); 
    void BoostDaughters(ATOOLS::Vec4D pold, ATOOLS::Vec4D pnew, 
			const ATOOLS::Vec4D & pmom, Knot * mo);

  };
}
#endif




