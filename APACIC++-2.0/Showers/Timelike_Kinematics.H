#ifndef Timelike_Kinematics_H
#define Timelike_Kinematics_H

#include "Knot.H"
#include "Jet_Finder.H"

namespace ATOOLS {
  class Data_Read;
}

namespace APACIC {
  class Timelike_Kinematics {
    double m_pt2min;
    bool   m_jetveto, m_nolosejetveto;
    int    m_pt_scheme; 
    int    m_zrange_scheme; 
    int    m_mass_scheme;
    int    m_type;
    int    m_losejet_veto;

    ATOOLS::Jet_Finder * p_jf;
  public:
    //-----------------------------------------------------------------------
    //----------------------- Constructors ----------------------------------
    //----------------------------------------------------------------------- 
    //! Constructor to initialise the Jet_Finder
     Timelike_Kinematics(double pt2min, ATOOLS::Data_Read * const dataread);
    //! Destructor to delete the Jet_Finder
    ~Timelike_Kinematics() {
      if (p_jf)   delete p_jf;
    };
    //-----------------------------------------------------------------------
    //------------------- Checks for kinematics and its evaluation ----------
    //----------------------------------------------------------------------- 
    // ! check if daughts fit and choose one daughter to be diced again if neccessary
    // 1 - nothing changed
    // 0 - selected a daughter to be diced again 
    bool CheckZRange(Knot const * const mo, ATOOLS::Flavour const * const d1_flav, 
		     ATOOLS::Flavour const * const d2_flav) const;
    // check for extra jets and throw veto if neccessary
    bool JetVeto(double t, double e2, double z, double t1, double t2);

    /*!
      Wrapper for both ShuffleZ and ShuffleMoms. ShuffleZ is called for
      first==0 and ShuffleMoms is called for first==1 (e.g. when FillBranch is
      in used from InitializeJets).
    */
    int Shuffle(Knot * const mo,const int first) const;
    /*! 
      Determines recursively the momenta of partons, using the linked knots
      and the kinematical information within to reconstruct momenta.
      We'll have to redo this for an AO Shower.
    */
    bool DoKinematics(Knot * const mo) const;

    void SetJetVeto(bool jetveto) { m_jetveto = jetveto; }
    void SetNoLoseJetVeto(bool nolosejetveto) { m_nolosejetveto = nolosejetveto; }

    bool ExtraJetCheck(Knot const * const mo, Knot const * const d1, Knot const * const d2) const;

    bool ArrangeColourPartners(ATOOLS::Particle const * const aup, Knot const * const d1, Knot const * const d2) const;
    void SetJetvetoPt2(const double pt2) const { p_jf->SetShowerPt2(pt2); }

    static bool   CheckZRange(double z, double E2, double t, double t1, double t2);
    static double CalcZShift(double z, double t, double t1, double t2, double t01=0., double t02=0.);
    static double CalcKt2(double z, double E2, double t, double t1, double t2);
  private:
    //-----------------------------------------------------------------------
    //------------------- Service for the Checks ----------------------------
    //----------------------------------------------------------------------- 
    /*!
      Given the off-shell masses of the - before on-shell -  daughters,
      the z of the mother, i.e. the energy fractions of the daughters are altered
      in here and the kinematics is checked then.
    */
    int ShuffleZ(Knot * const mo) const;
    /*!
      This is like ShuffleZ, only with daughters with known momenta, i.e. from a ME.
      The same checks etc. applies anyhow
    */
    int ShuffleMoms(Knot * const mo) const; 
    /*! 
      Check basic kinematic of daughters (1==OK, 0==NO):
       -# Timelike daughters at all ?
       -# Physical opening angle (which is then set) ?
       -# Physical deflection angle ?
       -# New jet from shower (test only if first==1)?
       -# Daughters to be checked ? 
    */
    int KinCheck(const int first, Knot * const mo) const;
    /*!
      Check individual momenta :
      nan values, negative energies ?
    */
    bool CheckVector(const ATOOLS::Vec4D mom) const;
    /*!
      Boost/Rotate daughters of given knot, according to difference in pold and pnew.
     */
    void BoostDaughters(Knot * const mo) const; 
    void BoostDaughters(ATOOLS::Vec4D pold, ATOOLS::Vec4D pnew, 
			const ATOOLS::Vec4D & pmom, Knot * const mo) const;

  };


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




