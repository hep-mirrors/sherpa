#ifndef APACIC_Showers_Jet_Veto_H
#define APACIC_Showers_Jet_Veto_H

#include "ATOOLS/Math/Cluster_Algorithm.H"
#include "AddOns/Apacic++/Showers/Timelike_Kinematics.H"
#include "ATOOLS/Org/Exception.H"
#include "AddOns/Apacic++/Main/Tree.H"
/*!
  \file Jet_Veto.H
  \brief Declares the class Jet_Veto
*/

namespace PHASIC { class Jet_Finder; }

namespace APACIC {

  class Jet_Veto {
  private:

    PHASIC::Jet_Finder  *p_jf;

  public:

    // constructor
    Jet_Veto();

    // destructor
    ~Jet_Veto();

    //member functions 
    int TestFSKinematics(Knot *const mo);
    int TestISKinematics(Knot *const mo,Knot *const partner);

    // inline functions
    inline void SetJetFinder(PHASIC::Jet_Finder *const jf) { p_jf=jf; }

    inline PHASIC::Jet_Finder *const JetFinder() const { return p_jf; }

  };// end of class Jet_Veto

}//end of namespace APACIC

#endif
