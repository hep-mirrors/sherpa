#ifndef Hadron_Decay_Handler_h
#define Hadron_Decay_Handler_h

#include <set>
#include <iterator>

#include "CXXFLAGS.H"

#include "Return_Value.H"
#include "Flavour.H"

#ifdef USING__Hadrons
namespace HADRONS { class Hadrons; }
#endif
namespace ATOOLS { class Blob; class Particle; class Vec4D; }

namespace SHERPA {
  class Lund_Interface;

  class Hadron_Decay_Handler {
  protected :
    std::string                   m_decmodel;
    int                           m_mode;
    bool                          m_mass_smearing;
#ifdef USING__Hadrons
    HADRONS::Hadrons            * p_hadrons;
#endif
    Lund_Interface              * p_lund; 
    std::set<ATOOLS::kf::code>  * p_cans;

  public :
#ifdef USING__Hadrons
    Hadron_Decay_Handler(HADRONS::Hadrons *);
#endif
    Hadron_Decay_Handler(Lund_Interface *);

    ~Hadron_Decay_Handler();

    ATOOLS::Return_Value::code FillHadronDecayBlob(ATOOLS::Blob *, const ATOOLS::Vec4D&);
    bool CanDealWith(ATOOLS::kf::code);

#ifdef USING__Hadrons
    HADRONS::Hadrons * GetHadrons() { return p_hadrons; }
#endif
    Lund_Interface   * GetLund()    { return p_lund; }
    std::string        Name()       { return m_decmodel; }
    std::set<ATOOLS::kf::code>  * GetCans() { return p_cans; }
    void               SetMassSmearing(bool ms)      { m_mass_smearing = ms; }
    bool               GetMassSmearing()             { return m_mass_smearing; }
    bool               DiceMass(ATOOLS::Particle*, double min, double max);
    void               SetSignalProcessBlob(ATOOLS::Blob* spblob);
  };

  typedef std::map<std::string,Hadron_Decay_Handler *> HDHandlersMap;
  typedef HDHandlersMap::iterator                      HDHandlersIter;
}

  /*!
    \file Hadron_Decay_Handler.H
    \brief Declares the class SHERPA::Hadron_Decay_Handler
    */
  /*!
    \class Hadron_Decay_Handler
    \brief A wrapper class around all implemented hadron decay modules.

    Knows all different hadron decay handlers and can tell whether a specific
    one can deal with an unstable particle (CanDealWith()). It can then pass
    the to be decayed blob to it (FillHadronDecayBlob()).
    */
  /*!
    \fn FillHadronDecayBlob(ATOOLS::Blob *)
    \brief Passes the to be decayed blob to the hadron decay handler
    incorporated.

    It calls the decay method of its hadron decay handler for the given blob,
    and returns the \c Return_Value::code that it receives from it.
    */
  /*!
    \fn EraseTreated(std::set<int> *)
    \brief Erases all kfcodes, which the HADRONS module is
    able to decay, from the given set

    This is used for finding out whether we need an additional decay module
    (like Lund).
    */

#endif
