#ifndef Simple_String_H
#define Simple_String_H
/*!
  \file Simple_String.H
  \brief Declares the class Simple_String
*/

#include "MI_Base.H"
#include "Environment.H"

namespace ATOOLS { class Object; }

namespace AMISIC {

  class Reggeon_Trajectory;

  class Simple_String: public MI_Base {
  private:

    ATOOLS::Object *p_remnants[2];

    std::vector<Reggeon_Trajectory*> m_reggeons;

    void CleanUp();
    void PrepareTerminate();

    bool CreateMomenta();

  public:

    // constructor
    Simple_String();
    /*!
      \fn Simple_String()
      \brief The default constructor
    */

    // destructor
    ~Simple_String();

    // member functions
    bool Initialize();
    /*!
      \fn bool Initialize()
      \brief Initalizes Simple_String with the desired
      Reggeon trajectories.
    */
    void Update(const MI_Base *mibase);
    void Reset();

    bool VetoProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool VetoProcess(ATOOLS::Blob *blob)
      \brief Dummy method
    */
    bool DiceProcess();
    /*!
      \fn bool DiceProcess()
      \brief Generates one process according to 
      \f${\rm d}\sigma\,/\,{\rm d}t\f$ employing 
      the class Reggeon_Trajectory.
    */

  }; // end of class Simple_String
  /*!
    \class Simple_String
    \brief Generates soft underlying events according to a given
    Reggeon trajectory.

    This class is employed to generate soft underlying events 
    according to a given Reggeon trajectory. In particular the 
    Pomeron trajectory is used to translate the interaction 
    between two gluons from the initial hadrons into a Reggeon 
    exchange, which yields two final state gluons which are 
    colour-coupled only to their respective beam remnants.
  */

} // end of namespace AMISIC

#endif
