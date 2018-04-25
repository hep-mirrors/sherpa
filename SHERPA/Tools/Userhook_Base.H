#ifndef SHERPA_Tools_Userhook_Base_H
#define SHERPA_Tools_Userhook_Base_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Return_Value.H"
#include <string>

namespace SHERPA {

  class Sherpa;

  struct Userhook_Arguments {
    Sherpa* p_sherpa;
    Userhook_Arguments(Sherpa* const sherpa):
      p_sherpa(sherpa) {}
  };

  class Userhook_Base {
  protected:

    std::string m_name;

  public:

    typedef ATOOLS::Getter_Function
    <Userhook_Base,Userhook_Arguments> Getter_Function;

  public:

    /**
     * Constructor
     *
     * Any classes deriving from Userhook_Base should use a constructor like:
     *   Userhook_Example(const Userhook_Arguments args) :
     *     Userhook_Base("Example")
     *     { ... }
     */
    Userhook_Base(const std::string &name);

    /**
     * Destructor
     */
    virtual ~Userhook_Base();

    /**
     * This is where the action happens. Allowed return values are:
     * - Return_Value::Nothing if the event was not changed.
     * - Return_Value::Success if the event was successfully changed.
     * - Return_Value::New_Event if the event is to be discarded (and NTrials increased)
     * - Return_Value::Error if an error happened
     */
    virtual ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs, double &weight) = 0;

    /**
     * Optional function which runs at the end of the run.
     */
    virtual void Finish();

    /**
     * Optional function which runs when an event is being retried.
     */
    virtual void CleanUp();

    inline const std::string Name() { return m_name; }
  };

  typedef std::vector<Userhook_Base*> Userhook_Vector;

}

#endif
