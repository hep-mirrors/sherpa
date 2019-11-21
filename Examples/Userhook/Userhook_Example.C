#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Main/Sherpa.H"

using namespace ATOOLS;
using namespace SHERPA;

class Userhook_Example : public Userhook_Base {

  Sherpa* p_sherpa;

public:

  Userhook_Example(const Userhook_Arguments args) :
    Userhook_Base("Example"), p_sherpa(args.p_sherpa)
  {
    PRINT_INFO("We are using a user hook within Sherpa and are using PDF "<<p_sherpa->PDFInfo());
  }

  ~Userhook_Example() {}

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs, double &weight) {
    DEBUG_INFO("Let's do something with the bloblist for each event:");

    if(blobs->FourMomentumConservation()) {
      return Return_Value::Nothing;
    }
    else {
      return Return_Value::Error;
    }
  }

  void Finish() {
    PRINT_INFO("Printing something at the end of the run...");
  }

};

DECLARE_GETTER(Userhook_Example,"Example",
	       Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Userhook_Example>::
operator()(const Userhook_Arguments &args) const
{
  return new Userhook_Example(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Userhook_Example>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Example userhook";
}
