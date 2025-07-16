#include "MODEL/UFO/UFO_Color_Functions.H"
#include "ATOOLS/Org/Exception.H"
#include <algorithm>

using namespace MODEL;
using std::stringstream;

namespace UFO{

  void UFO_Color_Function::SetNext(UFO_Color_Function* next) {
    if (p_next) delete p_next;
    p_next = next;
  }

  UFO_Color_Function UFO_Color_Function::operator * (const UFO_Color_Function & other) const {

    UFO_Color_Function ret(*this);
    UFO_Color_Function* next = new UFO_Color_Function(other);
    // ownership transferred here
    ret.SetNext(next);
    return ret;
  }

}
