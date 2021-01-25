#include "ATOOLS/Phys/Fragmentation_Base.H"

#define COMPILE__Getter_Function
#define OBJECT_TYPE ATOOLS::Fragmentation_Base
#define PARAMETER_TYPE ATOOLS::Fragmentation_Getter_Parameters
#define EXACTMATCH false
#include "ATOOLS/Org/Getter_Function.C"

using namespace ATOOLS;
using namespace std;

Fragmentation_Base::Fragmentation_Base()
{
}
   
Fragmentation_Base::~Fragmentation_Base() 
{
}


namespace ATOOLS {
  class No_Fragmentation : public Fragmentation_Base {
  public:
    No_Fragmentation(const std::string& shower) {}
    ~No_Fragmentation() {}
    
    Return_Value::code Hadronize(ATOOLS::Blob_List *) {
      return Return_Value::Nothing;
    }
  };
}

DEFINE_FRAGMENTATION_GETTER(No_Fragmentation, "None");
