%{
#include <string>
%}

namespace ATOOLS {

  class Exception;

  class Terminator_Object {
    
  protected:
    
    virtual bool ReadInStatus(const std::string &path)
    {return true;}
    
    virtual void PrepareTerminate() {}
    
    friend class Terminator_Object_Handler;
    
  public:

    virtual ~Terminator_Object() {}

  };

}

