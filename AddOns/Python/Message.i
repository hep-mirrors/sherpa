//%module Message
%{
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <ATOOLS/Org/STL_Tools.H>
#include <ATOOLS/Org/CXXFLAGS.H>
#include <ATOOLS/Org/Message.H>
  %}


namespace ATOOLS {

  class Message {
  public:

    // constructor
    Message();

    // destructor
    ~Message(); 

    // member functions
    std::ostream &Error() const;      
    std::ostream &Info() const;       
  };// end of class Message

  %immutable;
  extern Message *msg;
  %mutable;

}
