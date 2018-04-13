%include <std_string.i>
%{
#include <ATOOLS/Org/Message.H>
#include <ATOOLS/Org/Exception.H>
#include <ATOOLS/Org/MyStrStream.H>
#include <iostream>
%}

namespace ATOOLS {

  %rename(SherpaException) Exception;

  class Exception {
  public:
    %extend {
      PyObject* __str__() {
	MyStrStream conv;
	conv<<*self;
	return PyString_FromString(conv.str().c_str());
      };
    };

  };
  
}
