//%module iostream
%{
#include <iostream>
#include <ATOOLS/Org/Exception.H>
  %}

namespace std{
  
  class ostream{
  protected:
    ostream();
    ~ostream();
  public:
    // In order to be able to flush the stream in python, the "ostream" 
    // is extended with the "Flush()" method. Similarly, "Endl()" is added.
    %extend{
    void Flush(){
      *$self<<std::flush;
    }
    void Endl(){
      *$self<<std::endl;
    }

    }
  };
  %rename(Stream_String) &operator<<(std::ostream &str,const char*);
  std::ostream &operator<<(std::ostream &str,const char* );
}

// The declaration of the input operator std::ostream &operator<<(std::ostream &str,const ATOOLS::Exception &exception);
// resides in ATOOLS/Org/Exception.H in Sherpa. 


