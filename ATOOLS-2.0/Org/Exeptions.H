#ifndef Exeptions_H
#define Exeptions_H

#include <string>
#include <iostream>

namespace ATOOLS{

  class Exeption {
    std::string mesg;
  public:
    Exeption() { mesg=" default exeption thrown!"; }
    Exeption(std::string error) { mesg=error;}
    void Print() const { 
      std::cout<<mesg<<std::endl;
      std::cerr<<mesg<<std::endl;
    }
  };
  
  class TestExeption : public Exeption {
    TestExeption():Exeption("This is a test exeption") {}; 
  };

} // ATOOLS
#endif // Exeptions_H
