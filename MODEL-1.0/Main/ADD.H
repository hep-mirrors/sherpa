#ifndef ADD_h
#define ADD_h

#include "Model_Base.H"

namespace MODEL {
  class ADD : public Model_Base {
  private :
    void   ReadInFile();
  public :
    ADD(std::string,std::string);
    ~ADD()                     {}
    bool   RunSpectrumGenerator(){}
    void   FillMasses()        {}
    void   FillWidths()        {}
  };
}

#endif


