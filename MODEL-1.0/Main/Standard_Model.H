#ifndef Standard_Model_h
#define Standard_Model_h

#include "Model_Base.H"

namespace MODEL {
  class Standard_Model : public Model_Base {
  private :
    std::string m_generator;
    int         m_ewscheme,m_ckmorder;

    void ReadInFile();
    void FixEWParameters();
    void FixCKM();
  public :
    Standard_Model(std::string,std::string);
    ~Standard_Model()                   {}
    bool RunSpectrumGenerator();
    void FillMasses() {}
    void FillWidths() {}
  };
}

#endif


