#ifndef CFPShower_Tools_CFP_Parameters_H
#define CFPShower_Tools_CFP_Parameters_H

#include "ATOOLS/Org/Default_Reader.H"

#include <map>
#include <string>

namespace CFPSHOWER {
  class CFP_Parameters {
  private:
    std::map<std::string,int> m_switches;
    std::map<std::string,double> m_parameters;

    void Output();
  public:
    CFP_Parameters();
    ~CFP_Parameters() {}

    bool Init(ATOOLS::Default_Reader *const defaultreader);
    int    operator[](std::string keyword);
    double operator()(std::string keyword);
  };
  
  extern CFP_Parameters * cfp_pars;
}
#endif
