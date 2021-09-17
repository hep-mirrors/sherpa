#ifndef CFPShower_Tools_CFP_Parameters_H
#define CFPShower_Tools_CFP_Parameters_H

#include <map>
#include <string>

namespace CFPSHOWER {
  class CFP_Parameters {
  private:
    std::map<std::string,int>    m_switches;
    std::map<std::string,double> m_parameters;
  public:
    CFP_Parameters();
    ~CFP_Parameters() {}

    bool Init();
    int    operator[](std::string keyword);
    double operator()(std::string keyword);
    
    void Output();
  };
  
  extern CFP_Parameters * cfp_pars;
}
#endif
