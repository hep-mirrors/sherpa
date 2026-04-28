#ifndef BEAM_Ions_Ion_Parameters_H
#define BEAM_Ions_Ion_Parameters_H

#include <map>
#include <string>


namespace BEAM {
  class Ion_Parameters {
  private:
    std::map<std::string, double> m_parametermap;
    std::map<std::string, int>    m_integermap;
  public:
    Ion_Parameters()   {}
    ~Ion_Parameters() {}

    void Initialise();
    const double Get(std::string keyword);
    const int    GetInt(std::string keyword);
  };
  extern Ion_Parameters * ionpars;

}

#endif
