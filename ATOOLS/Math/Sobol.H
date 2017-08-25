#ifndef ATOOLS__Math__Sobol_H
#define ATOOLS__Math__Sobol_H

#include <vector>
#include <string>

namespace ATOOLS {

  class Sobol {
  private:

    unsigned n, D, L;
    
    std::vector<unsigned> d, s, a, X;
    std::vector<std::vector<unsigned> > m, V;
    
    unsigned C(const unsigned i);

  public:

    Sobol(unsigned _D,std::string file="21201");

    std::vector<double> Point();

  };// end of class Sobol

}// end of namespace ATOOLS

#endif
