#ifndef STL_Tools_H
#define STL_Tools_H

#include <iostream>
#include <vector>

namespace std {

  template <typename __Tp>
  std::ostream &operator<<(std::ostream &str,const std::vector<__Tp> &v);

}

#endif

