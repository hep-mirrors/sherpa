#ifndef ATOOLS_Org_STL_Tools_H
#define ATOOLS_Org_STL_Tools_H

#include <iostream>
#include <vector>

namespace std {

  template <typename __Tp>
  std::ostream &operator<<(std::ostream &str,const std::vector<__Tp> &v);

}

namespace ATOOLS {

  class String_Sort {
  private:

    std::less<std::string> m_less;

  public:

    bool operator()(const std::string &a,const std::string &b) const;
    
  };// end of class String_Sort

  std::vector<int> ID(size_t id);

  size_t IdCount(size_t id);

}

#endif

