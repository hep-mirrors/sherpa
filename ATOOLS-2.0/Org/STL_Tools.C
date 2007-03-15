#include "STL_Tools.H"

namespace std {

  template <typename __Tp>
  std::ostream &operator<<
    (std::ostream &str,const std::vector<__Tp> &v)
  {
    str<<"(";
    if (v.size()>0) str<<v[0];
    else str<<"<no entry>";
    for (size_t i=1;i<v.size();++i) str<<","<<v[i];
    return str<<")";
  }

  template std::ostream &operator<<
    (std::ostream &str,const std::vector<int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<long int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<unsigned int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<long unsigned int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<float> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<double> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<char> &v);

}
