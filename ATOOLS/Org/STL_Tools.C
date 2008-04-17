#include "STL_Tools.H"

#include <string>

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

  template <unsigned char> std::ostream &operator<<
    (std::ostream &str,const std::vector<unsigned char> &v)
  {
    str<<"(";
    if (v.size()>0) str<<v[0];
    else str<<"<no entry>";
    for (size_t i=1;i<v.size();++i) str<<","<<(unsigned short int)(v[i]);
    return str<<")";
  }

  template <char> std::ostream &operator<<
    (std::ostream &str,const std::vector<char> &v)
  {
    str<<"(";
    if (v.size()>0) str<<v[0];
    else str<<"<no entry>";
    for (size_t i=1;i<v.size();++i) str<<","<<(short int)(v[i]);
    return str<<")";
  }

  template std::ostream &operator<<
    (std::ostream &str,const std::vector<short int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<short unsigned int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<unsigned int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<long int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<long unsigned int> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<float> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<double> &v);
  template std::ostream &operator<<
    (std::ostream &str,const std::vector<std::string> &v);

}

namespace ATOOLS {

  bool String_Sort::operator()
    (const std::string &a,const std::string &b) const
  {
    if (a.length()<b.length()) return true;
    if (a.length()>b.length()) return false;
    return m_less(a,b);
  }

}
