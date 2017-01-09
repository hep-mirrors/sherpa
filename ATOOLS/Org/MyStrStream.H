#ifndef ATOOLS_Org_MyStrStream_H
#define ATOOLS_Org_MyStrStream_H

// Note: new compilers name "stringstream" in <sstream>
//       instead of "strstream" in <strstream>  !


#if defined __GNUC__
#if __GNUC__ > 2 
// GNU gcc 3.x.x C++ Compiler
#include <sstream>
#include <iomanip>
typedef std::stringstream MyStrStream;
typedef std::ios_base::fmtflags MyFmtFlags;
#define IOS_BASE std::ios_base
#else
// GNU gcc 2.95.x C++ Compiler
#include <strstream>
#include <string>
#include <iomanip>
// namespace std {
//   typedef strstream stringstream;
// }
typedef std::strstream MyStrStream;
typedef std::ios::fmtflags MyFmtFlags;
#define IOS_BASE std::ios
#endif
#endif


#if defined(__sgi) && !defined(__GNUC__)
#include <sstream>
#include <iomanip>
// SGI IRIX C++ Compiler
typedef std::stringstream MyStrStream;
typedef std::ios_base::fmtflags MyFmtFlags;
#define IOS_BASE std::ios_base
#endif

#include <vector>

namespace ATOOLS {
  
  template <class Value_Type>
  std::string ToString(const Value_Type &value,
		       const size_t precision=12) {
    MyStrStream converter;
    std::string converted;
    converter.precision(precision);
    converter<<value;
    converter>>converted;
    return converted;
  }

  template <class Value_Type>
  std::string VectorToString(const std::vector<Value_Type> &values,
                             const size_t precision=12) {
    MyStrStream converter;
    converter.precision(precision);
    for (typename std::vector<Value_Type>::const_iterator it = values.begin();
         it != values.end();
         ++it) {
      if (it != values.begin()) {
        converter << " ";
      }
      converter << ToString<Value_Type>(*it);
    }
    return converter.str();
  }

  template <class Value_Type>
  std::string MatrixToString(const std::vector<std::vector<Value_Type> > &values,
                             const size_t precision=12) {
    MyStrStream converter;
    converter.precision(precision);
    for (typename std::vector<std::vector<Value_Type> >::const_iterator it = values.begin();
         it != values.end();
         ++it) {
      if (it != values.begin()) {
        converter << "\n";
      }
      converter << VectorToString<Value_Type>(*it);
    }
    return converter.str();
  }

  template <class Value_Type>
  Value_Type ToType(const std::string &value,
		    const size_t precision=12) {
    MyStrStream converter;
    Value_Type converted;
    converter.precision(precision);
    converter<<value;
    converter>>converted;
    return converted;
  }

  std::string StringReplace(const std::string &original,
                            const std::string &from, const std::string &to);

  std::string ReplaceUnits(std::string v);

}// end of namespace ATOOLS;

#endif
