#ifndef MyStrStream_H
#define MyStrStream_H

// Note: new compilers name "stringstream" in <sstream>
//       instead of "strstream" in <strstream>  !


#if defined __GNUC__
#if __GNUC__ > 2 
// GNU gcc 3.x.x C++ Compiler
#include <sstream>
typedef std::stringstream MyStrStream;
#else
// GNU gcc 2.95.x C++ Compiler
#include <strstream>
// namespace std {
//   typedef strstream stringstream;
// }
typedef std::strstream MyStrStream;
#endif
#endif

#if defined(__sgi) && !defined(__GNUC__)
#include <sstream>
// SGI IRIX C++ Compiler
typedef std::stringstream MyStrStream;
#endif

namespace ATOOLS {
  
  template <class Value_Type>
  std::string ToString(const Value_Type value) {
    MyStrStream converter;
    std::string converted;
    converter<<value;
    converter>>converted;
    return converted;
  }

  template <class Value_Type>
  Value_Type ToType(const std::string &value) {
    MyStrStream converter;
    Value_Type converted;
    converter<<value;
    converter>>converted;
    return converted;
  }

}// end of namespace ATOOLS;

#endif
