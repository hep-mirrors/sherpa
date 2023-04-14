#ifndef ATOOLS_Org_MyStrStream_H
#define ATOOLS_Org_MyStrStream_H

#include "ATOOLS/Org/Exception.H"

// Note: new compilers name "stringstream" in <sstream>
//       instead of "strstream" in <strstream>  !


#include <sstream>
#include <iomanip>
typedef std::stringstream MyStrStream;
typedef std::ios_base::fmtflags MyFmtFlags;
#define IOS_BASE std::ios_base

#include <map>
#include <vector>

namespace ATOOLS {

  typedef std::map<std::string,std::string> String_Map;

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
                             const size_t precision=12,
                             const std::string& separator=" ") {
    MyStrStream converter;
    converter.precision(precision);
    for (typename std::vector<Value_Type>::const_iterator it = values.begin();
         it != values.end();
         ++it) {
      if (it != values.begin()) {
        converter << separator;
      }
      converter << ToString<Value_Type>(*it);
    }
    return converter.str();
  }

  /**
   * this specialization is necessary, because the template code above will
   * truncate the input string after the first whitespace
   */
  template<> std::string ToString<std::string>(const std::string& value,
                                               const size_t precision);


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
    if (converter.fail())
      THROW(fatal_error, "Failed to parse " + value);
    return converted;
  }

  template <class Value_Type>
  std::vector<Value_Type> ToVector(const std::string &values,
                                   const char separator=' ')
  {
    std::vector<Value_Type> converted;
    MyStrStream converter(values);
    std::string token;
    while(std::getline(converter, token, separator)) {
      if (token != "")
        converted.push_back(ToType<Value_Type>(token));
    }
    return converted;
  }

  /**
   * this specialization is necessary, because the template code above will
   * truncate the input string after the first whitespace
   */
  template<> std::string ToType<std::string>(const std::string& value,
                                             const size_t precision);

  /// specialisation to handle conversion of strings like "true" etc.
  template<> bool ToType<bool>(const std::string& value,
                               const size_t precision);

  template<> double ToType<double>(const std::string& value,
                                   const size_t precision);

  std::string StringTrim(const std::string&);

  std::string StringReplace(const std::string &original,
                            const std::string &from, const std::string &to);

  std::string SubString(const std::string& original,
                        const std::string& begin, const std::string& end);


  std::string ReplaceUnits(const std::string& v);

}// end of namespace ATOOLS;

#endif
