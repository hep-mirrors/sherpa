#ifndef MyString_H
#define MyString_H

#ifdef __GNUC__
#include <string>
#else
#include <string>
#endif

namespace ATOOLS {

  const std::string nullstring = std::string("");

}

#endif
