#ifndef ATOOLS_Org_Stacktrace_H
#define ATOOLS_Org_Stacktrace_H

#include <iostream>

namespace ATOOLS {

  void GenerateStackTrace(std::ostream &ostr,
			  const bool endline=true,
			  const std::string &comment="");
}

#endif
