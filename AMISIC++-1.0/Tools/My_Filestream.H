#ifndef MyFileStream_H
#define MyFileStream_H

#ifdef __GNUC__
#if __GNUC__ > 2
#include <fstream>
#else
#include <fstream.h>
namespace std {

  class ios_base {
  public:

    typedef std::ios::fmtflags fmtflags;

    enum openmode { in  =  1,
		    out =  2 };

  };

}
#endif
#else 
#include <fstream.h>
namespace std {

  class ios_base {
  public:

    typedef std::ios::fmtflags fmtflags;

    enum openmode { in  = 1,
		    out = 2 };

  };

}
#endif

#endif
