#ifndef MyStringStream_H
#define MyStringStream_H

#ifdef __GNUC__
#if __GNUC__ > 2
#include <sstream>
#else
#include <strstream.h>
typedef std::strstream std::stringstream;
#endif
#endif

#endif
