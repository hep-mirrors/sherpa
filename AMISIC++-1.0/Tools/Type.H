#ifndef Type_H
#define Type_H

#include "My_String.H"
#include "My_Stringstream.H"
#include "My_Filestream.H"

namespace ATOOLS {

  class Type {
  public:

    enum ID { TChar     =  1,
	      TShort    =  2,
	      TInt      =  3,
	      TLong     =  4,
	      TFloat    =  5,
	      TDouble   =  6,
	      TString   =  7,
	      TIStream  = 11,
	      TOStream  = 12,
	      TIFStream = 13,
	      TOFStream = 14,
	      TFStream  = 15,
	      TUnknown  = 99 };

  public:

    // member functions
    ID GetType();
    ID GetType(const char&);
    ID GetType(const short int&);
    ID GetType(const int&);
    ID GetType(const long int&);
    ID GetType(const float&);
    ID GetType(const double&);
    ID GetType(const std::string&);
    ID GetType(const std::istream&);
    ID GetType(const std::ostream&);
    ID GetType(const std::ifstream&);
    ID GetType(const std::ofstream&);
    ID GetType(const std::fstream&);

  }; // end of class Type

} // end of namespace ATOOLS

#include "Type.IC"

#endif
