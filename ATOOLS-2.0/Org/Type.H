#ifndef Type_H
#define Type_H

#include "Data_Read.H"
#include "Flavour.H"
#include "MyStrStream.H"
#include <vector>
#include <fstream>

namespace ATOOLS {

  class Type {
  public:

    enum ID { TChar       =  1,
	      TShort      =  2,
	      TInt        =  3,
	      TUInt       =  4,
	      TLong       =  5,
	      TFloat      =  6,
	      TDouble     =  7,
	      TString     =  8,
	      TIStream    = 11,
	      TOStream    = 12,
	      TIFStream   = 13,
	      TOFStream   = 14,
	      TFStream    = 15,
	      TSwitch     = 16,
	      TBeamType   = 17,
	      TBeamShape  = 18,
	      TBeamGen    = 19,
	      TStringType = 20,
	      TModelType  = 21,
	      TISRType    = 22,
	      TFlavour    = 23,
 	      TUnknown    = 99 };

  public:

    // member functions
    static ID GetType();
    static ID GetType(const char&);
    static ID GetType(const short int&);
    static ID GetType(const int&);
    static ID GetType(const unsigned int&);
    static ID GetType(const long int&);
    static ID GetType(const float&);
    static ID GetType(const double&);
    static ID GetType(const std::string&);
    static ID GetType(const std::istream&);
    static ID GetType(const std::ostream&);
    static ID GetType(const std::ifstream&);
    static ID GetType(const std::ofstream&);
    static ID GetType(const std::fstream&);
    static ID GetType(const Switch::code);
    static ID GetType(const Beam_Type::code);
    static ID GetType(const Beam_Shape::code);
    static ID GetType(const Beam_Generator::code);
    static ID GetType(const String_Type::code);
    static ID GetType(const Model_Type::code);
    static ID GetType(const ISR_Type::code);
    static ID GetType(const Flavour &);

  }; // end of class Type

} // end of namespace ATOOLS

#include "Type.inl.H"

#endif
