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
    template <class ValueType>
    static ID GetType(const ValueType &) { return Type::TUnknown; }

  }; // end of class Type

} // end of namespace ATOOLS

#endif
