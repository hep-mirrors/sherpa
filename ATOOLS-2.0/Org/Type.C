#include "Type.H"

using namespace ATOOLS;

template <> Type::ID Type::GetType(const char &type)          { return Type::TChar;     }
template <> Type::ID Type::GetType(const short int &type)     { return Type::TShort;    }
template <> Type::ID Type::GetType(const int &type)           { return Type::TInt;      }
template <> Type::ID Type::GetType(const unsigned int &type)  { return Type::TUInt;     }
template <> Type::ID Type::GetType(const long int &type)      { return Type::TLong;     }
template <> Type::ID Type::GetType(const float &type)         { return Type::TFloat;    }
template <> Type::ID Type::GetType(const double &type)        { return Type::TDouble;   }
template <> Type::ID Type::GetType(const std::string &type)   { return Type::TString;   }
template <> Type::ID Type::GetType(const std::istream &type)  { return Type::TIStream;  }
template <> Type::ID Type::GetType(const std::ostream &type)  { return Type::TOStream;  }
template <> Type::ID Type::GetType(const std::ifstream &type) { return Type::TIFStream; }
template <> Type::ID Type::GetType(const std::ofstream &type) { return Type::TOFStream; }
template <> Type::ID Type::GetType(const std::fstream &type)  { return Type::TFStream;  }

template <> Type::ID Type::GetType(const ATOOLS::Switch::code &code)         { return Type::TSwitch;     }
template <> Type::ID Type::GetType(const ATOOLS::Beam_Type::code &code)      { return Type::TBeamType;   }
template <> Type::ID Type::GetType(const ATOOLS::Beam_Shape::code &code)     { return Type::TBeamShape;  }
template <> Type::ID Type::GetType(const ATOOLS::Beam_Generator::code &code) { return Type::TBeamGen;    }
template <> Type::ID Type::GetType(const ATOOLS::Model_Type::code &code)     { return Type::TModelType;  }
template <> Type::ID Type::GetType(const ATOOLS::String_Type::code &code)    { return Type::TStringType; }
template <> Type::ID Type::GetType(const ATOOLS::ISR_Type::code &code)       { return Type::TISRType;    }

template <> Type::ID Type::GetType(const ATOOLS::Flavour &flavour) { return Type::TFlavour; }
