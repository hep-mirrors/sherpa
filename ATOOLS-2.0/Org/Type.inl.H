#ifndef Type_inl_H
#define Type_inl_H

namespace ATOOLS {

  inline Type::ID Type::GetType()                          { return Type::TUnknown;  }
  inline Type::ID Type::GetType(const char& type)          { return Type::TChar;     }
  inline Type::ID Type::GetType(const short int& type)     { return Type::TShort;    }
  inline Type::ID Type::GetType(const int& type)           { return Type::TInt;      }
  inline Type::ID Type::GetType(const unsigned int& type)  { return Type::TUInt;     }
  inline Type::ID Type::GetType(const long int& type)      { return Type::TLong;     }
  inline Type::ID Type::GetType(const float& type)         { return Type::TFloat;    }
  inline Type::ID Type::GetType(const double& type)        { return Type::TDouble;   }
  inline Type::ID Type::GetType(const std::string& type)   { return Type::TString;   }
  inline Type::ID Type::GetType(const std::istream& type)  { return Type::TIStream;  }
  inline Type::ID Type::GetType(const std::ostream& type)  { return Type::TOStream;  }
  inline Type::ID Type::GetType(const std::ifstream& type) { return Type::TIFStream; }
  inline Type::ID Type::GetType(const std::ofstream& type) { return Type::TOFStream; }
  inline Type::ID Type::GetType(const std::fstream& type)  { return Type::TFStream;  }

  inline Type::ID Type::GetType(const ATOOLS::Switch::code code)         { return Type::TSwitch;     }
  inline Type::ID Type::GetType(const ATOOLS::Beam_Type::code code)      { return Type::TBeamType;   }
  inline Type::ID Type::GetType(const ATOOLS::Beam_Shape::code code)     { return Type::TBeamShape;  }
  inline Type::ID Type::GetType(const ATOOLS::Beam_Generator::code code) { return Type::TBeamGen;    }
  inline Type::ID Type::GetType(const ATOOLS::Model_Type::code code)     { return Type::TModelType;  }
  inline Type::ID Type::GetType(const ATOOLS::String_Type::code code)    { return Type::TStringType; }
  inline Type::ID Type::GetType(const ATOOLS::ISR_Type::code code)       { return Type::TISRType;    }

  inline Type::ID Type::GetType(const ATOOLS::Flavour &flavour) { return Type::TFlavour; }

} // end of namespace ATOOLS

#endif
