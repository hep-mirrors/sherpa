#include "ATOOLS/Phys/EWSudakov_Info.H"

#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

EWSudakov_Log_Type ATOOLS::EWSudakovLogTypeFromString(const std::string& logt)
{
  if (logt == "Ls")
    return EWSudakov_Log_Type::Ls;
  else if (logt == "lZ")
    return EWSudakov_Log_Type::lZ;
  else if (logt == "lSSC")
    return EWSudakov_Log_Type::lSSC;
  else if (logt == "lC")
    return EWSudakov_Log_Type::lC;
  else if (logt == "lYuk")
    return EWSudakov_Log_Type::lYuk;
  else if (logt == "lPR")
    return EWSudakov_Log_Type::lPR;
  else
    THROW(fatal_error,
          "Can not convert " + logt + " to EW Sudakov log type.");
}

std::ostream& ATOOLS::operator<<(std::ostream& os, const EWSudakov_Log_Type& t)
{
  switch (t) {
    case EWSudakov_Log_Type::Ls:
      return os << "LSC";
    case EWSudakov_Log_Type::lZ:
      return os << "Z";
    case EWSudakov_Log_Type::lSSC:
      return os << "SSC";
    case EWSudakov_Log_Type::lC:
      return os << "C";
    case EWSudakov_Log_Type::lYuk:
      return os << "Yuk";
    case EWSudakov_Log_Type::lPR:
      return os << "PR";
  }
}
