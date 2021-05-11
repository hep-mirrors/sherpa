#include "ATOOLS/Phys/EWSudakov_Info.H"

#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

EWSudakov_Log_Type ATOOLS::EWSudakovLogTypeFromString(const std::string& logt)
{
  if (logt == "LSC")
    return EWSudakov_Log_Type::Ls;
  else if (logt == "Z")
    return EWSudakov_Log_Type::lZ;
  else if (logt == "SSC")
    return EWSudakov_Log_Type::lSSC;
  else if (logt == "C")
    return EWSudakov_Log_Type::lC;
  else if (logt == "Yuk")
    return EWSudakov_Log_Type::lYuk;
  else if (logt == "PR")
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

double EWSudakov_Log_Corrections_Map::KFactor() const
{
  double kfac = 1.0;
  for (const auto &kv : *this) {
    kfac += kv.second;
  }
  return kfac;
}

std::ostream& operator<<(std::ostream& o,
                         const EWSudakov_Log_Corrections_Map& m)
{
  o << "1 - K_EWSud = " << m.KFactor() << " (";
  bool is_first {true};
  for (const auto &kv : m) {
    o << kv.first << ": " << (is_first ? "" : ", ") << kv.second;
    is_first = false;
  }
  return o << ')';
}
