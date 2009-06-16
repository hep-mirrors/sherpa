#ifndef ATOOLS_Org_My_Limits_H
#define ATOOLS_Org_My_Limits_H

#ifdef __GNUC__
#if __GNUC__ < 3
#else
#define USE_LIMITS
#endif
#endif

#if defined(__sgi) && !defined(__GNUC__)
#define USE_LIMITS
#endif

#ifdef USE_LIMITS
#include <limits>
#else
namespace std {

  template <typename _Tp>
  struct numeric_limits {

    static _Tp max() { return static_cast<_Tp>(0); }

  };// end of class numeric_limits

  template <> int          numeric_limits<int>::max()          { return 32767;  }
  template <> unsigned int numeric_limits<unsigned int>::max() { return 65535;  }
  template <> long int     numeric_limits<long int>::max()     { return 32767;  }
  template <> float        numeric_limits<float>::max()        { return 3.4e38; }
  template <> double       numeric_limits<double>::max()       { return 3.4e38; }

}// end of namespace std
#endif

#endif
