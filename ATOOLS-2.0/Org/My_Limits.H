#ifndef My_Limits_H
#define My_Limits_H

#ifdef __GNUC__
#if __GNUC__ > 2
#include <limits>
#else
namespace std {

  template <> numeric_limits<int>::max()          { return 32767;  }
  template <> numeric_limits<unsigned int>::max() { return 65535;  }
  template <> numeric_limits<long int>::max()     { return 32767;  }
  template <> numeric_limits<float>::max()        { return 3.4e38; }
  template <> numeric_limits<double>::max()       { return 3.4e38; }

}// end of namespace std
#endif
#else
namespace std {

  template <> numeric_limits<int>::max()          { return 32767;  }
  template <> numeric_limits<unsigned int>::max() { return 65535;  }
  template <> numeric_limits<long int>::max()     { return 32767;  }
  template <> numeric_limits<float>::max()        { return 3.4e38; }
  template <> numeric_limits<double>::max()       { return 3.4e38; }

}// end of namespace std
#endif

#endif
