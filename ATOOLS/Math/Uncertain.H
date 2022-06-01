#ifndef ATOOLS_Math_Uncertain_H
#define ATOOLS_Math_Uncertain_H

namespace ATOOLS {

  template <typename T>
  struct Uncertain {
    T value;
    T error;

    int PercentError() const { return int(error / value * 10000) / 100; }
  };

}

#endif
