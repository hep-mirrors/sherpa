#ifndef MyTiming_h
#define MyTiming_h
#include <sys/times.h>
#include <time.h>
#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

namespace AORGTOOLS {
  class MyTiming {
    tms starttms, currenttms, stoptms;
    clock_t startclock, currentclock, stopclock;
    int status;
    void SetCurrent();
  public:
    MyTiming();
    void Start();
    void Stop();
    void PrintTime();
    ~MyTiming();
  };
}
#endif
