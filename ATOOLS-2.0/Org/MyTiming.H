#ifndef MyTiming_H
#define MyTiming_H
#include <sys/times.h>
#include <time.h>

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

namespace ATOOLS {
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
  };

  /*!
    \file
    \brief Contains the class MyTiming
  */

  /*!
    \class MyTiming
    \brief Simple class for making time measurements
  */

  /*!
    \fn MyTiming::MyTiming()
    \brief Standard constructor
  */

  /*!
    \fn void MyTiming::Start();
    \brief start watch
  */

  /*!
    \fn void MyTiming::Stop();
    \brief stop watch
  */

  /*!
    \fn void MyTiming::PrintTime();
    \brief print results on screen
  */


}
#endif
