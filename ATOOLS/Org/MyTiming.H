#ifndef ATOOLS_Org_MyTiming_H
#define ATOOLS_Org_MyTiming_H
#include <sys/times.h>
#include <time.h>
#include <string>

#ifndef CLK_TCK
#define CLK_TCK CLOCKS_PER_SEC
#endif

namespace ATOOLS {
  class MyTiming {
    tms starttms, currenttms, stoptms;
    clock_t startclock, currentclock, stopclock;
    int status;
    double clk_tck;
    void SetCurrent();
  public:
    MyTiming();
    void Start();
    void Stop();
    void PrintTime();
    double SystemTime();
    double UserTime();
    double RealTime();
    std::string TimeString(const int format=0);
    std::string StrFTime(const std::string &format,
			 const time_t &offset=0);
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
