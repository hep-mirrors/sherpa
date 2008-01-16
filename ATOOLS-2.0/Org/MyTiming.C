#include "MyTiming.H"
#include "Message.H"
#include <iostream>
using std::endl;
using namespace ATOOLS;

MyTiming::MyTiming()
{
  clk_tck=sysconf(_SC_CLK_TCK);
  status=3; //never started or stopped 
}

void MyTiming::SetCurrent()
{
  currentclock = times(&currenttms);
}

void MyTiming::Start()
{
  if (status==1) { 
  } 
  else { 
    status=1;
    SetCurrent();
    startclock=currentclock;
    starttms=currenttms;
  }
}

void MyTiming::Stop()
{
  if ((status==0)||(status==3)) { 
  } 
  else { 
    status=0;
    SetCurrent();
    stopclock=currentclock;
    stoptms=currenttms;
  }
}

void MyTiming::PrintTime()
{
  if (status==3) {
  } 
  else {
    if (status==1) SetCurrent();
    double clocks=currentclock-startclock;
    double secs=clocks/clk_tck;
    msg_Info()<<"Time: "<<secs<<" s (clocks="<<clocks<<") on "
	      <<TimeString()<<"\n";
    double utime=(currenttms.tms_utime-starttms.tms_utime)/clk_tck;
    double stime=(currenttms.tms_stime-starttms.tms_stime)/clk_tck;
    double cutime=(currenttms.tms_cutime-starttms.tms_cutime)/clk_tck;
    double cstime=(currenttms.tms_cstime-starttms.tms_cstime)/clk_tck;
    msg_Info()<<" (User: "<<utime<<" s ,System: "<<stime<<" s ,Children User: "
	      <<cutime<<" s ,Children System: "<<cstime<<")\n";
  }
}

double MyTiming::SystemTime()
{
  SetCurrent();
  return (currenttms.tms_stime-starttms.tms_stime)/clk_tck;
}

double MyTiming::UserTime()
{
  SetCurrent();
  return (currenttms.tms_utime-starttms.tms_utime)/clk_tck;
}

double MyTiming::RealTime()
{
  SetCurrent();
  return (currentclock-startclock)/clk_tck;
}

std::string MyTiming::TimeString(const int format)
{
  time_t t(time(NULL));
  std::string tstring(ctime(&t));
  tstring.erase(tstring.length()-1,1);
  for (size_t i(0);i<tstring.length();++i) {
    if ((format&1) && (tstring[i]==' ')) tstring[i]='_';
    if ((format&2) && (tstring[i]==':')) tstring[i]='-';
  }
  return tstring;
}
