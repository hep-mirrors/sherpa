#include "ATOOLS/Org/MyTiming.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"
#include <iostream>
#include <unistd.h>
#include <time.h>

using std::endl;
using namespace ATOOLS;

namespace ATOOLS {
std::string FormatTime(const int &intime,const int mode)
{
  size_t in = (intime<0 ? 0 : intime);
  int days(in/86400), hrs((in%86400)/3600);
  int mins(((in%86400)%3600)/60), secs(((in%86400)%3600)%60);
  std::string out;
  if (days) out+=ToString(days)+"d ";
  if (hrs) out+=ToString(hrs)+"h ";
  if (mins) out+=ToString(mins)+"m ";
  if (secs) out+=ToString(secs)+"s ";
  if (out.length()) out.erase(out.length()-1,1);
  else out="0s";
  return out;
}
}

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
    msg_Info()<<"Time: "<<FormatTime((size_t)secs)<<" on "
	      <<TimeString()<<"\n";
    double utime=(currenttms.tms_utime-starttms.tms_utime)/clk_tck;
    double stime=(currenttms.tms_stime-starttms.tms_stime)/clk_tck;
    double cutime=(currenttms.tms_cutime-starttms.tms_cutime)/clk_tck;
    double cstime=(currenttms.tms_cstime-starttms.tms_cstime)/clk_tck;
    msg_Info()<<" (User: "<<FormatTime((size_t)utime)<<", System: "
	      <<FormatTime((size_t)stime)<<", Children User: "
	      <<FormatTime((size_t)cutime)<<", Children System: "
	      <<FormatTime((size_t)cstime)<<")\n";
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

std::string MyTiming::StrFTime
(const std::string &format,const time_t &offset)
{
  time_t t(time(NULL));
  t+=offset;
  std::string rv(100,' ');
  if (strftime(&rv[0],rv.length(),format.c_str(),localtime(&t))==0) {
    msg_Error()<<METHOD<<"(): Error converting time string."<<std::endl;
    return "";
  }
  while (rv[0]==' ') rv.erase(0,1);
  while (rv[rv.length()-1]==' ') rv.erase(rv.length()-1,1);
  return rv;
}
