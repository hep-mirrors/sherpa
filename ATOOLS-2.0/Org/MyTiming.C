#include "MyTiming.H"
#include "Message.H"
#include <iostream>
using std::endl;
using namespace AORGTOOLS;

MyTiming::MyTiming()
{
  status=3; //never started or stopped 
  //  starttms=new tms;
  //  currenttms=new tms;
  //  endtms=new tms;
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
    double secs=clocks/double(CLK_TCK);
    msg.Out()<<"Time: "<<secs<<" s   (clocks="<<clocks<<")\n";
    double utime=(currenttms.tms_utime-starttms.tms_utime)/double(CLK_TCK);
    double stime=(currenttms.tms_stime-starttms.tms_stime)/double(CLK_TCK);
    double cutime=(currenttms.tms_cutime-starttms.tms_cutime)/double(CLK_TCK);
    double cstime=(currenttms.tms_cstime-starttms.tms_cstime)/double(CLK_TCK);
    msg.Out()<<" (User: "<<utime<<" s ,System: "<<stime<<" s ,Children User: "
	<<cutime<<" s ,Children System: "<<cstime<<")\n";
  }
}

