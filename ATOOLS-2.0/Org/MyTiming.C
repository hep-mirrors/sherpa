#include "MyTiming.H"
#include "Message.H"
#include <iostream>
using std::endl;
using namespace ATOOLS;

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
    //    double clk_tck=double(CLK_TCK);  // does not work properly!
    double clk_tck=100.;
    double secs=clocks/clk_tck;
    msg.Out()<<"Time: "<<secs<<" s   (clocks="<<clocks<<")\n";
    double utime=(currenttms.tms_utime-starttms.tms_utime)/clk_tck;
    double stime=(currenttms.tms_stime-starttms.tms_stime)/clk_tck;
    double cutime=(currenttms.tms_cutime-starttms.tms_cutime)/clk_tck;
    double cstime=(currenttms.tms_cstime-starttms.tms_cstime)/clk_tck;
    msg.Out()<<" (User: "<<utime<<" s ,System: "<<stime<<" s ,Children User: "
	<<cutime<<" s ,Children System: "<<cstime<<")\n";
//     msg.Out()<<" CLK_TCK="<<CLK_TCK<<endl;
//     msg.Out()<<" CLK/s="<<CLOCKS_PER_SEC<<endl;
  }
}

