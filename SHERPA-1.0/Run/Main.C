#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"

extern "C" {
  void apainit_();
  void aparun_();
}

using namespace SHERPA;

int main(int argc,char* argv[]) {  
  int nevt=10;

  Sherpa Generator;
  Generator.InitializeTheRun(std::string("./"));
  Generator.InitializeTheEventHandler();
  for (int i=0;i<nevt;i++) {
    if (Generator.GenerateOneEvent()) AORGTOOLS::msg.Events()<<"Sherpa : Passed "<<i<<" events."<<endl;
  }
  AORGTOOLS::msg.Events()<<"Sherpa did "<<nevt<<" with "<<Generator.NumberOfErrors()<<" errors."<<endl;
}



