#include "Sherpa.H"
#include "Message.H"
#include "prof.hh"
#include "Random.H"

extern "C" {
  void apainit_();
  void aparun_();
}

using namespace SHERPA;

int main(int argc,char* argv[]) {  

  Sherpa Generator;
  Generator.InitializeTheRun(std::string("./"));
  AORGTOOLS::msg.Out()<<"generate "<<AORGTOOLS::rpa.gen.NumberOfEvents()<<" events"<<std::endl;
  Generator.InitializeTheEventHandler();

  int nevt=AORGTOOLS::rpa.gen.NumberOfEvents();
  if (nevt>0) AORGTOOLS::msg.Out()<<"Starting event generation now. "<<std::endl;
  for (int i=1;i<=nevt;i++) {
    if (i%500==0) {
      AORGTOOLS::msg.Out()<<" Event "<<i<<endl;      
    }
    if (Generator.GenerateOneEvent()) AORGTOOLS::msg.Events()<<"Sherpa : Passed "<<i<<" events."<<std::endl;
  }
  Generator.SummarizeRun();
  AORGTOOLS::msg.Events()<<"Sherpa did "<<nevt<<" with "<<Generator.NumberOfErrors()<<" errors."<<std::endl;
}



