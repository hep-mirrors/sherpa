#include "Primitive_Integrator.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "My_Root.H"
#ifdef ROOT_SUPPORT
#include "TH2D.h"
#endif

using namespace ATOOLS;

class Camel: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point) const 
  {
    const double dx1=0.25, dy1=0.25, w1=1./0.004;
    const double dx2=0.75, dy2=0.75, w2=1./0.004;
    double weight=exp(-w1*((point[0]-dx1)*(point[0]-dx1)+
			   (point[1]-dy1)*(point[1]-dy1)))+
      exp(-w2*((point[0]-dx2)*(point[0]-dx2)+
	       (point[1]-dy2)*(point[1]-dy2)));
#ifdef ROOT_SUPPORT
    static TH2D *ps=NULL, *wt=NULL;
    if (ps==NULL) {
      ps = new TH2D("psg","psg",100,0.,1.,100,0.,1.);
      wt = new TH2D("wtg","wtg",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"psg");
      MYROOT::myroot->AddObject(wt,"wtg");
    }
    ps->Fill(point[0],point[1],1.0);
    wt->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Camel

class Line: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point) const 
  {
    const double w1=1./0.004;
    double weight=exp(-w1*sqr(point[1]+point[0]-1.0));
#ifdef ROOT_SUPPORT
    static TH2D *ps=NULL, *wt=NULL;
    if (ps==NULL) {
      ps = new TH2D("psl","psl",100,0.,1.,100,0.,1.);
      wt = new TH2D("wtl","wtl",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"psl");
      MYROOT::myroot->AddObject(wt,"wtl");
    }
    ps->Fill(point[0],point[1],1.0);
    wt->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Line

class Circle: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point) const 
  {
    const double dx1=0.6, dy1=0.6, rr=0.25, w1=1./0.004, ee=3.0;
    double weight=pow(point[0],ee)*
      exp(-w1*dabs(sqr(point[1]-dy1)+
		   sqr(point[0]-dx1)-sqr(rr)));
    weight+=pow(1.0-point[0],ee)*
      exp(-w1*dabs(sqr(point[1]-1.0+dy1)+
		   sqr(point[0]-1.0+dx1)-sqr(rr)));
#ifdef ROOT_SUPPORT
    static TH2D *ps=NULL, *wt=NULL;
    if (ps==NULL) {
      ps = new TH2D("psc","psc",100,0.,1.,100,0.,1.);
      wt = new TH2D("wtc","wtc",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"psc");
      MYROOT::myroot->AddObject(wt,"wtc");
    }
    ps->Fill(point[0],point[1],1.0);
    wt->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Circle

int main(int argc,char **argv)
{
#ifdef ROOT_SUPPORT
  MYROOT::myroot = new MYROOT::My_Root(argc,argv);
  Exception_Handler::AddTerminatorObject(MYROOT::myroot);
#endif
  std::set_terminate(Exception_Handler::Terminate);
  std::set_unexpected(Exception_Handler::Terminate);
  signal(SIGSEGV,Exception_Handler::SignalHandler);
  signal(SIGINT,Exception_Handler::SignalHandler);
  signal(SIGBUS,Exception_Handler::SignalHandler);
  signal(SIGFPE,Exception_Handler::SignalHandler);
  signal(SIGABRT,Exception_Handler::SignalHandler);
  signal(SIGTERM,Exception_Handler::SignalHandler);
  signal(SIGXCPU,Exception_Handler::SignalHandler);
  try {
    msg.SetLevel(2);
    msg.SetModifiable(true);
    rpa.gen.SetTimeOut(1000);
    PRINT_INFO("Initialize integrator");
    Primitive_Integrator integrator;
    integrator.SetMode(0);
    integrator.SetShuffleMode(1);
    integrator.SetDimension(2);
    integrator.SetNCells(500);
    integrator.SetNOpt(1000);
    integrator.SetNMax(1000000);
    integrator.SetError(1.e-3);
    integrator.Initialize();
    PRINT_INFO("Integrate camel");
    Camel camel;
    integrator.Integrate(&camel);
    integrator.Initialize();
    PRINT_INFO("Integrate line");
    Line line;
    integrator.Integrate(&line);
    integrator.Initialize();
    PRINT_INFO("Integrate circle");
    Circle circle;
    integrator.Integrate(&circle);
#ifdef ROOT_SUPPORT
    MYROOT::myroot->SetDrawOption("lego2");
    delete MYROOT::myroot;
#endif
    return 0;
  }
  catch (Exception exception) {
    exception.UpdateLogFile();
    msg.Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}
