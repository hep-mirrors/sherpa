#include "Primitive_Integrator.H"
#include "Exception.H"
#include "My_Root.H"
#ifdef ROOT_SUPPORT
#include "TH2D.h"
#endif

using namespace ATOOLS;

#ifdef ROOT_SUPPORT
class PS_Histogram: public TH2D {
private:
  TH2D *p_ps;
public:
  // constructor
  PS_Histogram(const char *name,const char *title,
	       const size_t nbinsx,const double xmin,const double xmax,
	       const size_t nbinsy,const double ymin,const double ymax):
    TH2D(name,title,nbinsx,xmin,xmax,nbinsy,ymin,ymax),
    p_ps(new TH2D((std::string(name)+"_ps").c_str(),
		  (std::string(title)+"_ps").c_str(),
		  nbinsx,xmin,xmax,nbinsy,ymin,ymax)) {}
  // destructor
  ~PS_Histogram() {}
  // member functions
  Int_t Fill(const Double_t x,const Double_t y,
	     const Double_t weight)
  {
    p_ps->Fill(x,y,1.0);
    return TH2D::Fill(x,y,weight);
  }
  void Draw(Option_t *option="")
  {
    for (Int_t i=0;i<GetNbinsX();++i)
      for (Int_t j=0;j<GetNbinsY();++j)
	SetBinContent(i,j,p_ps->GetBinContent(i,j)==0.0?0.0:
		      GetBinContent(i,j)/
		      p_ps->GetBinContent(i,j));
#ifdef USING__Distinct_Canvas
    TH2D::Draw(option);
    Int_t logx=gPad->GetLogx();
    Int_t logy=gPad->GetLogy();
    Int_t logz=gPad->GetLogz();
    TCanvas *psc = new TCanvas(p_ps->GetName(),
			       p_ps->GetTitle());
    psc->SetLogx(logx);
    psc->SetLogy(logy);
    psc->SetLogz(logz);
    p_ps->Draw(option);
#else
    TVirtualPad *psc=gPad;
    Int_t logx=psc->GetLogx();
    Int_t logy=psc->GetLogy();
    Int_t logz=psc->GetLogz();
    psc->Divide(2,1);
    psc->cd(1);
    gPad->SetLogx(logx);
    gPad->SetLogy(logy);
    gPad->SetLogz(logz);
    TH2D::Draw(option);
    psc->cd(2);
    gPad->SetLogx(logx);
    gPad->SetLogy(logy);
    gPad->SetLogz(logz);
    p_ps->Draw(option);
#endif
  }
};// end of class PS_Histogram
#endif

class Camel: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point)
  {
    const double dx1=0.25, dy1=0.25, w1=1./0.004;
    const double dx2=0.75, dy2=0.75, w2=1./0.004;
    double weight=exp(-w1*((point[0]-dx1)*(point[0]-dx1)+
			   (point[1]-dy1)*(point[1]-dy1)))+
      exp(-w2*((point[0]-dx2)*(point[0]-dx2)+
	       (point[1]-dy2)*(point[1]-dy2)));
#ifdef ROOT_SUPPORT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("camel","camel",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"camel");
    }
    ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Camel

class Line: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point)
  {
    const double w1=1./0.004, mm=0.01;
    double weight=exp(-w1*sqr(point[1]+point[0]-1.0));
    if (point[0]<mm || point[0]>1.0-mm ||
	point[1]<mm || point[1]>1.0-mm) weight=0.0;
#ifdef ROOT_SUPPORT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("line","line",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"line");
    }
    ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Line

class Circle: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point)
  {
    const double dx1=0.4, dy1=0.6, rr=0.25, w1=1./0.004, ee=3.0;
    double weight=pow(point[1],ee)*
      exp(-w1*dabs(sqr(point[1]-dy1)+
		   sqr(point[0]-dx1)-sqr(rr)));
    weight+=pow(1.0-point[1],ee)*
      exp(-w1*dabs(sqr(point[1]-1.0+dy1)+
		   sqr(point[0]-1.0+dx1)-sqr(rr)));
#ifdef ROOT_SUPPORT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("circle","circle",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"circle");
    }
    ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Circle

class Fatal: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point)
  {
    const double cx=0.5, cy=0.5, mm=0.01;
    double weight=1.0;
    if (point[0]<mm || point[0]>1.0-mm ||
	point[1]<mm || point[1]>1.0-mm) weight=0.0;
    if ((point[0]<cx && point[1]<cy) ||
	(point[0]>cx && point[1]>cy)) weight=0.0;
#ifdef ROOT_SUPPORT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("fatal","fatal",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"fatal");
    }
    ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Fatal

class Chess: public Primitive_Integrand {
public:
  double operator()(const std::vector<double> &point)
  {
    const double mm=0.01;
    double weight=1.0;
    if (point[0]<mm || point[0]>1.0-mm ||
	point[1]<mm || point[1]>1.0-mm) weight=0.0;
    size_t ix=0, iy=0;
    for (double i=0.0;i<1.0;i+=0.125) {
      if (point[0]<=i) break;
      ++ix;
    }
    for (double i=0.0;i<1.0;i+=0.125) {
      if (point[1]<=i) break;
      ++iy;
    }
    if ((ix+iy)%2==0) weight=0.0;
#ifdef ROOT_SUPPORT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("chess","chess",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"chess");
    }
    ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Chess

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
    msg.Init(2,"");
    msg.SetModifiable(true);
    PRINT_INFO("Initialize integrator");
    Primitive_Integrator integrator;
    integrator.SetShuffleMode(1);
    integrator.SetDimension(2);
    integrator.SetNCells(4000);
    integrator.SetNOpt(1000);
    integrator.SetNMax(4000000);
    integrator.SetError(5.0e-4);
    // set default mode for
    // variance optimization
    integrator.SetMode(0);
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
    // these distributions cannot be integrated 
    // using variance optimization 
    // -> set maximum minimization
    integrator.SetMode(1);
    integrator.Initialize();
    PRINT_INFO("Integrate fatal");
    Fatal fatal;
    integrator.Integrate(&fatal);
    integrator.Initialize();
    PRINT_INFO("Integrate chess");
    Chess chess;
    integrator.Integrate(&chess);
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
