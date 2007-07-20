#include "Primitive_Integrator.H"
#include "Exception.H"
#include "My_Root.H"
#include <termios.h>
#include <unistd.h>
#ifdef USING__ROOT
#include "TH2D.h"
#endif

using namespace ATOOLS;

#ifdef USING__ROOT
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
  }
};// end of class PS_Histogram
#endif

class Camel: public Primitive_Integrand {
public:
  bool m_fill;
  Camel(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double dx1=0.25, dy1=0.25, w1=1./0.004;
    const double dx2=0.75, dy2=0.75, w2=1./0.004;
    double weight=exp(-w1*((point[0]-dx1)*(point[0]-dx1)+
			   (point[1]-dy1)*(point[1]-dy1)))+
      exp(-w2*((point[0]-dx2)*(point[0]-dx2)+
	       (point[1]-dy2)*(point[1]-dy2)));
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("camel","camel",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"camel");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Camel

class Line: public Primitive_Integrand {
public:
  bool m_fill;
  Line(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double w1=1./0.004, mm=0.01;
    double weight=exp(-w1*sqr(point[1]+point[0]-1.0));
    if (point[0]<mm || point[0]>1.0-mm ||
   	point[1]<mm || point[1]>1.0-mm) weight=0.0;
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("line","line",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"line");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Line

class Cross: public Primitive_Integrand {
public:
  bool m_fill;
  Cross(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double w1=1./0.004, mm=0.01;
    double weight(exp(-w1*Min(sqr(point[1]-point[0]),
			      sqr(point[1]+point[0]-1.0))));
    if (point[0]<mm || point[0]>1.0-mm ||
 	point[1]<mm || point[1]>1.0-mm) weight=0.0;
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("cross","cross",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"cross");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Cross

class Circle: public Primitive_Integrand {
public:
  bool m_fill;
  Circle(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double dx1=0.4, dy1=0.6, rr=0.25, w1=1./0.004, ee=3.0;
    double weight=exp(-w1*dabs(sqr(point[1]-dy1)+
			       sqr(point[0]-dx1)-sqr(rr)));
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("circle","circle",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"circle");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Circle

class Circles: public Primitive_Integrand {
public:
  bool m_fill;
  Circles(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double dx1=0.4, dy1=0.6, rr=0.25, w1=1./0.004, ee=1.0;
    double weight=Max(pow(point[1],ee)*
		      exp(-w1*dabs(sqr(point[1]-dy1)+
				   sqr(point[0]-dx1)-sqr(rr))),
		      pow(1.0-point[1],ee)*
		      exp(-w1*dabs(sqr(point[1]-1.0+dy1)+
				   sqr(point[0]-1.0+dx1)-sqr(rr))));
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("circles","circles",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"circles");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Circles

class Fatal: public Primitive_Integrand {
public:
  bool m_fill;
  Fatal(): m_fill(false) {}
  double operator()(const std::vector<double> &point)
  {
    const double cx=0.5, cy=0.5, mm=0.01;
    double weight=1.0;
    if (point[0]<mm || point[0]>1.0-mm ||
 	point[1]<mm || point[1]>1.0-mm) weight=0.0;
    if ((point[0]<cx && point[1]<cy) ||
	(point[0]>cx && point[1]>cy)) weight=0.0;
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("fatal","fatal",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"fatal");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Fatal

class Chess: public Primitive_Integrand {
public:
  bool m_fill;
  Chess(): m_fill(false) {}
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
#ifdef USING__ROOT
    static PS_Histogram *ps=NULL;
    if (ps==NULL) {
      ps = new PS_Histogram("chess","chess",100,0.,1.,100,0.,1.);
      MYROOT::myroot->AddObject(ps,"chess");
    }
    if (m_fill) ps->Fill(point[0],point[1],weight);
#endif
    return weight;
  }
};// end of class Chess

int main(int argc,char **argv)
{
#ifdef USING__ROOT
  MYROOT::myroot = new MYROOT::My_Root(argc,argv);
  exh->AddTerminatorObject(MYROOT::myroot);
#endif
  std::set_terminate(exh->Terminate);
  std::set_unexpected(exh->Terminate);
  signal(SIGSEGV,exh->SignalHandler);
  signal(SIGINT,exh->SignalHandler);
  signal(SIGBUS,exh->SignalHandler);
  signal(SIGFPE,exh->SignalHandler);
  signal(SIGABRT,exh->SignalHandler);
  signal(SIGTERM,exh->SignalHandler);
  signal(SIGXCPU,exh->SignalHandler);
  try {
    msg_Init(2,"");
    termios testos;
    if (tcgetattr(STDOUT_FILENO,&testos)==0) msg_SetModifiable(true);
    else msg_SetModifiable(false);
    PRINT_INFO("Initialize integrator");
    Primitive_Integrator integrator;
    integrator.SetShuffleMode(1);
    integrator.SetDimension(2);
    integrator.SetNCells(1000);
    integrator.SetNOpt(1000);
    integrator.SetNMax(2000000);
    integrator.SetError(5.0e-4);
    // set default mode for
    // variance optimization
    integrator.SetMode(imc::varopt);

#define CAMEL
#define LINE
#define CIRCLES
#define CHESS

#ifdef CAMEL
    integrator.Initialize();
    PRINT_INFO("Integrate camel");
    Camel camel;
    integrator.Integrate(&camel);
    camel.m_fill=true;
    for (size_t i(0);i<1000000;++i) integrator.Point();
#endif
#ifdef CROSS
    integrator.Initialize();
    PRINT_INFO("Integrate cross");
    Cross cross;
    integrator.Integrate(&cross);
    cross.m_fill=true;
    for (size_t i(0);i<1000000;++i) integrator.Point();
#endif
#ifdef LINE
    integrator.Initialize();
    PRINT_INFO("Integrate line");
    Line line;
    integrator.Integrate(&line);
    line.m_fill=true;
    for (size_t i(0);i<1000000;++i) integrator.Point();
#endif
#ifdef SCIRCLE
    integrator.Initialize();
    PRINT_INFO("Integrate circle");
    Circle circle;
    integrator.Integrate(&circle);
    circle.m_fill=true;
    for (size_t i(0);i<100000;++i) integrator.Point();
#endif
#ifdef CIRCLES
    integrator.Initialize();
    PRINT_INFO("Integrate circles");
    Circles circles;
    integrator.Integrate(&circles);
    circles.m_fill=true;
    for (size_t i(0);i<100000;++i) integrator.Point();
#endif
    // these distributions cannot be integrated 
    // using variance optimization 
    // -> set maximum minimization
    integrator.SetMode(imc::maxopt);
#ifdef FATAL
    integrator.Initialize();
    PRINT_INFO("Integrate fatal");
    Fatal fatal;
    integrator.Integrate(&fatal);
    fatal.m_fill=true;
    for (size_t i(0);i<1000000;++i) integrator.Point();
#endif
#ifdef CHESS
    integrator.Initialize();
    PRINT_INFO("Integrate chess");
    Chess chess;
    integrator.Integrate(&chess);
    chess.m_fill=true;
    for (size_t i(0);i<1000000;++i) integrator.Point();
#endif
#ifdef USING__ROOT
    MYROOT::myroot->SetDrawOption("lego2");
    delete MYROOT::myroot;
#endif
    return 0;
  }
  catch (Exception exception) {
    exception.UpdateLogFile();
    msg_Error()<<exception<<std::endl;
    std::terminate();
  }
  catch (std::exception exception) {
    std::cout<<"Sherpa: throws std::exception "<<exception.what()<<" ..."<<std::endl;
    std::terminate();
  }
}
