#include "Doubly_Unintegrated_PDF.H"

#include "Primitive_Integrator.H"
#include "NLL_Sudakov_Base.H"
#include "MyStrStream.H"
#include "MathTools.H"
#include "Random.H"
#include "Exception.H"
#include "Run_Parameter.H"

#ifdef PROFILE__All
#include "prof.hh"
#else
#ifdef PROFILE__DUPDF
#include "prof.hh"
#else
#define PROFILE_HERE
#endif
#endif

#if defined TEST__Sudakov || defined TEST__Weights || defined TEST__DPDF || defined TEST__IDPDF
#ifndef USING__ROOT 
#error ROOT is disabled
#endif
#include "My_Root.H"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#endif

using namespace PDF;

static size_t m_gridpoints=25;

class DGLAP_PDF: public ATOOLS::Primitive_Integrand {
private:

  PDF_Base *p_pdf;
  double    m_q0, m_x, m_z, m_kp2, m_kp;

  ATOOLS::Flavour m_fl;

public:

  DGLAP_PDF(PDF_Base *const pdf,const double &mu02); 

  double Integrate(const double &x,const double &z,const double &kp2);

  double operator()(const std::vector<double> &point);

};// end of class DGLAP_PDF

DGLAP_PDF::DGLAP_PDF(PDF_Base *const pdf,const double &mu02): 
  p_pdf(pdf), m_q0(sqrt(mu02)), m_fl(ATOOLS::kf::gluon) {}

double DGLAP_PDF::operator()(const std::vector<double> &point)
{
  double q=pow(10.0,point[0]), phi=point[1];
  double kpp2=ATOOLS::sqr((1.0-m_z)*q)+m_kp2+2.0*m_kp*(1.0-m_z)*q*cos(phi);
  if (kpp2<p_pdf->Q2Min()) return 0.0;
  // d x*pdf(x,q^2) / d log(q^2)
  p_pdf->Calculate(m_x/m_z,0.0,0.0,kpp2*(1.0+1.0e-9));
  double cur=p_pdf->GetXPDF(m_fl);
  p_pdf->Calculate(m_x/m_z,0.0,0.0,kpp2);
  cur=(cur-p_pdf->GetXPDF(m_fl))/1.0e-9;
  // integrand
  cur*=m_kp2/kpp2;
  msg_Debugging()<<"DGLAP_PDF::Value(): kp = "<<m_kp<<" GeV, q = "<<q
		 <<", phi = "<<(phi/M_PI*180.0)<<" ° -> "
		 <<cur<<"     "<<ATOOLS::bm::cr<<std::flush;
#ifdef TEST__IDPDF
  static TH2D *pint, *wint;
  static bool initi=false;
  if (!initi) {
    initi=true;
    pint = new TH2D("pint","pint",50,0.0,2.0*M_PI,100,-1.0,4.0);
    wint = new TH2D("wint","wint",50,0.0,2.0*M_PI,100,-1.0,4.0);
    MYROOT::myroot->AddObject(pint,"wint");
    MYROOT::myroot->AddObject(wint,"rint");
  }
  pint->Fill(phi,log10(q),1.0);
  wint->Fill(phi,log10(q),cur);
#endif
  return cur;
}

double DGLAP_PDF::Integrate(const double &x,const double &z,const double &kp2)
{
  msg_Tracking()<<ATOOLS::bm::cr<<"DGLAP_PDF::Integrate("
		<<x<<","<<z<<","<<kp2<<"):";
  m_x=x;
  m_z=z;
  m_kp=sqrt(m_kp2=kp2);
  ATOOLS::Primitive_Integrator integrator;
  std::vector<double> bound(2);
  bound[0]=log10(m_q0);
  bound[1]=0.0;
  integrator.SetMin(bound);
  bound[0]=log10(ATOOLS::rpa.gen.Ecms()/2.0);
  bound[1]=2.0*M_PI;
  integrator.SetMax(bound);
  integrator.SetMode(ATOOLS::imc::varopt);
  integrator.SetShuffleMode(0);
  integrator.SetNOpt(1000);
  integrator.SetNCells(1000);
  integrator.SetNMax(100000);
  integrator.SetError(0.01);
  int old=ATOOLS::msg->Level();
  ATOOLS::msg->SetLevel(0);
  double cur=integrator.Integrate(this);
  ATOOLS::msg->SetLevel(old);
  msg_Tracking()<<" Result = "<<cur<<"      "<<ATOOLS::bm::cr<<std::flush;
  return cur;
}

Doubly_Unintegrated_PDF::
Doubly_Unintegrated_PDF(PDF_Base *_p_pdf,MODEL::Running_AlphaS *_p_alphas,
			const double mu02):
  p_pdf(_p_pdf),
  p_alphas(_p_alphas),
  m_mu02(mu02),
  m_epsilon(1.e-2),
  m_kperpscheme(kps::function),
  m_fixedktexponent(0.125),
  m_mode(0), m_sudmode(1), m_splitmode(0),
  p_integral(new ATOOLS::Grid(3,m_gridpoints))
{
  m_type=std::string("DUPDF(")+p_pdf->Type()+std::string(")");
  m_xmin=p_pdf->XMin();
  m_xmax=p_pdf->XMax();
  m_q2min=0.0;
  m_q2max=ATOOLS::sqr(p_pdf->Q2Max());
  m_bunch=p_pdf->Bunch();
  m_partons=p_pdf->Partons();
}

Doubly_Unintegrated_PDF::~Doubly_Unintegrated_PDF()
{
#ifdef TEST__Weights
  static TH2D *p_weights;
  static bool initw=false;
  if (!initw) {
    initw=true;
    p_weights = new TH2D("DUPDF_Weights","DUPDF_Weights",
			 100,-2.0,3.0,100,-2.0,3.0);
    MYROOT::myroot->AddObject(p_weights,"DUPDF_Weights");
  }
  for (double i=-2.0;i<3.0;i+=0.05) {
    for (double j=-2.0;j<3.0;j+=0.05) {
      Calculate(0.01,0.1,exp(i),exp(j));
      std::cout<<i<<" "<<j<<" "<<exp(i)<<" "<<exp(j)<<std::endl;
      p_weights->Fill(i,j,GetXPDF(ATOOLS::kf::gluon));
    }
  }
  THROW(normal_exit,"finished histogram");
#endif
#ifdef TEST__Sudakov
  static bool inits=false;
  if (!inits) {
    inits=true;
    double min=0.082, max=91.0, step=0.1;
    TH1D *sudu = new TH1D(ATOOLS::ToString((long int)this).c_str(),
			  "Quark_Sudakov",(int)((log(max)-log(min))/step),
			  log(min/max),0.);
    TH1D *sudg = new TH1D(ATOOLS::ToString(1+(long int)this).c_str(),
			  "Gluon_Sudakov",(int)((log(max)-log(min))/step),
			  log(min/max),0.);
    TH2D *sudu2 = new TH2D(ATOOLS::ToString(2+(long int)this).c_str(),
			   "Quark_Sudakov_2D",
			   (int)((log(max)-log(min))/step),
			   log(min/max),0.,(int)((log(max)-log(min))/step),
			   log(min/max),0.);
    TH2D *sudg2 = new TH2D(ATOOLS::ToString(3+(long int)this).c_str(),
			   "Gluon_Sudakov_2D",
			   (int)((log(max)-log(min))/step),
			   log(min/max),0.,(int)((log(max)-log(min))/step),
			   log(min/max),0.);
    MYROOT::myroot->AddObject(sudu,"sudu");
    MYROOT::myroot->AddObject(sudg,"sudg");
    MYROOT::myroot->AddObject(sudu2,"sudu2");
    MYROOT::myroot->AddObject(sudg2,"sudg2");
    for (double logkt=log(min)+step;logkt<log(max);logkt+=step) {
      double kt=exp(logkt);
      sudu->Fill(log(kt/91.),p_sudakov->Delta(ATOOLS::kf::u)(91.,kt));
      sudg->Fill(log(kt/91.),p_sudakov->Delta(ATOOLS::kf::gluon)(91.,kt));
    }
    for (double logmu=log(min)+step;logmu<log(max);logmu+=step) {
      for (double logkt=log(min)+step;logkt<log(max);logkt+=step) {
	double kt=exp(logkt), mu=exp(logmu);
	if (logkt<logmu) {
	  sudu2->Fill(log(kt/91.),log(mu/91.),
		      p_sudakov->Delta(ATOOLS::kf::u)(mu,kt));
	  sudg2->Fill(log(kt/91.),log(mu/91.),
		      p_sudakov->Delta(ATOOLS::kf::gluon)(mu,kt));
	}
      }
    }
    double linstep=.002;
    TH2D *jetrate = new TH2D(ATOOLS::ToString(4+(long int)this).c_str(),
			     "Jet_Rate",
			     (int)(1./linstep),0.,1.,
			     (int)(1./linstep),0.,1.);
    MYROOT::myroot->AddObject(jetrate,"jetrate");
    for (double x1=linstep;x1<1.;x1+=linstep) {
      for (double x2=linstep;x2<1.;x2+=linstep) {
	double x3=2.-x1-x2;
	if (x3>linstep && x1>linstep && x2>linstep &&
	    x3<1.-linstep) {
	  jetrate->Fill(x1,x2,1./(1.-x1)/(1.-x2)*
			((1.-x1)/x3*(1.+pow((x1/(2.-x2)),2)) 
			 +(1.-x2)/x3*(1.+pow((x2/(2.-x1)),2))));
	}
      }
    }
    THROW(normal_exit,"Finished Sudakov histograms");
  }
#endif
  while (m_branching.size()>0) {
    delete (*m_branching.begin()).second;
    m_branching.erase(m_branching.begin());
  }
  delete p_sudakov;
  delete p_integral;
}

bool Doubly_Unintegrated_PDF::Initialize()
{
  for (size_t i=0;i<m_partons.size();++i) {
    if (m_partons[i].Size()>1) {
      for (int j=0;j<m_partons[i].Size();++j) {
	m_branching[m_partons[i][j]] = 
	  new LL_Branching(m_partons[i][j],p_alphas);
      }
    }
    else {
      m_branching[m_partons[i]] = new LL_Branching(m_partons[i],p_alphas);
    }
  }
  p_sudakov = new LL_Sudakov(p_alphas);
  p_sudakov->SetOutPath(std::string("DUPDF_")+
			ATOOLS::ToString(ATOOLS::rpa.gen.Ecms()));
  p_sudakov->Initialize();
  if (m_mode>0) return InitializeGluonPDF();
  return true;
}

bool Doubly_Unintegrated_PDF::InitializeGluonPDF()
{
  const std::string gridname="GDUPDF_"+p_pdf->Type()+"_"+
    ATOOLS::ToString(m_mu02)+"_"+ATOOLS::ToString(ATOOLS::rpa.gen.Ecms());
  if (!p_integral->ReadIn(gridname)) {
    msg_Info()<<"Doubly_Unintegrated_PDF::Doubly_Unintegrated_PDF(..):"
	      <<"Initializing DUPDF grid ... "<<std::flush;
    DGLAP_PDF pdf(p_pdf,m_mu02);
    std::vector<size_t> i(3);
    std::vector<double> xz(m_gridpoints), q2(m_gridpoints);
    double dxz=6.0/(m_gridpoints+1);
    double dq=4.0/m_gridpoints;
    for (i[0]=0;i[0]<m_gridpoints;++i[0]) {
      xz[i[0]]=-6.0+i[0]*dxz;
      for (i[1]=i[0]+1;i[1]<m_gridpoints;++i[1]) {
	double z=-6.0+i[1]*dxz;
	for (i[2]=0;i[2]<m_gridpoints;++i[2]) {
	  q2[i[2]]=0.1+i[2]*dq;
	  double cur=pdf.Integrate(pow(10.0,xz[i[0]]),pow(10.0,z),
				   pow(10.0,2.0*q2[i[2]]));
	  p_integral->SetY(i,cur);
#ifdef TEST__DPDF
	  static TH3D *dpdfh;
	  static bool ini=false;
	  if (!ini) {
	    dpdfh = new TH3D("GUPDF","GUPDF",
			     m_gridpoints,-6.0,-6.0+dxz*(m_gridpoints-1),
			     m_gridpoints,-6.0,-6.0+dxz*(m_gridpoints-1),
			     m_gridpoints,0.1,0.1+dq*m_gridpoints);
	    MYROOT::myroot->AddObject(dpdfh,"GUPDF");
	    ini=true;
	  }
	  dpdfh->Fill(xz[i[0]],z,q2[i[2]],cur);
#endif
	}
      }
    }
    p_integral->SetX(0,xz);
    p_integral->SetX(1,xz);
    p_integral->SetX(2,q2);
    p_integral->WriteOut(gridname);
    msg_Info()<<"done."<<std::endl;
  }
  return true;
}

void Doubly_Unintegrated_PDF::AssignKeys(ATOOLS::Integration_Info *const info)
{
  p_sudakov->AssignKeys(info);
}

double Doubly_Unintegrated_PDF::ConstantIntegrated(ATOOLS::Flavour flavour)
{
  p_pdf->Calculate(m_x,0.,0.,m_mu02);                                         
  double integrated=p_pdf->GetXPDF(flavour);
  integrated*=p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_mu02)); 
  return integrated/((1.-m_x)*m_mu02);
}

double Doubly_Unintegrated_PDF::SmoothIntegrated(ATOOLS::Flavour flavour)
{
  switch (m_kperpscheme) {
  case kps::function: {
    double savekp2=m_kperp2;
    this->Calculate(m_x,m_z,m_mu02*(1.+m_epsilon),m_mu2);
    double fprime=this->GetXPDF(flavour)*m_mu02;
    this->Calculate(m_x,m_z,m_mu02,m_mu2);
    double f=this->GetXPDF(flavour)*m_mu02;
    fprime-=f;
    fprime/=m_mu02*m_epsilon;
    m_unintegrated=0.;
    m_kperp2=savekp2;
    p_pdf->Calculate(m_x,0.,0.,m_mu02);
    m_integrated=p_pdf->GetXPDF(flavour);
    m_integrated*=p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_mu02));
    m_integrated/=(1.-m_x);
// #define USING__Variable_Exponent
#ifdef USING__Variable_Exponent
    double m1=1.+m_fixedktexponent, m2=m1*2.;
    double b=(m1+1)*(f-m1*m_integrated);
    double a=(f-b);
    double d1=(m1*(a+b)+b)/m_mu02, d2=d1;
    double accuracy=ATOOLS::Max(ATOOLS::rpa.gen.Accu(),1.0e-6*fprime);
    int step=0;
    do {
      b=(m2+1)*(f-m2*m_integrated);
      a=(f-b);
      d2=(m2*(a+b)+b)/m_mu02;
      double d=(d2-d1)/(m2-m1);
      m1=m2;
      d1=d2;
      m2=m2-(d2-fprime)/d;
      if (++step>100 || (!(d2>=0.) && !(d2<0.))) {
	msg_Error()<<"Doubly_Unintegrated_PDF::SmoothIntegrated(..): "
			   <<"Integrate failed.\n   x = "<<m_x<<", z = "<<m_z
			   <<" k_\\perp^2 = "<<m_kperp2<<" \\mu^2 = "<<m_mu2
			   <<"( "<<step<<" steps m = "<<m2<<" )"<<std::endl;
	return 0.;
      }
    } while (ATOOLS::dabs(d2-fprime)>accuracy);
    if (m2<=ATOOLS::rpa.gen.Accu()) {
      msg_Error()<<"Doubly_Unintegrated_PDF::SmoothIntegrated(..): "
			 <<"Integrate failed.\n   x = "<<m_x<<", z = "<<m_z
			 <<" k_\\perp^2 = "<<m_kperp2<<" \\mu^2 = "<<m_mu2
			 <<"( "<<step<<" steps m = "<<m2<<" )"<<std::endl;
      return 0.;
    }
    return (a+b*m_kperp2/m_mu02)*pow(m_kperp2/m_mu02,m2-1)/m_mu02;
#else
    double m=1.+m_fixedktexponent;
    double b=(m+1)*(f-m*m_integrated);
    double a=(f-b);
    return (a+b*m_kperp2/m_mu02)*pow(m_kperp2/m_mu02,m-1)/m_mu02;
#endif
  }
  case kps::derivative: {
    double savekp2=m_kperp2;
    this->Calculate(m_x,m_z,m_mu02*(1.+m_epsilon),m_mu2);
    double fprime=this->GetXPDF(flavour)*m_mu02;
    this->Calculate(m_x,m_z,m_mu02,m_mu2);
    double f=this->GetXPDF(flavour)*m_mu02;
    fprime-=f;
    fprime/=m_mu02*m_epsilon;
    m_unintegrated=0.;
    m_kperp2=savekp2;
    p_pdf->Calculate(m_x,0.,0.,m_mu02);
    m_integrated=p_pdf->GetXPDF(flavour);
    m_integrated*=p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_mu02));
    m_integrated/=(1.-m_x);
    double m=1.+m_fixedktexponent;
    double Ft=m_integrated/pow(m_mu02,m), ft=f/pow(m_mu02,m);
    double fpt=fprime/pow(m_mu02,m-1.);
    double c=(m+2.)*(m*(m+1.)*Ft-(2.*m+1.)*ft+fpt);
    double b=fpt-m*ft-c;
    c/=2.;
    double a=ft-b-c;
    c/=m_mu02*m_mu02;
    b/=m_mu02;
    return (a+b*m_kperp2+c*m_kperp2*m_kperp2)*pow(m_kperp2,m-1.);
  }
  default: return 0.;
  }
  return 0.;
}
		 
bool Doubly_Unintegrated_PDF::Unintegrate(ATOOLS::Flavour flavour)
{
  PROFILE_HERE;
  m_unintegrated=m_integrated=0.;
  if (Collinear(m_kperp2)) {
    switch (m_kperpscheme) {
    case kps::function: 
    case kps::derivative:
      m_integrated=SmoothIntegrated(flavour); 
      return true;
    case kps::constant: 
      m_integrated=ConstantIntegrated(flavour); 
      return true;
    }
    return false;     
  }
  switch (m_mode) {
  case 1: {
    LL_Branching::SF_Set::iterator sfit=
      LL_Branching::AllSplittings().begin();
    for (;sfit!=LL_Branching::AllSplittings().end();++sfit) {
      if ((*sfit)->GetFlB()==flavour) {
	if ((*sfit)->GetFlC().IsGluon() && 
	    m_z*(1.+sqrt(m_kperp2/m_mu2))>1.) continue;
	double cur=(*(*sfit))(m_z);
	if (flavour.IsGluon() && (*sfit)->GetFlA().IsGluon()) {
 	  cur-=6.0/m_z;
	  cur*=2.0;
	}
	m_unintegrated+=cur*p_pdf->GetXPDF((*sfit)->GetFlA());
      }
    }
    if (flavour==ATOOLS::kf::gluon && 
	m_z*(1.+sqrt(m_kperp2/m_mu2))<=1.) {
      std::vector<double> x(3,0.5*log10(m_kperp2));
      x[0]=log10(m_x);
      x[1]=log10(m_z);
      p_pdf->Calculate(m_x/m_z,0.0,0.0,m_kperp2*(1.0+1.0e-9));
      double cur=p_pdf->GetXPDF(flavour);
      p_pdf->Calculate(m_x/m_z,0.0,0.0,m_kperp2);
      cur=(cur-p_pdf->GetXPDF(flavour))/1.0e-9;
      m_unintegrated+=2.0*6.0/m_z*((*p_integral)(x)/M_PI-
				   cur*(log(m_kperp2)-log(m_mu02)));
    }
    break;
  }
  default: {
    LL_Branching::SF_Set::iterator sfit=
      LL_Branching::AllSplittings().begin();
    for (;sfit!=LL_Branching::AllSplittings().end();++sfit) {
      if ((*sfit)->GetFlB()==flavour) {
  	if ((*sfit)->GetFlC().IsGluon() && 
  	    m_z*(1.+sqrt(m_kperp2/m_mu2))>1.) continue;
  	if (m_splitmode==1 &&
	    !(*sfit)->GetFlA().IsGluon()) continue;
	m_unintegrated+=(*(*sfit))(m_z)*p_pdf->GetXPDF((*sfit)->GetFlA());
	if (flavour.IsGluon() && (*sfit)->GetFlA().IsGluon())
	  m_unintegrated+=(*(*sfit))(m_z)*p_pdf->GetXPDF((*sfit)->GetFlA());
      }
    }
    break;
  }
  }
  m_unintegrated*=(*p_alphas)(m_kperp2)/(2.0*M_PI);
  if (m_sudmode>0) 
    m_unintegrated*=p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_kperp2));
#ifdef TEST__Doubly_Unintegrated_PDF
  std::cout<<flavour<<" sud / pdf : "<<p_sudakov->
    Delta(flavour)(sqrt(m_mu2),sqrt(m_kperp2))
   	   <<" / "<<m_unintegrated<<std::endl;
#endif
  m_unintegrated/=m_kperp2;
  return true;
}

void Doubly_Unintegrated_PDF::Calculate(double x,double z,
					double kperp2,double mu2)
{
  PROFILE_HERE;
  m_x=x;
  m_z=z;
  m_kperp2=kperp2;
  m_mu2=mu2;
  m_calculate=true;
  if (m_z<m_x || m_kperp2>p_pdf->Q2Max()) {
    m_calculate=false;
    return; 
  }
  if (m_kperp2<p_pdf->Q2Min()) return;
  p_pdf->Calculate(m_x/m_z,0.,0.,m_kperp2);
}

double Doubly_Unintegrated_PDF::GetXPDF(const ATOOLS::Flavour flavour)
{
  if (!m_calculate) return 0.;
  if (!Unintegrate(flavour)) {
    msg_Error()<<"Doubly_Unintegrated_PDF::GetXPDF("<<flavour<<"): "
		       <<"Cannot unintegrate PDF."<<std::endl;
  }
  return m_integrated+m_unintegrated;
}

void Doubly_Unintegrated_PDF::Output()
{ 
  return p_pdf->Output(); 
}

PDF_Base *Doubly_Unintegrated_PDF::GetCopy()
{ 
  return p_pdf->GetCopy(); 
}

bool Doubly_Unintegrated_PDF::Collinear(const double kp2) const
{
  return kp2<m_mu02;
}

PDF_Base *Doubly_Unintegrated_PDF::GetBasicPDF() 
{
  return p_pdf;
}

double Doubly_Unintegrated_PDF::Cut(const std::string &type) 
{
  if (type=="kp") return m_mu02;
  return PDF_Base::Cut(type);
}
