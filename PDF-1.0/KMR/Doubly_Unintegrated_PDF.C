#include "Doubly_Unintegrated_PDF.H"

#include "NLL_Sudakov_Base.H"
#include "MyStrStream.H"

#ifdef PROFILE__All
#include "prof.hh"
#else
#ifdef PROFILE__DUPDF
#include "prof.hh"
#else
#define PROFILE_HERE
#endif
#endif

#include "Exception.H"
#ifdef TEST__Sudakov
#ifdef ROOT_SUPPORT
#include "Scaling.H"
#include "My_Root.H"
#include "TH1D.h"
#include "TH2D.h"
static bool initialized=false;
#else
#error ROOT is disabled
#endif
#endif

using namespace PDF;

Doubly_Unintegrated_PDF::Doubly_Unintegrated_PDF(PDF_Base *_p_pdf,MODEL::Running_AlphaS *_p_alphas,
						 const double mu02,const kps::type kperpscheme):
  p_pdf(_p_pdf),
  p_alphas(_p_alphas),
  m_mu02(mu02),
  m_epsilon(1.e-2),
  m_kperpscheme(kperpscheme)
{
  m_type=std::string("DUPDF(")+p_pdf->Type()+std::string(")");
  m_xmin=m_xmax=0.0;
  m_xmax=p_pdf->XMax();
  m_q2min=0.0;
  m_q2max=ATOOLS::sqr(p_pdf->Q2Max());
  m_bunch=p_pdf->Bunch();
  m_partons=p_pdf->Partons();
  for (size_t i=0;i<m_partons.size();++i) {
    if (m_partons[i].Size()>1) {
      for (int j=0;j<m_partons[i].Size();++j) {
	m_branching[m_partons[i][j]] = new LL_Branching(m_partons[i][j],p_alphas);
      }
    }
    else {
      m_branching[m_partons[i]] = new LL_Branching(m_partons[i],p_alphas);
    }
  }
  p_sudakov = new LL_Sudakov(p_alphas);
  std::string ecms;
  MyStrStream sstr;
  sstr<<ATOOLS::rpa.gen.Ecms();
  sstr>>ecms;
  p_sudakov->SetOutPath(std::string("DUPDF_")+ecms);
  p_sudakov->Initialize();
}

Doubly_Unintegrated_PDF::~Doubly_Unintegrated_PDF()
{
#ifdef TEST__Sudakov
#ifdef ROOT_SUPPORT
  if (!initialized) {
    initialized=true;
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
    throw(ATOOLS::Exception(ATOOLS::ex::normal_exit,
			    "Finished Sudakov histograms","DUPDF","DUPDF"));
  }
#else
#error ROOT is disabled
#endif
#endif
  while (m_branching.size()>0) {
    delete (*m_branching.begin()).second;
    m_branching.erase(m_branching.begin());
  }
  delete p_sudakov;
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
  double c=3.*(fprime+(2.*m_integrated-3.*f)/m_mu02);
  double b=fprime-f/m_mu02-c;
  c/=2.;
  double a=f/m_mu02-b-c;
  c/=m_mu02*m_mu02;
  b/=m_mu02;
  return a+b*m_kperp2+c*m_kperp2*m_kperp2;
}

bool Doubly_Unintegrated_PDF::Unintegrate(ATOOLS::Flavour flavour)
{
  PROFILE_HERE;
  m_unintegrated=m_integrated=0.;
  if (m_kperp2<m_mu02) {
    switch (m_kperpscheme) {
    case kps::smooth: m_integrated=SmoothIntegrated(flavour); 
      return true;
    case kps::constant: m_integrated=ConstantIntegrated(flavour); 
      return true;
    }
    return false;     
  }
  LL_Branching::SF_Set::iterator sfit=LL_Branching::AllSplittings().begin();
  for (;sfit!=LL_Branching::AllSplittings().end();++sfit) {
    if ((*sfit)->GetFlB()==flavour) {
      if (((flavour.IsGluon() && (*sfit)->GetFlA().IsGluon()) ||
	   (flavour.IsQuark() && (*sfit)->GetFlA().IsQuark())) && 
	  m_z*(1.+sqrt(m_kperp2/m_mu2))>1.) continue;
      m_unintegrated+=(*(*sfit))(m_z)*p_pdf->GetXPDF((*sfit)->GetFlA());
      if (flavour.IsGluon() && (*sfit)->GetFlA().IsGluon())
	m_unintegrated+=(*(*sfit))(m_z)*p_pdf->GetXPDF((*sfit)->GetFlA());
    }
  }
  m_unintegrated*=(*p_alphas)(m_kperp2)/(2.0*M_PI);
  m_unintegrated*=p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_kperp2));
#ifdef TEST__Doubly_Unintegrated_PDF
  std::cout<<flavour<<" sud / pdf : "<<p_sudakov->Delta(flavour)(sqrt(m_mu2),sqrt(m_kperp2))
   	   <<" / "<<m_unintegrated<<std::endl;
#endif
  m_unintegrated/=m_kperp2;
  return true;
}

void Doubly_Unintegrated_PDF::Calculate(double x,double z,double kperp2,double mu2)
{
  PROFILE_HERE;
  m_x=x;
  m_z=z;
  m_kperp2=kperp2;
  m_mu2=mu2;
  m_calculate=true;
  if (m_z<m_x || m_kperp2>p_pdf->Q2Max() ||
      (m_mu2<p_pdf->Q2Min() && m_mu02<m_kperp2)) {
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
    ATOOLS::msg.Error()<<"Doubly_Unintegrated_PDF::GetXPDF("<<flavour<<"): "
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
