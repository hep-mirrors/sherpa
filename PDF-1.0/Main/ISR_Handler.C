#include "ISR_Handler.H"

#include "Intact.H"
#include "Structure_Function.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H" 
#include "Info_Key.H"
#include "Exception.H"
#include "Random.H"
#include "ISR_Info.H"
#include "Blob.H"
#include <stdio.h>

#define USING__No_Exact_Kinematics

#ifdef TEST__ISR_Handler
#include "Gauss_Integrator.H"
#endif

#ifdef NO_ANALYSIS__all
#define NO_ANALYSIS__ISR_Handler
#endif

#ifdef PROFILE__All
#include "prof.hh"
#else
#ifdef PROFILE__ISR_Handler
#include "prof.hh"
#else
#define PROFILE_HERE
#endif
#endif

using namespace ATOOLS;
using namespace PDF;
using namespace std;

#ifdef TEST__ISR_Handler
class TestPDF: public Function_Base {
public:
  ATOOLS::Flavour m_fl;
  PDF_Base *p_pdf;
  double m_x, m_z, m_kt2, m_mu2;
  TestPDF(): m_fl(ATOOLS::kf::u) {};
  double operator()(double z) {
    p_pdf->Calculate(m_x,z,m_kt2,m_mu2);
    return p_pdf->GetXPDF(m_fl);
  };
  double operator()() {
    p_pdf->Calculate(m_x,m_z,m_kt2,m_mu2);
    return p_pdf->GetXPDF(m_fl);
  };
};
#endif

double Lambda2(double sp,double sp1,double sp2) 
{ 
  return (sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
}

ISR_Handler::ISR_Handler(ISR_Base **isrbase,const double *splimits,const double *kplimits):
  p_isrbase(isrbase),
  p_info(new ATOOLS::Integration_Info()),
  m_kperpscheme((int)kps::constant),
  m_weight(1.),
  m_info_lab(8),
  m_info_cms(8)
{
  m_mode=0;
  m_kmrmode=0;
  for (short int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      p_isrbase[i]->AssignKeys(p_info);
      m_mode += i+1;
      if (p_isrbase[i]->KMR()) m_kmrmode += i+1;
    }
  }
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());
  m_x[1]=m_x[0]=1.; 
  Init(splimits,kplimits);
  Doubly_Unintegrated_PDF *dupdf=
    dynamic_cast<Doubly_Unintegrated_PDF*>(p_isrbase[0]->PDF());
  if (dupdf!=NULL) m_kperpscheme=dupdf->KPerpScheme();
#ifdef TEST__ISR_Handler
  TestPDF testpdf;
  testpdf.p_pdf=p_isrbase[0]->PDF();
  testpdf.m_mu2=8100.;
  Gauss_Integrator gauss(&testpdf);
  std::ofstream *out = new std::ofstream("pdftest");
  for (double test=-1.;test<=4.;test+=0.1) {
    testpdf.m_z=0.8;
    testpdf.m_kt2=exp(test);
    testpdf.m_x=0.1;
    *out<<exp(test)<<" "<<exp(test)*gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.01;
    *out<<" "<<exp(test)*gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.001;
    *out<<" "<<exp(test)*gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.0001;
    *out<<" "<<exp(test)*gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.1;
    *out<<exp(test)<<" "<<exp(test)*testpdf();
    testpdf.m_x=0.01;
    *out<<" "<<exp(test)*testpdf();
    testpdf.m_x=0.001;
    *out<<" "<<exp(test)*testpdf();
    testpdf.m_x=0.0001;
    *out<<" "<<exp(test)*testpdf()<<std::endl;
  }
  delete out;
  throw(ATOOLS::Exception(ATOOLS::ex::normal_exit,"finished integration"));
#endif
}

ISR_Handler::~ISR_Handler() 
{
  if (p_isrbase) {
    for (int i=0;i<2;i++) {
      if (p_isrbase[i]) delete p_isrbase[i];  
    }
    delete[] p_isrbase; p_isrbase = 0;
  }
  delete p_info;
}

void ISR_Handler::Init(const double *splimits,const double *kplimits) 
{
  m_type = p_isrbase[0]->Type()+std::string("*")+p_isrbase[1]->Type();
  double s = sqr(ATOOLS::rpa.gen.Ecms());
  m_splimits[0] = s*splimits[0];
  m_splimits[1] = ATOOLS::Min(s*splimits[1],s*Upper1()*Upper2());
  m_splimits[2] = s;
  m_fixed_smin = m_splimits[0];
  m_fixed_smax = m_splimits[1];
  m_ylimits[0] = -10.;
  m_ylimits[1] = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * p_isrbase[0]->Exponent() * p_isrbase[1]->Exponent();
  m_zlimits[0] = 0.;
  m_zlimits[1] = 1.;
  if (p_isrbase[0]->PDF()!=NULL) {
    m_kplimits[0] = kplimits[0]*s;
    m_kplimits[1] = p_isrbase[0]->Cut("kp");
    m_kplimits[2] = kplimits[1]*s;
  }
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());
  double E=ATOOLS::rpa.gen.Ecms();
  double x=1./2.+(m_mass2[0]-m_mass2[1])/(2.*E*E);
  double E1=x*E;
  double E2=E-E1;
  m_fixvecs[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass2[0]));
  m_fixvecs[1]=Vec4D(E2,0.,0.,-m_fixvecs[0][3]);
}

void ISR_Handler::SetSprimeMin(const double spmin)       
{ 
  m_splimits[0] = Max(m_fixed_smin,spmin); 
}

void ISR_Handler::SetSprimeMax(const double spmax)       
{ 
  m_splimits[1] = Min(m_fixed_smax,spmax); 
}

void ISR_Handler::SetFixedSprimeMin(const double spmin)  
{ 
  m_fixed_smin = Max(m_fixed_smin,spmin);
  m_splimits[0] = Max(m_splimits[0],spmin);
}

void ISR_Handler::SetFixedSprimeMax(const double spmax)  
{
  m_fixed_smax = Min(m_fixed_smax,spmax);
  m_splimits[1] = Min(m_splimits[1],spmax);
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *bunches,ATOOLS::Flavour *partons) 
{
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      if (bunches[i] != PDF(i)->Bunch()) { fit = 0; break; }
      fit = 0;
      for (unsigned int j = 0;j<(PDF(i)->Partons()).size();j++) {
	if (partons[i] == (PDF(i)->Partons())[j]) {
	  fit = 1;
	  break; 
	}
      }
      if (fit == 0) break;
    }
    else {
      if (partons[i]!=p_isrbase[i]->Flavour()) {
	fit = 0;
	break;
      }
    }
  }
  return fit;
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *partons) 
{
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      fit = 0;
      for (unsigned int j = 0;j<(PDF(i)->Partons()).size();j++) {
	if (partons[i] == (PDF(i)->Partons())[j]) {
	  fit = 1;
	  break; 
	}
      }
      if (fit == 0) break;
    }
    else {
      if (partons[i]!=p_isrbase[i]->Flavour()) {
	fit = 0;
	break;
      }
    }
  }
  return fit;
}

void ISR_Handler::SetPartonMasses(Flavour *fl) 
{
  m_mass2[0]=sqr(fl[0].Mass());
  m_mass2[1]=sqr(fl[1].Mass());
  double E=ATOOLS::rpa.gen.Ecms();
  double x=1./2.+(m_mass2[0]-m_mass2[1])/(2.*E*E);
  double E1=x*E;
  double E2=E-E1;
  m_fixvecs[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass2[0]));
  m_fixvecs[1]=Vec4D(E2,0.,0.,-m_fixvecs[0][3]);
}

bool ISR_Handler::MakeISR(Vec4D *const p,const size_t n,
			  const ATOOLS::Flavour *flavs,const size_t nflavs) 
{
  PROFILE_HERE;
  p_info->ResetAll();
  m_weight=1.;
  if (m_mode==0) {
    m_kpkey[1][3]=m_kpkey[0][3]=0.;
    m_x[1]=m_x[0]=1.;
    m_zkey[1][2]=m_zkey[0][2]=1.;
    m_flux=.25;
    m_flux/=sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
    return true;
  }
  if (m_spkey[3]<m_splimits[0] || m_spkey[3]>m_splimits[1]) {
    ATOOLS::msg.Error()<<"ISR_Handler::MakeISR(..): "<<std::endl
		       <<om::red<<" sprime out of bounds :"<<om::reset
		       <<" s'_{min}, s'_{max 1,2} vs. s': "<<m_splimits[0]
		       <<", "<<m_splimits[1]<<", "<<m_splimits[2]
		       <<" vs. "<<m_spkey[3]<<endl;
    return false;
  }
  double Q=sqrt(m_splimits[2]), E=sqrt(m_spkey[3]);
  double E1=E*(1./2.+(m_mass2[0]-m_mass2[1])/(2.*m_spkey[3]));
  p_cms[0]=p[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass2[0]));
  p_cms[1]=p[1]=Vec4D(E-E1,(-1.)*Vec3D(p[0]));
  Vec4D plab[2]; plab[0]=p[0]; plab[1]=p[1];
  m_cmsboost=Poincare(Vec4D(cosh(m_ykey[2]),0.,0.,sinh(m_ykey[2])));
  m_cmsboost.BoostBack(p_cms[0]);
  m_cmsboost.BoostBack(p_cms[1]);
  m_x[0]=2.*p_cms[0][0]/Q;
  m_x[1]=2.*p_cms[1][0]/Q;
  m_flux=.25;
  m_flux/=sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
  if (!m_kmrmode) {
#ifndef NO_ANALYSIS__ISR_Handler
    m_info_lab[iic::E_1]=m_x[0];
    m_info_lab[iic::E_2]=m_x[1];
#endif
    return true;
  }
  double phi=0.0;
  for (size_t i=0;i<2;++i) {
    phi+=2.0*M_PI*ran.Get();
    double kp=sqrt(m_kpkey[i][3]); 
#ifndef USING__No_Exact_Kinematics
    if (m_kperpscheme==(int)kps::constant &&
    	p_isrbase[i]->Collinear(m_kpkey[i][3])) m_zkey[i][2]=0.;
#endif
    m_kp[i]=Vec4D(0.0,cos(phi)*kp,sin(phi)*kp,0.0);
  }
  E=sqrt(m_spkey[3]-(m_kp[0]+m_kp[1]).Abs2());
  double xi=exp(m_ykey[2]);
  double D1=E/Q*xi;
  double D2=E/Q/xi;
  double C1=m_zkey[0][2]/(1.-m_zkey[0][2])*m_kpkey[0][3]/m_splimits[2];
  double C2=m_zkey[1][2]/(1.-m_zkey[1][2])*m_kpkey[1][3]/m_splimits[2];
  m_x[0]=.5*(D1+(C2-C1)/D2);
  m_x[0]=m_x[0]+sqrt(m_x[0]*m_x[0]+C1*D1/D2);
  m_x[1]=(D2*m_x[0]+C1-C2)/D1;
  double b1=C1/m_x[0], b2=C2/m_x[1];
  if (m_zkey[0][2]>1. || m_zkey[1][2]>1.) return false;
  if (m_x[0]>m_zkey[0][2] || m_x[1]>m_zkey[1][2]) return false;
  if (b1>m_x[1] || b2>m_x[0]) return false;
  if (b1<0. || b2<0.) return false;
  p[0]=Vec4D((m_x[0]-b1)*Q/2.,m_kp[0][1],m_kp[0][2],(m_x[0]+b1)*Q/2.);
  p[1]=Vec4D((m_x[1]-b2)*Q/2.,m_kp[1][1],m_kp[1][2],-(m_x[1]+b2)*Q/2.);
  double a1=m_zkey[0][2]==0.0?0.0:m_x[0]*(1./m_zkey[0][2]-1.);
  p_kmrlast[0]=Vec4D((a1+b1)*Q/2.,-m_kp[0][1],-m_kp[0][2],(a1-b1)*Q/2.);
  double a2=m_zkey[1][2]==0.0?0.0:m_x[1]*(1./m_zkey[1][2]-1.);
  p_kmrlast[1]=Vec4D((a2+b2)*Q/2.,-m_kp[1][1],-m_kp[1][2],-(a2-b2)*Q/2.);
  double min=0.0;
  for (size_t i=2;i<nflavs;++i) min+=flavs[i].Mass();
  if ((p[0]+p[1]).Abs2()<min*min) return false;
#ifndef NO_ANALYSIS__ISR_Handler
  m_info_lab[iic::E_1]=p[0][0];
  m_info_lab[iic::t_1]=p[0].Abs2();
  m_info_lab[iic::Em_1]=p[0][0]/p[0].Mass();		
  m_info_lab[iic::z_1]=m_zkey[0][2];
  m_info_lab[iic::E_2]=p[1][0];
  m_info_lab[iic::t_2]=p[1].Abs2();
  m_info_lab[iic::Em_2]=p[1][0]/p[0].Mass();		
  m_info_lab[iic::z_2]=m_zkey[1][2];
#endif
  m_kmrboost=Poincare(p[0]+p[1]);
  if (!m_kmrboost.CheckBoost()) return false;
  m_kmrboost.Boost(p[0]);
  m_kmrboost.Boost(p[1]);
  if (p[1][3]>0.0) m_kmrrot=Poincare(Vec4D::ZVEC,p[1]);
  else m_kmrrot=Poincare(Vec4D::ZVEC,p[0]);
  m_kmrrot.RotateBack(p[0]);
  m_kmrrot.RotateBack(p[1]);
  xi*=xi;
#ifndef NO_ANALYSIS__ISR_Handler
  m_info_cms[iic::E_1]=p[0][0];
  m_info_cms[iic::t_1]=p[0].Abs2();
  m_info_cms[iic::Em_1]=p[0][0]/p[0].Mass();		
  m_info_cms[iic::E_2]=p[1][0];
  m_info_cms[iic::t_2]=p[1].Abs2();
  m_info_cms[iic::Em_2]=p[1][0]/p[0].Mass();		
#endif
  m_weight=(m_x[1]*m_x[0]-b2*b1)/(m_x[1]*m_x[0]);
  m_flux=.25;
  m_flux/=sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
  for (int i=0;i<2;++i) {
    p[n-2+i]=m_fixvecs[i];
    m_kmrboost.Boost(p[n-2+i]);
    m_kmrrot.RotateBack(p[n-2+i]);
    m_kmrboost.Boost(p_kmrlast[i]);
    m_kmrrot.RotateBack(p_kmrlast[i]);
  }
  return true;
}

void ISR_Handler::AssignKeys(ATOOLS::Integration_Info *const info)
{
  m_spkey.Assign("s' isr",4,0,info);
  m_ykey.Assign("y isr",3,0,info);
  m_xkey.Assign("x isr",5,0,info);
  m_zkey[0].Assign("z_1",3,0,info);
  m_zkey[1].Assign("z_2",3,0,info);
  m_kpkey[0].Assign("k_perp_1",4,0,info);
  m_kpkey[1].Assign("k_perp_2",4,0,info);
  m_mu2key[0].Assign("mu2_1",1,0,info);
  m_mu2key[1].Assign("mu2_2",1,0,info);
}

void ISR_Handler::Reset() 
{
  m_splimits[1]=m_fixed_smax*Upper1()*Upper2();
}

void ISR_Handler::SetLimits() 
{
  for (short int i=0;i<3;++i) {
    m_spkey[i]=m_splimits[i];
    if (i<2) m_ykey[i]=m_ylimits[i];
    for (short int j=0;j<2;++j) {
      m_kpkey[j][i]=m_kplimits[i];
      if (i<2) m_zkey[j][i]=m_zlimits[i];
    }
  }
  m_xkey[0]=-0.5*std::numeric_limits<double>::max();
  m_xkey[2]=-0.5*std::numeric_limits<double>::max();
  m_xkey[1]=log(Upper1());
  m_xkey[3]=log(Upper2());
}

bool ISR_Handler::CalculateWeight(const double scale) 
{
  if (!m_kmrmode) {
    m_mu2[0]=m_mu2[1]=scale;
  }
  else {
    m_mu2[0]=m_mu2key[0][0];
    m_mu2[1]=m_mu2key[1][0];
#ifndef NO_ANALYSIS__ISR_Handler
    m_info_cms[iic::mu_1]=m_mu2[0];
    m_info_cms[iic::mu_2]=m_mu2[1];
#endif
  }
  switch (m_mode) {
  case 3 :
#ifdef TEST__Doubly_Unintegrated_PDF
    {
      ATOOLS::Flavour flu(ATOOLS::kf::u);
      ATOOLS::Flavour fld(ATOOLS::kf::d);
      ATOOLS::Flavour fls(ATOOLS::kf::s);
      ATOOLS::Flavour flg(ATOOLS::kf::gluon);
      p_isrbase[0]->CalculateWeight(0.1,0.3,100.,8100.);
      p_isrbase[0]->PDF()->GetXPDF(flu);
      p_isrbase[0]->PDF()->GetXPDF(fld);
      p_isrbase[0]->PDF()->GetXPDF(fls);
      p_isrbase[0]->PDF()->GetXPDF(flu.Bar());
      p_isrbase[0]->PDF()->GetXPDF(fld.Bar());
      p_isrbase[0]->PDF()->GetXPDF(fls.Bar());
      p_isrbase[0]->PDF()->GetXPDF(flg);
      abort();
    }
#endif
    if (p_isrbase[0]->CalculateWeight(m_x[0],m_zkey[0][2],m_kpkey[0][3],m_mu2[0]) && 
	p_isrbase[1]->CalculateWeight(m_x[1],m_zkey[1][2],m_kpkey[1][3],m_mu2[1])) return 1;
    break;
  case 2 :
    if (p_isrbase[1]->CalculateWeight(m_x[1],m_zkey[1][2],m_kpkey[1][3],m_mu2[1])) return 1;
    break;
  case 1 :
    if (p_isrbase[0]->CalculateWeight(m_x[0],m_zkey[0][2],m_kpkey[0][3],m_mu2[0])) return 1;
    break;
  }
  return 0;
}

bool ISR_Handler::CalculateWeight2(const double scale) 
{
  if (m_mode != 3) { 
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"Called for one ISR only.",
			    "ISR_Handler","CalculateWeight2"));
  }
  if (!m_kmrmode) {
    m_mu2[0]=m_mu2[1]=scale;
  }
  else {
    m_mu2[0]=m_mu2key[0][0];
    m_mu2[1]=m_mu2key[1][0];
  }
  if (p_isrbase[0]->CalculateWeight(m_x[1],m_zkey[1][2],m_kpkey[1][3],m_mu2[1]) && 
      p_isrbase[1]->CalculateWeight(m_x[0],m_zkey[0][2],m_kpkey[0][3],m_mu2[0])) { 
    return 1;
  }
  return 0;
}

double ISR_Handler::Weight(const Flavour *const flin)
{
  if (m_mode!=3 || (CheckRemnantKinematics(flin[0],m_xkey[0],0) &&
		    CheckRemnantKinematics(flin[1],m_xkey[1],1))) 
    return p_isrbase[0]->Weight(flin[0])*p_isrbase[1]->Weight(flin[1])/m_weight;
  return 0.;
}

double ISR_Handler::Weight2(const Flavour *const flin)
{
  if (CheckRemnantKinematics(flin[0],m_xkey[0],1) &&
      CheckRemnantKinematics(flin[1],m_xkey[1],0)) 
    return p_isrbase[0]->Weight(flin[1])*p_isrbase[1]->Weight(flin[0])/m_weight;
  return 0.;
}

bool ISR_Handler::BoostInCMS(Vec4D *p,const size_t n) 
{
  for (size_t i=0; i<n; ++i) {
    if (!m_kmrmode) {
      m_cmsboost.Boost(p[i]);
    }
    else {
      m_kmrboost.Boost(p[i]);
      m_kmrrot.RotateBack(p[i]);
    }
  }
  if (m_kmrmode) {
    for (size_t i=0; i<2; ++i) {
      m_kmrboost.Boost(p_kmrlast[i]);
      m_kmrrot.RotateBack(p_kmrlast[i]);
    }
  }
  return true;
}

bool ISR_Handler::BoostInLab(Vec4D* p,const size_t n) 
{
  for (size_t i=0; i<n; ++i) {
    if (!m_kmrmode) {
      m_cmsboost.BoostBack(p[i]);
    }
    else {
      m_kmrrot.Rotate(p[i]);
      m_kmrboost.BoostBack(p[i]);
    }
  }
  if (m_kmrmode) {
    for (size_t i=0; i<2; ++i) {
      m_kmrrot.Rotate(p_kmrlast[i]);
      m_kmrboost.BoostBack(p_kmrlast[i]);
    }
  }
  return true;
}


const Flavour ISR_Handler::DiQuark(const Flavour & fl1,const Flavour & fl2) 
{
  // lightes flavour with that content
  int kf1=fl1.Kfcode();
  int kf2=fl2.Kfcode();
  Flavour diquark;
  if (kf1>kf2) diquark =(kf::code)(kf1*1000 + kf2*100 + 1);
  else if (kf1<kf2)  diquark =(kf::code)(kf2*1000 + kf1*100 + 1);
  else diquark =(kf::code)(kf2*1000 + kf1*100 + 3);
  if (fl1.IsAnti()) diquark=diquark.Bar();
  return diquark;    
}


bool ISR_Handler::CheckRemnantKinematics(const ATOOLS::Flavour & fl,double x,int nbeam)
{
  if (x<.99) return true;
  double mf   = fl.PSMass();
  double msum = 0.;
  double erem = (1. -x -1.e-6)*ATOOLS::rpa.gen.Ecms();
  ATOOLS::Flavour bunch=p_isrbase[nbeam]->Flavour();
  if (!bunch.IsHadron()) return true;
  int hadint=(bunch.Kfcode()-bunch.Kfcode()/10000)/10;
  if ((hadint<=100)||(hadint>=1000)) return true;
  Flavour constit[3];
  constit[0]=ATOOLS::Flavour(ATOOLS::kf::code(hadint)/100);
  constit[1]=ATOOLS::Flavour(ATOOLS::kf::code((hadint-(hadint/100)*100)/10));
  constit[2]=ATOOLS::Flavour(ATOOLS::kf::code(hadint-(hadint/10)*10));
  if (bunch.IsAnti()) {
    for (int i=0;i<3;i++) constit[i]=constit[i].Bar();
  }
  // valence quark
  for (int i=0;i<3;++i) if (constit[i]==fl) {
    for (int j=i+1;j<3;++j) {
      constit[j-1]=constit[j];
    }
    Flavour diquark=DiQuark(constit[0], constit[1]);
    msum+=diquark.PSMass();
    break;
  }
  if (msum==0) {
    // gluon
    if (fl.IsGluon()) {
      Flavour diquark=DiQuark(constit[0], constit[2]);
      msum+=constit[1].PSMass()+diquark.PSMass();    
    }
    else {
      // sea quark
      Flavour diquark=DiQuark(constit[0], constit[2]);
      msum+=constit[1].PSMass()+diquark.PSMass();
      msum+=fl.PSMass();
    }
  }
  if (erem < mf + msum) return false;
  return true;
}

void ISR_Handler::Extract(const ATOOLS::Flavour flavour,const double energy,
			  const size_t i) const 
{ 
  if (p_isrbase[i]->PDF()!=NULL) {
    p_isrbase[i]->Extract(flavour,2.*energy/sqrt(Pole())); 
  }
}

void ISR_Handler::Reset(const size_t i) const 
{ 
  if (p_isrbase[i]->PDF()!=NULL) p_isrbase[i]->Reset(); 
}

ATOOLS::Blob_Data_Base *const ISR_Handler::Info(const int frame) const
{
  if (frame==0) return new ATOOLS::Blob_Data<std::vector<double> >(m_info_cms);
  return new ATOOLS::Blob_Data<std::vector<double> >(m_info_lab);
}

