#include "ISR_Handler.H"

#include "Beam_Base.H"
#include "Intact.H"
#include "Structure_Function.H"
#include "Doubly_Unintegrated_PDF.H"
#include "Run_Parameter.H" 
#include "Hadron_Remnant.H"
#include "Electron_Remnant.H"
#include "Photon_Remnant.H"
#include "No_Remnant.H"
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

#ifdef PROFILE__all
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

ISR_Handler::ISR_Handler(ISR_Base **isrbase):
  p_isrbase(isrbase),
  p_info(new ATOOLS::Integration_Info()),
  m_kperpscheme((int)kps::constant),
  m_rmode(0),
  m_weight(1.),
  m_info_lab(8),
  m_info_cms(8)
{
  p_remnants[1]=p_remnants[0]=NULL;
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
  Doubly_Unintegrated_PDF *dupdf=
    dynamic_cast<Doubly_Unintegrated_PDF*>(p_isrbase[0]->PDF());
  if (dupdf!=NULL) m_kperpscheme=dupdf->KPerpScheme();
  for (size_t i=0;i<2;++i) {
    if (Flav(i).IsHadron()) {
      Hadron_Remnant *remnant = new Hadron_Remnant(this,i);
      remnant->SetStringDrawing(1.0,0);
      remnant->SetStringDrawing(0.0,1);
      p_remnants[i]=remnant;
    }
    else if (Flav(i).IsLepton()) 
      p_remnants[i] = new Electron_Remnant(this,i);
    else if (Flav(i).IsPhoton()) 
      p_remnants[i] = new Photon_Remnant(i);
    else p_remnants[i] = new No_Remnant(i);
  }
  for (size_t i=0;i<2;++i) p_remnants[i]->SetPartner(p_remnants[1-i]);
#ifdef TEST__ISR_Handler
  TestPDF testpdf;
  testpdf.p_pdf=p_isrbase[0]->PDF();
  testpdf.m_mu2=8100.;
  Gauss_Integrator gauss(&testpdf);
  std::ofstream *out = new std::ofstream("pdftest");
  for (double test=-1.;test<=4.;test+=0.1) {
    testpdf.m_z=0.8;
    testpdf.m_kt2=exp(test);
#ifdef USING__KPerp
#define PREFACTOR exp(test)*
#else
#define PREFACTOR 
#endif
    testpdf.m_x=0.1;
    *out<<exp(test)<<" "<<PREFACTOR gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.01;
    *out<<" "<<PREFACTOR gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.001;
    *out<<" "<<PREFACTOR gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.0001;
    *out<<" "<<PREFACTOR gauss.Integrate(testpdf.m_x,1.,1.e-3);
    testpdf.m_x=0.1;
    *out<<exp(test)<<" "<<PREFACTOR testpdf();
    testpdf.m_x=0.01;
    *out<<" "<<PREFACTOR testpdf();
    testpdf.m_x=0.001;
    *out<<" "<<PREFACTOR testpdf();
    testpdf.m_x=0.0001;
    *out<<" "<<PREFACTOR testpdf()<<std::endl;
  }
  delete out;
  THROW(normal_exit,"finished integration");
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

void ISR_Handler::Init(double *splimits,double *kplimits) 
{
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());

  double s=(p_beam[0]->OutMomentum()+
	    p_beam[1]->OutMomentum()).Abs2();
  ATOOLS::rpa.gen.SetEcms(sqrt(s));

  m_type = p_isrbase[0]->Type()+std::string("*")+p_isrbase[1]->Type();
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
  double E=ATOOLS::rpa.gen.Ecms();
  double x=1./2.+(m_mass2[0]-m_mass2[1])/(2.*E*E);
  double E1=x*E;
  double E2=E-E1;
  p_remnants[0]->SetBeam(p_beam[0]);
  p_remnants[1]->SetBeam(p_beam[1]);
  m_fixvecs[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass2[0]));
  m_fixvecs[1]=Vec4D(E2,0.,0.,-m_fixvecs[0][3]);
}

void ISR_Handler::SetSprimeMin(const double spmin)       
{ 
  m_spkey[0]=m_splimits[0]=Max(m_fixed_smin,spmin); 
}

void ISR_Handler::SetSprimeMax(const double spmax)       
{ 
  m_spkey[1]=m_splimits[1]=Min(m_fixed_smax,spmax); 
}

void ISR_Handler::SetFixedSprimeMin(const double spmin)  
{ 
  m_fixed_smin = spmin;
  m_splimits[0] = spmin;
}

void ISR_Handler::SetFixedSprimeMax(const double spmax)  
{
  m_fixed_smax = spmax;
  m_splimits[1] = spmax;
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

void ISR_Handler::SetMasses(const Flavour *fl,const size_t nout) 
{
  m_mass2[0]=sqr(fl[0].Mass());
  m_mass2[1]=sqr(fl[1].Mass());
  double emin=0.0;
  for (size_t i=0;i<nout;++i) emin+=fl[2+i].Mass();
  emin=ATOOLS::Max(emin,fl[0].Mass()+fl[1].Mass());
  m_splimits[0]=ATOOLS::Max(m_splimits[0],sqr(emin));
}

void ISR_Handler::SetPartonMasses(Flavour *fl) 
{
  SetMasses(fl);
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
    m_flux=0.25/sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
    return true;
  }
  if (m_spkey[3]<m_splimits[0] || m_spkey[3]>m_splimits[1]) {
    msg_Error()<<METHOD<<"(..): "<<om::red
		       <<"s' out of bounds.\n"<<om::reset
		       <<"  s'_{min}, s'_{max 1,2} vs. s': "<<m_splimits[0]
		       <<", "<<m_splimits[1]<<", "<<m_splimits[2]
		       <<" vs. "<<m_spkey[3]<<std::endl;
    return false;
  }
  if (m_ykey[2]<m_ykey[0] || m_ykey[2]>m_ykey[1]) {
    msg_Error()<<METHOD<<"(..): "<<om::red
		       <<"y out of bounds.\n"<<om::reset
		       <<"  y_{min}, y_{max} vs. y: "<<m_ykey[0]
		       <<", "<<m_ykey[1]<<" vs. "<<m_ykey[2]<<std::endl;
    return false;
  }
  double Q=sqrt(m_splimits[2]), E=sqrt(m_spkey[3]);
  double E1=(m_spkey[3]+sqr(flavs[0].Mass())-sqr(flavs[1].Mass()))/(2.0*E);
  p_cms[0]=p[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-sqr(flavs[0].Mass())));
  p_cms[1]=p[1]=Vec4D(E-E1,(-1.)*Vec3D(p[0]));
  Vec4D plab[2]; plab[0]=p[0]; plab[1]=p[1];
  m_cmsboost=Poincare(Vec4D(cosh(m_ykey[2]),0.,0.,sinh(m_ykey[2])));
  m_cmsboost.BoostBack(p_cms[0]);
  m_cmsboost.BoostBack(p_cms[1]);
  m_x[0]=p_cms[0].PPlus()/Q;
  m_x[1]=p_cms[1].PMinus()/Q;
  m_flux=0.25/sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
  if (!m_kmrmode) {
#ifndef NO_ANALYSIS__ISR_Handler
    m_info_lab[iic::E_1]=m_x[0];
    m_info_lab[iic::E_2]=m_x[1];
#endif
    return true;
  }
  E=sqrt(m_spkey[3]-(m_kpkey[0](0)+m_kpkey[1](0)).Abs2());
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
  p[0]=Vec4D((m_x[0]-b1)*Q/2.,m_kpkey[0](0)[1],
	     m_kpkey[0](0)[2],(m_x[0]+b1)*Q/2.);
  p[1]=Vec4D((m_x[1]-b2)*Q/2.,m_kpkey[1](0)[1],
	     m_kpkey[1](0)[2],-(m_x[1]+b2)*Q/2.);
  double a1=m_zkey[0][2]==0.0?0.0:m_x[0]*(1./m_zkey[0][2]-1.);
  p_kmrlast[0]=Vec4D((a1+b1)*Q/2.,-m_kpkey[0](0)[1],
		     -m_kpkey[0](0)[2],(a1-b1)*Q/2.);
  double a2=m_zkey[1][2]==0.0?0.0:m_x[1]*(1./m_zkey[1][2]-1.);
  p_kmrlast[1]=Vec4D((a2+b2)*Q/2.,-m_kpkey[1](0)[1],
		     -m_kpkey[1](0)[2],-(a2-b2)*Q/2.);
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
#ifndef NO_ANALYSIS__ISR_Handler
  m_info_cms[iic::E_1]=p[0][0];
  m_info_cms[iic::t_1]=p[0].Abs2();
  m_info_cms[iic::Em_1]=p[0][0]/p[0].Mass();		
  m_info_cms[iic::E_2]=p[1][0];
  m_info_cms[iic::t_2]=p[1].Abs2();
  m_info_cms[iic::Em_2]=p[1][0]/p[0].Mass();		
#endif
  m_weight=(m_x[1]*m_x[0]-b2*b1)/(m_x[1]*m_x[0]);
  // for on-shell kinematics only !
  //   m_flux=0.5/sqrt(sqr(m_spkey[3]-m_mass2[0]-m_mass2[1])
  // 		  -4.0*m_mass2[0]*m_mass2[1]);
  m_flux=sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2();
  if (m_flux>0.0) m_flux=0.25/sqrt(m_flux);
  else return false;
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
  m_beamspkey.Assign("s' beam",4,0,info);
  m_beamykey.Assign("y beam",3,0,info);
  m_xkey.Assign("x isr",5,0,info);
  m_zkey[0].Assign("z_1",3,0,info);
  m_zkey[1].Assign("z_2",3,0,info);
  m_kpkey[0].Assign("k_perp_1",4,1,info);
  m_kpkey[1].Assign("k_perp_2",4,1,info);
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
  m_kpkey[1](0)=m_kpkey[0](0)=Vec4D();
  m_xkey[0]=m_mass2[0]==0.0?-0.5*std::numeric_limits<double>::max():
    log(m_mass2[0]/sqr(p_beam[0]->OutMomentum().PPlus()));
  m_xkey[2]=m_mass2[1]==0.0?-0.5*std::numeric_limits<double>::max():
    log(m_mass2[1]/sqr(p_beam[1]->OutMomentum().PMinus()));
  double e1=p_beam[0]->OutMomentum()[0];
  m_xkey[1]=ATOOLS::Min(e1/p_beam[0]->OutMomentum().PPlus()*
			(1.0+sqrt(1.0-m_mass2[0]/sqr(e1))),Upper1());
  double e2=p_beam[1]->OutMomentum()[0];
  m_xkey[3]=ATOOLS::Min(e2/p_beam[1]->OutMomentum().PMinus()*
			(1.0+sqrt(1.0-m_mass2[1]/sqr(e2))),Upper2());
  m_spkey[1]=m_splimits[1]=Min(m_splimits[1],m_splimits[2]*m_xkey[1]*m_xkey[3]);
  m_xkey[1]=log(m_xkey[1]);
  m_xkey[3]=log(m_xkey[3]);
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
      double x=0.01, z=0.9999, kp2=4000.0, mu2=8100.0;
      p_isrbase[0]->CalculateWeight(x,z,kp2,mu2);
      PRINT_INFO("x = "<<x<<", z = "<<z<<", kt2 = "<<kp2<<", mu2 = "<<mu2);
      PRINT_INFO("f_u  = "<<(p_isrbase[0]->PDF()->GetXPDF(flu)*kp2));
      PRINT_INFO("f_d  = "<<(p_isrbase[0]->PDF()->GetXPDF(fld)*kp2));
      PRINT_INFO("f_s  = "<<(p_isrbase[0]->PDF()->GetXPDF(fls)*kp2));
      PRINT_INFO("f_ub = "<<(p_isrbase[0]->PDF()->GetXPDF(flu.Bar())*kp2));
      PRINT_INFO("f_db = "<<(p_isrbase[0]->PDF()->GetXPDF(fld.Bar())*kp2));
      PRINT_INFO("f_sb = "<<(p_isrbase[0]->PDF()->GetXPDF(fls.Bar())*kp2));
      PRINT_INFO("f_g  = "<<(p_isrbase[0]->PDF()->GetXPDF(flg)*kp2));
      PDF_Base *pdf=p_isrbase[0]->PDF()->GetBasicPDF();
      pdf->Calculate(x,0.0,0.0,mu2);
      PRINT_INFO("xu  = "<<pdf->GetXPDF(flu));
      PRINT_INFO("xd  = "<<pdf->GetXPDF(fld));
      PRINT_INFO("xs  = "<<pdf->GetXPDF(fls));
      PRINT_INFO("xub = "<<pdf->GetXPDF(flu.Bar()));
      PRINT_INFO("xdb = "<<pdf->GetXPDF(fld.Bar()));
      PRINT_INFO("xsb = "<<pdf->GetXPDF(fls.Bar()));
      PRINT_INFO("xg  = "<<pdf->GetXPDF(flg));
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
    THROW(fatal_error,"Called for one ISR only.");
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
  if (m_mode!=3 || (CheckRemnantKinematics(flin[0],m_x[0],0,false) &&
		    CheckRemnantKinematics(flin[1],m_x[1],1,false)))
    return p_isrbase[0]->Weight(flin[0])*p_isrbase[1]->Weight(flin[1])
      /m_weight;
  return 0.;
}

double ISR_Handler::Weight2(const Flavour *const flin)
{
  if (CheckRemnantKinematics(flin[0],m_x[0],1,true) &&
      CheckRemnantKinematics(flin[1],m_x[1],0,true)) 
    return p_isrbase[0]->Weight(flin[1])*p_isrbase[1]->Weight(flin[0])
      /m_weight;
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

bool ISR_Handler::CheckRemnantKinematics(const ATOOLS::Flavour &fl,
					 double &x,int beam,bool swaped)
{
  PROFILE_HERE;
  if (m_rmode==0) return true;
  p_remnants[beam]->QuickClear();
  double pp(beam==0?x*p_beam[0]->OutMomentum().PPlus():
	    x*p_beam[1]->OutMomentum().PMinus());
  double pm(sqr(fl.Mass()));
  if (m_kmrmode>0) pm=-m_kpkey[beam][3]/(1.0-m_zkey[beam][2])+m_kpkey[beam][3];
  pm/=pp;
  Vec4D mom((pp+pm)/2.0,m_kpkey[beam](0)[1],
	    m_kpkey[beam](0)[2],beam==0?(pp-pm)/2.0:(pm-pp)/2.0);
  if (m_kmrmode>0) {
    mom+=p_kmrlast[swaped?1-beam:beam];
    if (mom[0]<0.0 || mom[0]>p_beam[beam]->OutMomentum()[0]) return false;
  }
  return p_remnants[beam]->TestExtract(fl,mom);
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

