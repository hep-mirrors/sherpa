#include "PDF/Main/ISR_Handler.H"

#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Base.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "PDF/Remnant/Hadron_Remnant.H"
#include "PDF/Remnant/Electron_Remnant.H"
#include "PDF/Remnant/Photon_Remnant.H"
#include "PDF/Remnant/No_Remnant.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Blob.H"
#include <stdio.h>

using namespace ATOOLS;
using namespace PDF;
using namespace std;

double Lambda2(double sp,double sp1,double sp2) 
{ 
  return (sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
}

ISR_Handler::ISR_Handler(ISR_Base **isrbase):
  p_isrbase(isrbase),
  p_info(new ATOOLS::Integration_Info()),
  m_rmode(0),
  m_weight(1.),
  m_info_lab(8),
  m_info_cms(8)
{
  p_remnants[1]=p_remnants[0]=NULL;
  m_mode=0;
  for (short int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      p_isrbase[i]->AssignKeys(p_info);
      m_mode += i+1;
    }
  }
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());
  m_x[1]=m_x[0]=1.; 
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
#ifdef USING__Threading
  pthread_mutex_init(&m_mtx,NULL);
#endif
}

ISR_Handler::~ISR_Handler() 
{
#ifdef USING__Threading
  pthread_mutex_destroy(&m_mtx);
#endif
  if (p_isrbase) {
    for (int i=0;i<2;i++) {
      if (p_isrbase[i]) delete p_isrbase[i];  
    }
    delete[] p_isrbase; p_isrbase = 0;
  }
  delete p_info;
  for (size_t i(0);i<2;++i) 
    if (p_remnants[i]!=NULL) delete p_remnants[i];
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

void ISR_Handler::SetPDFMember() const
{
  for (int i=0;i<2;++i)
    if (p_isrbase[i]->On()) p_isrbase[i]->PDF()->SetPDFMember();
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *bunches,ATOOLS::Flavour *partons) 
{
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      if (bunches[i] != PDF(i)->Bunch()) { fit = 0; break; }
      fit = PDF(i)->Contains(partons[i]);
      if (fit == 0) break;
    }
    else {
      bool found(false);
      for (size_t j(0);j<p_isrbase[i]->Flavour().Size();++j)
	if (partons[i]==p_isrbase[i]->Flavour()[j]) {
	  found=true;
	  break;
	}
      if (!found) return false;
    }
  }
  return fit;
}

bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *partons) 
{
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) {
      fit = PDF(i)->Contains(partons[i]);
      if (fit == 0) break;
    }
    else {
      bool found(false);
      for (size_t j(0);j<p_isrbase[i]->Flavour().Size();++j)
	if (partons[i]==p_isrbase[i]->Flavour()[j]) {
	  found=true;
	  break;
	}
      if (!found) return false;
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
  p_info->ResetAll();
  m_weight=1.;
  if (m_mode==0) {
    m_x[1]=m_x[0]=1.;
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
  Vec4D pa(p_beam[0]->OutMomentum()), pb(p_beam[1]->OutMomentum());
  double papb(pa*pb), sa(pa.Abs2()), sb(pb.Abs2());
  double gam(papb+sqrt(sqr(papb)-sa*sb));
  double aa(sa/gam), ab(sb/gam), bet(1.0/(1.0-aa*ab));
  Vec4D pp(bet*(pa-aa*pb)), pm(bet*(pb-ab*pa));
  double s1(sqr(flavs[0].Mass())), s2(sqr(flavs[1].Mass()));
  double st(2.0*pp*pm), tau(0.5/st*(m_spkey[3]-s1-s2));
  if (tau*tau<s1*s2/(st*st)) {
    msg_Error()<<METHOD<<"(): s' out of range."<<std::endl;
    return false;
  }
  tau+=sqrt(tau*tau-s1*s2/(st*st));
  if (m_mode==1) {
    m_x[1]=pb.PMinus()/pm.PMinus();
    m_x[0]=tau/m_x[1];
  }
  else if (m_mode==2) {
    m_x[0]=pa.PPlus()/pp.PPlus();
    m_x[1]=tau/m_x[0];
  }
  else if (m_mode==3) {
    double yt(m_ykey[2]-0.5*log((tau+s2/st)/(tau+s1/st)));
    tau=sqrt(tau);
    yt=exp(yt);
    m_x[0]=tau*yt;
    m_x[1]=tau/yt;
  }
  else {
    THROW(fatal_error,"Invalid ISR mode");
  }
  if (PDF(0) && (m_x[0]<PDF(0)->XMin() || m_x[0]>PDF(0)->XMax())) return false;
  if (PDF(1) && (m_x[1]<PDF(1)->XMin() || m_x[1]>PDF(1)->XMax())) return false;
  p[0]=p_cms[0]=m_x[0]*pp+s1/st/m_x[0]*pm;
  p[1]=p_cms[1]=m_x[1]*pm+s2/st/m_x[1]*pp;
  m_cmsboost=Poincare(p_cms[0]+p_cms[1]);
  m_cmsboost.Boost(p[0]);
  m_cmsboost.Boost(p[1]);
  m_flux=0.25/sqrt(sqr(p[0]*p[1])-p[0].Abs2()*p[1].Abs2());
  // m_x[0]=p_cms[0].PPlus()/pa.PPlus();
  // m_x[1]=p_cms[1].PMinus()/pb.PMinus();
  if (m_x[0]>=1.0) m_x[0]=1.0-1.0e-12;
  if (m_x[1]>=1.0) m_x[1]=1.0-1.0e-12;
  return true;
}

void ISR_Handler::AssignKeys(ATOOLS::Integration_Info *const info)
{
  m_spkey.Assign("s' isr",4,0,info);
  m_ykey.Assign("y isr",3,0,info);
  m_beamspkey.Assign("s' beam",4,0,info);
  m_beamykey.Assign("y beam",3,0,info);
  m_xkey.Assign("x isr",5,0,info);
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
  }
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
  if (scale>0.0) {
    m_mu2[0]=m_mu2[1]=scale;
  }
  else {
    m_mu2[0]=m_mu2key[0][0];
    m_mu2[1]=m_mu2key[1][0];
  }
  if (PDF(0) && (m_mu2[0]<PDF(0)->Q2Min() || m_mu2[0]>PDF(0)->Q2Max()))
    return 0.;
  if (PDF(1) && (m_mu2[1]<PDF(1)->Q2Min() || m_mu2[1]>PDF(1)->Q2Max()))
    return 0.;
  switch (m_mode) {
  case 3 :
    if (m_x[0]>p_isrbase[0]->PDF()->RescaleFactor() ||
	m_x[1]>p_isrbase[1]->PDF()->RescaleFactor()) return 0;
    if (p_isrbase[0]->CalculateWeight(m_x[0],0.0,0.0,m_mu2[0]) && 
	p_isrbase[1]->CalculateWeight(m_x[1],0.0,0.0,m_mu2[1])) return 1;
    break;
  case 2 :
    if (m_x[1]>p_isrbase[1]->PDF()->RescaleFactor()) return 0;
    if (p_isrbase[1]->CalculateWeight(m_x[1],0.0,0.0,m_mu2[1])) return 1;
    break;
  case 1 :
    if (m_x[0]>p_isrbase[0]->PDF()->RescaleFactor()) return 0;
    if (p_isrbase[0]->CalculateWeight(m_x[0],0.0,0.0,m_mu2[0])) return 1;
    break;
  case 0 : return 1; 
  }
  return 0;
}

bool ISR_Handler::CalculateWeight2(const double scale) 
{
  if (m_mode != 3) { 
    THROW(fatal_error,"Called for one ISR only.");
  }
  if (scale>0.0) {
    m_mu2[0]=m_mu2[1]=scale;
  }
  else {
    m_mu2[0]=m_mu2key[0][0];
    m_mu2[1]=m_mu2key[1][0];
  }
  if (PDF(0) && (m_mu2[1]<PDF(0)->Q2Min() || m_mu2[1]>PDF(0)->Q2Max()))
    return 0.;
  if (PDF(1) && (m_mu2[0]<PDF(1)->Q2Min() || m_mu2[0]>PDF(1)->Q2Max()))
    return 0.;
  if (m_x[0]>p_isrbase[1]->PDF()->RescaleFactor() ||
      m_x[1]>p_isrbase[0]->PDF()->RescaleFactor()) return 0;
  if (p_isrbase[0]->CalculateWeight(m_x[1],0.0,0.0,m_mu2[1]) && 
      p_isrbase[1]->CalculateWeight(m_x[0],0.0,0.0,m_mu2[0])) { 
    return 1;
  }
  return 0;
}

double ISR_Handler::Weight(const Flavour *const flin)
{
  if (m_mode!=3 || (CheckRemnantKinematics(flin[0],m_x[0],0,false) &&
		    CheckRemnantKinematics(flin[1],m_x[1],1,false))) {
    m_xf1[0]=p_isrbase[0]->Weight(flin[0]);
    m_xf2[0]=p_isrbase[1]->Weight(flin[1]);
    return m_xf1[0]*m_xf2[0]/m_weight;
  }
  return 0.;
}

double ISR_Handler::Weight2(const Flavour *const flin)
{
  if (CheckRemnantKinematics(flin[0],m_x[0],1,true) &&
      CheckRemnantKinematics(flin[1],m_x[1],0,true)) {
    m_xf1[1]=p_isrbase[0]->Weight(flin[1]);
    m_xf2[1]=p_isrbase[1]->Weight(flin[0]);
    return m_xf1[1]*m_xf2[1]/m_weight;
  }
  return 0.;
}

bool ISR_Handler::BoostInCMS(Vec4D *p,const size_t n) 
{
  for (size_t i=0; i<n; ++i) {
    m_cmsboost.Boost(p[i]);
  }
  return true;
}

bool ISR_Handler::BoostInLab(Vec4D* p,const size_t n) 
{
  for (size_t i=0; i<n; ++i) {
    m_cmsboost.BoostBack(p[i]);
  }
  return true;
}

bool ISR_Handler::CheckRemnantKinematics(const ATOOLS::Flavour &fl,
					 double &x,int beam,bool swaped)
{
  if (x>p_isrbase[beam]->PDF()->RescaleFactor()) return false;
  if (m_rmode==0) return true;
  p_remnants[beam]->QuickClear();
  double pp(beam==0?x*p_beam[0]->OutMomentum().PPlus():
	    x*p_beam[1]->OutMomentum().PMinus());
  double pm(sqr(fl.Mass()));
  pm/=pp;
  Vec4D mom((pp+pm)/2.0,0.0,0.0,beam==0?(pp-pm)/2.0:(pm-pp)/2.0);
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

void ISR_Handler::MtxLock()
{
#ifdef USING__Threading
  pthread_mutex_lock(&m_mtx);
#endif
}

void ISR_Handler::MtxUnLock()
{
#ifdef USING__Threading
  pthread_mutex_unlock(&m_mtx);
#endif
}
