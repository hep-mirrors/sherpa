#include "PDF/Main/ISR_Handler.H"

#include "BEAM/Main/Beam_Base.H"
#include "PDF/Main/ISR_Base.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "PDF/Remnant/Hadron_Remnant.H"
#include "PDF/Remnant/Electron_Remnant.H"
#include "PDF/Remnant/Photon_Remnant.H"
#include "PDF/Remnant/No_Remnant.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;
using namespace PDF;
using namespace std;

static int s_nozeropdf=-1;

double Lambda2(double sp,double sp1,double sp2) 
{ 
  return (sp-sp1-sp2)*(sp-sp1-sp2)-4.0*sp1*sp2;
}

ISR_Handler::ISR_Handler(ISR_Base **isrbase):
  p_isrbase(isrbase),
  m_rmode(0),
  m_info_lab(8),
  m_info_cms(8)
{
  if (s_nozeropdf<0) {
    Data_Reader dr(" ",";","!","=");
    dr.AddComment("#");
    dr.AddWordSeparator("\t");
    dr.SetInputPath(rpa->GetPath());
    dr.SetInputFile(rpa->gen.Variable("INTEGRATION_DATA_FILE"));
    s_nozeropdf=dr.GetValue<int>("NO_ZERO_PDF",0);
  }
  m_mu2[0]=m_mu2[1]=0.0;
  m_xf1[0]=m_xf2[0]=m_xf1[1]=m_xf2[1]=1.0;
  p_remnants[1]=p_remnants[0]=NULL;
  m_mode=0;
  for (short int i=0;i<2;i++) {
    if (p_isrbase[i]->On()) m_mode += i+1;
  }
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());
  m_x[1]=m_x[0]=1.; 
  m_mu2[1]=m_mu2[0]=0.0; 
  m_xf1[1]=m_xf1[0]=0.0;
  m_xf2[1]=m_xf2[0]=0.0;
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
  for (size_t i(0);i<2;++i) 
    if (p_remnants[i]!=NULL) delete p_remnants[i];
}

void ISR_Handler::Init(double *splimits) 
{
  m_mass2[0]=sqr(p_isrbase[0]->Flavour().Mass());
  m_mass2[1]=sqr(p_isrbase[1]->Flavour().Mass());

  double s=(p_beam[0]->OutMomentum()+
	    p_beam[1]->OutMomentum()).Abs2();
  ATOOLS::rpa->gen.SetEcms(sqrt(s));

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
  double E=ATOOLS::rpa->gen.Ecms();
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
  m_splimits[0]=Max(m_fixed_smin,spmin); 
}

void ISR_Handler::SetSprimeMax(const double spmax)       
{ 
  m_splimits[1]=Min(m_fixed_smax,spmax); 
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
      if (fit == 0) {
	for (size_t j(0);j<partons[i].Size();++j)
	  fit |= PDF(i)->Contains(partons[i][j]);
      }
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

void ISR_Handler::SetMasses(const Flavour_Vector &fl) 
{
  m_mass2[0]=sqr(fl[0].Mass());
  m_mass2[1]=sqr(fl[1].Mass());
  double emin=0.0;
  for (size_t i=2;i<fl.size();++i) emin+=fl[i].Mass();
  emin=ATOOLS::Max(emin,fl[0].Mass()+fl[1].Mass());
  m_splimits[0]=ATOOLS::Max(m_splimits[0],sqr(emin));
}

void ISR_Handler::SetPartonMasses(const Flavour_Vector &fl) 
{
  SetMasses(fl);
  double E=ATOOLS::rpa->gen.Ecms();
  double x=1./2.+(m_mass2[0]-m_mass2[1])/(2.*E*E);
  double E1=x*E;
  double E2=E-E1;
  m_fixvecs[0]=Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass2[0]));
  m_fixvecs[1]=Vec4D(E2,0.,0.,-m_fixvecs[0][3]);
}

bool ISR_Handler::MakeISR(const double &sp,const double &y,
			  Vec4D_Vector &p,const Flavour_Vector &flavs) 
{
  if (m_mode==0) {
    m_x[1]=m_x[0]=1.;
    return true;
  }
  if (sp<m_splimits[0] || sp>m_splimits[1]) {
    msg_Error()<<METHOD<<"(..): "<<om::red
		       <<"s' out of bounds.\n"<<om::reset
		       <<"  s'_{min}, s'_{max 1,2} vs. s': "<<m_splimits[0]
		       <<", "<<m_splimits[1]<<", "<<m_splimits[2]
		       <<" vs. "<<sp<<std::endl;
    return false;
  }
  Vec4D pa(p_beam[0]->OutMomentum()), pb(p_beam[1]->OutMomentum());
  double papb(pa*pb), sa(pa.Abs2()), sb(pb.Abs2());
  double gam(papb+sqrt(sqr(papb)-sa*sb));
  double aa(sa/gam), ab(sb/gam), bet(1.0/(1.0-aa*ab));
  Vec4D pp(bet*(pa-aa*pb)), pm(bet*(pb-ab*pa));
  double s1(sqr(flavs[0].Mass())), s2(sqr(flavs[1].Mass()));
  double st(2.0*pp*pm), tau(0.5/st*(sp-s1-s2));
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
    double yt(y-0.5*log((tau+s2/st)/(tau+s1/st)));
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
  // m_x[0]=p_cms[0].PPlus()/pa.PPlus();
  // m_x[1]=p_cms[1].PMinus()/pb.PMinus();
  if (m_x[0]>=1.0) m_x[0]=1.0-1.0e-12;
  if (m_x[1]>=1.0) m_x[1]=1.0-1.0e-12;
  return true;
}

void ISR_Handler::Reset() 
{
  m_splimits[1]=m_fixed_smax*Upper1()*Upper2();
}

void ISR_Handler::SetLimits(Double_Vector &spkey,Double_Vector &ykey,
			    Double_Vector &xkey) 
{
  for (short int i=0;i<3;++i) {
    spkey[i]=m_splimits[i];
    if (i<2) ykey[i]=m_ylimits[i];
  }
  xkey[0]=m_mass2[0]==0.0?-0.5*std::numeric_limits<double>::max():
    log(m_mass2[0]/sqr(p_beam[0]->OutMomentum().PPlus()));
  xkey[2]=m_mass2[1]==0.0?-0.5*std::numeric_limits<double>::max():
    log(m_mass2[1]/sqr(p_beam[1]->OutMomentum().PMinus()));
  double e1=p_beam[0]->OutMomentum()[0];
  xkey[1]=ATOOLS::Min(e1/p_beam[0]->OutMomentum().PPlus()*
		      (1.0+sqrt(1.0-m_mass2[0]/sqr(e1))),Upper1());
  double e2=p_beam[1]->OutMomentum()[0];
  xkey[3]=ATOOLS::Min(e2/p_beam[1]->OutMomentum().PMinus()*
		      (1.0+sqrt(1.0-m_mass2[1]/sqr(e2))),Upper2());
  spkey[1]=m_splimits[1]=Min(m_splimits[1],m_splimits[2]*xkey[1]*xkey[3]);
  xkey[1]=log(xkey[1]);
  xkey[3]=log(xkey[3]);
}

double ISR_Handler::Weight(const int mode,Vec4D p1,Vec4D p2,
			   double Q12,double Q22,Flavour fl1,Flavour fl2,int warn)
{
  if (m_mode==0) return 0.25/sqrt(sqr(p1*p2)-p1.Abs2()*p2.Abs2());
  msg_IODebugging()<<METHOD<<"(mode = "<<mode<<")\n";
  if (fl1.Size()>1 || fl2.Size()>1)
    THROW(fatal_error,"Do not try to calculate an ISR weight with containers.");
  double x1(0.),x2(0.);
  if (mode) {
    p1[3]=-p1[3];
    p2[3]=-p2[3];
    std::swap<Flavour>(fl1,fl2);
    std::swap<Vec4D>(p1,p2);
    std::swap<double>(Q12,Q22);
  }
  x1=Min(1.0,p1.PPlus()/p_beam[0]->OutMomentum().PPlus());
  x2=Min(1.0,p2.PMinus()/p_beam[1]->OutMomentum().PMinus());
  msg_IODebugging()<<"  "<<p1<<" from "<<p_beam[0]->OutMomentum()<<" -> "
		 <<p1.PPlus()<<" / "<<p_beam[0]->
    OutMomentum().PPlus()<<" = "<<x1<<std::endl;
  msg_IODebugging()<<"  "<<p2<<" from "<<p_beam[1]->OutMomentum()<<" -> "
		 <<p2.PMinus()<<" / "<<p_beam[1]->
    OutMomentum().PMinus()<<" = "<<x2<<std::endl;
  if (PDF(0) && (Q12<PDF(0)->Q2Min() || Q12>PDF(0)->Q2Max()))
    return 0.;
  if (PDF(1) && (Q22<PDF(1)->Q2Min() || Q22>PDF(1)->Q2Max()))
    return 0.;
  MtxLock();
  m_mu2[mode]=Q12;
  m_mu2[1-mode]=Q22;
  switch (m_mode) {
    case 3 :
      if (!p_isrbase[0]->PDF()->Contains(fl1) ||
          !p_isrbase[1]->PDF()->Contains(fl2)) { MtxUnLock(); return 0.; }
      if (x1>p_isrbase[0]->PDF()->RescaleFactor() ||
          x2>p_isrbase[1]->PDF()->RescaleFactor()) { MtxUnLock(); return 0.; }
      if (p_isrbase[0]->CalculateWeight(x1,0.0,0.0,Q12,warn) &&
          p_isrbase[1]->CalculateWeight(x2,0.0,0.0,Q22,warn)) break;
    case 2 :
      if (!p_isrbase[1]->PDF()->Contains(fl2)) { MtxUnLock(); return 0.; }
      if (m_x[1]>p_isrbase[1]->PDF()->RescaleFactor()) { MtxUnLock(); return 0.; }
      if (p_isrbase[1]->CalculateWeight(x2,0.0,0.0,Q22,warn)) break;
    case 1 :
      if (!p_isrbase[0]->PDF()->Contains(fl1)) { MtxUnLock(); return 0.; }
      if (m_x[0]>p_isrbase[0]->PDF()->RescaleFactor()) { MtxUnLock(); return 0.; }
      if (p_isrbase[0]->CalculateWeight(x1,0.0,0.0,Q12,warn)) break;
    case 0 : break;
    default : MtxUnLock(); return 0.;
  }
  if (m_mode!=3 || (CheckRemnantKinematics(fl1,x1,0,false) &&
                    CheckRemnantKinematics(fl2,x2,1,false))) {
    double f1=p_isrbase[0]->Weight(fl1), f2=p_isrbase[1]->Weight(fl2);
    m_xf1[mode]=x1*f1;
    m_xf2[mode]=x2*f2;
    MtxUnLock();
    msg_IODebugging()<<"  PDF1: "<<rpa->gen.Beam1()<<" -> "<<fl1<<" at ("<<x1
		   <<","<<sqrt(Q12)<<") -> "<<om::bold<<f1<<om::reset<<"\n";
    msg_IODebugging()<<"  PDF2: "<<rpa->gen.Beam2()<<" -> "<<fl2<<" at ("<<x2
		   <<","<<sqrt(Q22)<<") -> "<<om::bold<<f2<<om::reset<<"\n";
    double flux=0.25/sqrt(sqr(p1*p2)-p1.Abs2()*p2.Abs2());
    msg_IODebugging()<<"  Flux: "<<flux<<", Weight: "<<f1*f2*flux<<std::endl;
    if (IsBad(f1*f2)) return 0.0;
    if (s_nozeropdf && f1*f2==0.0)
      return pow(std::numeric_limits<double>::min(),0.25);
    return f1*f2*flux;
  }
  MtxUnLock();
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

ATOOLS::Blob_Data_Base* ISR_Handler::Info(const int frame) const
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
