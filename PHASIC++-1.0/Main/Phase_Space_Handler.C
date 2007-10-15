#include "Phase_Space_Handler.H"

#include "Phase_Space_Integrator.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "CXXFLAGS.H"
#include "Rambo.H"
#include "RamboKK.H"
#include "Sarge.H"
#include "VHAAG.H"
#include "Leading_Log_Z.H"
#include "LL_KPerp.H"
#include "FSR_Channel.H"
#include "ISR_Vegas.H"
#include "PI_Interface.H"
#include "Running_AlphaS.H"
#include "Permutation.H"

#include "Run_Parameter.H"
#include "Blob.H"
#include "Message.H"  
#include "Random.H"
#include "Shell_Tools.H"
#include "MyStrStream.H"
#include "Data_Reader.H"

#include <dlfcn.h>

#ifdef PROFILE__all
#define PROFILE__Phase_Space_Handler
#endif
#ifdef PROFILE__Phase_Space_Handler
#include "prof.hh"
#else
#define PROFILE_HERE 
#define PROFILE_LOCAL(LOCALNAME)
#endif

#ifdef USING__ROOT
// #define ANALYSE__Phase_Space_Handler
#endif
#ifdef ANALYSE__Phase_Space_Handler
#include "My_Root.H"
#include "TH2D.h"
#include "TH3D.h"
static ATOOLS::Info_Key m_isrzkey[2];
class PS_Histogram_2D: public TH2D {
private:
  TH2D *p_ps;
public:
  // constructor
  PS_Histogram_2D(const char *name,const char *title,
		  const size_t nbinsx,const double xmin,const double xmax,
		  const size_t nbinsy,const double ymin,const double ymax):
    TH2D(name,title,nbinsx,xmin,xmax,nbinsy,ymin,ymax),
    p_ps(new TH2D((std::string(name)+"_ps").c_str(),
		  (std::string(title)+"_ps").c_str(),
		  nbinsx,xmin,xmax,nbinsy,ymin,ymax)) {}
  // destructor
  ~PS_Histogram_2D() {}
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
};// end of class PS_Histogram_2D
class PS_Histogram_3D: public TH3D {
private:
  TH3D *p_ps;
public:
  // constructor
  PS_Histogram_3D(const char *name,const char *title,
	       const size_t nbinsx,const double xmin,const double xmax,
	       const size_t nbinsy,const double ymin,const double ymax,
	       const size_t nbinsz,const double zmin,const double zmax):
    TH3D(name,title,nbinsx,xmin,xmax,nbinsy,ymin,ymax,nbinsz,zmin,zmax),
    p_ps(new TH3D((std::string(name)+"_ps").c_str(),
		  (std::string(title)+"_ps").c_str(),
		  nbinsx,xmin,xmax,nbinsy,ymin,ymax,nbinsz,zmin,zmax)) {}
  // destructor
  ~PS_Histogram_3D() {}
  // member functions
  Int_t Fill(const Double_t x,const Double_t y,const Double_t z,
	     const Double_t weight)
  {
    p_ps->Fill(x,y,z,1.0);
    return TH3D::Fill(x,y,z,weight);
  }
  void Draw(Option_t *option="")
  {
    for (Int_t i=0;i<GetNbinsX();++i)
      for (Int_t j=0;j<GetNbinsY();++j)
	for (Int_t k=0;k<GetNbinsZ();++k)
	  SetBinContent(i,j,k,p_ps->GetBinContent(i,j,k)==0.0?0.0:
			GetBinContent(i,j,k)/
			p_ps->GetBinContent(i,j,k));
    TVirtualPad *psc=gPad;
    Int_t logx=psc->GetLogx();
    Int_t logy=psc->GetLogy();
    Int_t logz=psc->GetLogz();
    psc->Divide(2,1);
    psc->cd(1);
    gPad->SetLogx(logx);
    gPad->SetLogy(logy);
    gPad->SetLogz(logz);
    TH3D::Draw(option);
    psc->cd(2);
    gPad->SetLogx(logx);
    gPad->SetLogy(logy);
    gPad->SetLogz(logz);
    p_ps->Draw(option);
  }
};// end of class PS_Histogram_3D
typedef std::map<PHASIC::Integrable_Base*,
		 PS_Histogram_2D*> Analysis_Map;
static Analysis_Map s_psspy, s_psz, s_pskp;
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace BEAM;
using namespace PDF;

Integration_Info *PHASIC::Phase_Space_Handler::p_info=NULL;

Phase_Space_Handler::Phase_Space_Handler(Integrable_Base *proc,
					 ISR_Handler *ih,Beam_Spectra_Handler *bh, double error): 
  m_name(proc->Name()), p_process(proc), p_active(proc), p_integrator(NULL), p_cuts(NULL),
  p_beamhandler(bh), p_isrhandler(ih), p_fsrchannels(NULL), p_zchannels(NULL), p_kpchannels(NULL), 
  p_isrchannels(NULL), p_beamchannels(NULL), p_flavours(NULL), p_cms(NULL), p_lab(NULL), p_massboost(NULL),
  m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(0), m_use_pi(0), m_initialized(0),
  m_maxtrials(1000000), m_sumtrials(0), m_events(0), m_E(ATOOLS::rpa.gen.Ecms()), m_s(m_E*m_E), 
  m_weight(1.)
{
  p_pi=NULL;
  Data_Read dr(rpa.GetPath()+"/Integration.dat");
  m_error    = dr.GetValue<double>("ERROR",0.01);
  m_inttype  = dr.GetValue<int>("INTEGRATOR",3);
  m_fin_opt  = dr.GetValue<Switch::code>("FINISH_OPTIMIZATION");
  if (m_fin_opt==NotDefined<Switch::code>()) m_fin_opt=Switch::Off;
  if (error>0.) {
    m_error   = error;
    m_fin_opt = Switch::Off;
  }
  p_flavours = new Flavour[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) p_flavours[i] = proc->Flavours()[i];
  p_channellibnames = new std::list<std::string>;
  p_fsrchannels = new Multi_Channel("fsr_"+proc->Name());
  m_m[0] = p_flavours[0].Mass(); m_m2[0] = m_m[0]*m_m[0];
  if (m_nin==2) {
    m_m[1] = p_flavours[1].Mass(); m_m2[1] = m_m[1]*m_m[1]; 
    if (p_beamhandler) {
      if (p_beamhandler->On()>0) p_beamchannels = new Multi_Channel("beam_"+proc->Name());
    }
    if (p_isrhandler) {
      if (p_isrhandler->On()>0) {
	p_isrchannels = new Multi_Channel("isr_"+proc->Name());
 	if (p_isrhandler->KMROn()>0) {
 	  p_zchannels = new Multi_Channel("kmr_z_"+proc->Name());
 	  p_kpchannels = new Multi_Channel("kmr_kp_"+proc->Name());
	}
      }
    }
  }
  if (m_nin==2) {
    m_isrspkey.Assign("s' isr",4,0,p_info);
    m_isrykey.Assign("y isr",3,0,p_info);
    m_beamspkey.Assign("s' beam",4,0,p_info);
    m_beamykey.Assign("y beam",3,0,p_info);
    m_mu2key[0].Assign("mu2_1",1,0,p_info);
    m_mu2key[1].Assign("mu2_2",1,0,p_info);
#ifdef ANALYSE__Phase_Space_Handler
    m_isrzkey[0].Assign("z_1",3,0,p_info);
    m_isrzkey[1].Assign("z_2",3,0,p_info);
#endif
    m_isrkpkey[0].Assign("k_perp_1",4,1,p_info);
    m_isrkpkey[1].Assign("k_perp_2",4,1,p_info);
    p_beamhandler->AssignKeys(p_info);
    p_isrhandler->AssignKeys(p_info);
  }
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p_fsrchannels)  { delete p_fsrchannels;  p_fsrchannels = 0;   }
  if (p_kpchannels)   { delete p_kpchannels;   p_kpchannels = 0;    }
  if (p_zchannels)    { delete p_zchannels;    p_zchannels = 0;     }
  if (p_isrchannels)  { delete p_isrchannels;  p_isrchannels = 0;   }
  if (p_beamchannels) { delete p_beamchannels; p_beamchannels  = 0; }
  if (p_cuts)         { delete p_cuts;         p_cuts = 0;          }
  if (p_pi) delete p_pi;
  delete [] p_cms;
  delete [] p_lab;
  delete [] p_flavours;
  delete p_integrator;
  delete p_channellibnames;
}

void Phase_Space_Handler::InitCuts() 
{
  if (p_cuts!=NULL) delete p_cuts;
  p_cuts = new ATOOLS::Cut_Data();
  p_cuts->Init(m_nin+m_nout,p_flavours);
  if (p_process->Selector()) (p_process->Selector())->BuildCuts(p_cuts);
}

bool Phase_Space_Handler::InitIncoming(const double _mass) 
{
  if (m_nvec==0) {
    m_nvec=p_process->NVector();
    p_cms = new Vec4D[m_nvec];  
    p_lab = new Vec4D[m_nvec];  
  }
  if (!(MakeIncoming(p_lab)) ) {
    msg_Error()<<"Phase_Space_Handler::Integrate : Error !"<<std::endl
	       <<"  Either too little energy for initial state"
	       <<"  ("<<m_E<<" vs "<<m_m[0]+m_m[1]<<") or "<<std::endl
	       <<"  bad number of incoming particles ("<<m_nin<<")."<<std::endl;
    return 0;
  } 
  if (m_nin>1) {
    InitCuts();
    m_smin=ATOOLS::Max(sqr(p_process->ISRThreshold()),p_cuts->Smin());
  }
  m_initialized=1;
  return 1;
}

double Phase_Space_Handler::Integrate() 
{
  if (p_process->Points()>0 && p_process->TotalError()<m_error*p_process->TotalXS()) 
    return p_process->TotalXS()*rpa.Picobarn();
  p_integrator = new Phase_Space_Integrator();
  if (!InitIncoming()) return 0;
  if (rpa.gen.ModelName()==std::string("ADD") && p_isrhandler->On()==0 && p_beamhandler->On()==0) {
    if (rpa.gen.Ecms()>rpa.gen.ScalarConstant(std::string("M_cut"))) {
      msg_Error()<<"Warning in Phase_Space_Handler::Integrate() :"<<std::endl
		 <<"   Use of model ADD at a c.m. energy of "<<rpa.gen.Ecms()<<" GeV,"<<std::endl
		 <<"   but internal string/cut-off scale of model is "
		 <<rpa.gen.ScalarConstant(std::string("M_cut"))<<" GeV."<<std::endl
		 <<"   Return 0 pb as cross section for process "<<p_process->Name()<<std::endl;
      return 0.;
    }
  }
  msg_Debugging()<<"Phase_Space_Handler::Integrate with : "<<std::endl;
  if (m_nin>1) {
    if (p_beamchannels) 
      msg_Debugging()<<"  Beam   : "<<p_beamchannels->Name()<<" ("<<p_beamchannels<<") "
		     <<"  ("<<p_beamchannels->Number()<<","<<p_beamchannels->N()<<")"<<std::endl;
    if (p_isrchannels) 
      msg_Debugging()<<"  ISR    : "<<p_isrchannels->Name()<<" ("<<p_isrchannels<<") "
		     <<"  ("<<p_isrchannels->Number()<<","<<p_isrchannels->N()<<")"<<std::endl;
    if (p_zchannels) 
      msg_Debugging()<<"  KMR Z  : "<<p_zchannels->Name()<<" ("<<p_zchannels<<") "
 		     <<"  ("<<p_zchannels->Number()<<","<<p_zchannels->N()<<")"<<std::endl;
    if (p_kpchannels) 
      msg_Debugging()<<"  KMR kp : "<<p_kpchannels->Name()<<" ("<<p_kpchannels<<") "
 		     <<"  ("<<p_kpchannels->Number()<<","<<p_kpchannels->N()<<")"<<std::endl;
  }
  msg_Debugging()<<"  FSR    : "<<p_fsrchannels->Name()<<" ("<<p_fsrchannels<<") "
		 <<"  ("<<p_fsrchannels->Number()<<","<<p_fsrchannels->N()<<")"<<std::endl;
  if (m_nin==2) return p_integrator->Calculate(this,m_error,m_fin_opt);
  if (m_nin==1) return p_integrator->CalculateDecay(this,sqrt(p_lab[0].Abs2()),m_error);
  return 0.;
}

bool Phase_Space_Handler::MakeIncoming(ATOOLS::Vec4D *const p,const double mass) 
{
  if (m_nin == 1) {
    if (mass<0.) m_E = m_m[0];
    else m_E = mass;  
    m_flux = 1./(2.*m_E);
    m_s = m_E*m_E;
    p[0] = Vec4D(m_E,0.,0.,0.);
    return 1;
  }
  if (m_nin == 2) {
    if (m_isrspkey[3]==0.) m_isrspkey[3] = sqr(ATOOLS::rpa.gen.Ecms());
    double Eprime = sqrt(m_isrspkey[3]);
    if ((m_E<m_m[0]+m_m[1])) return 0;
    double x = 1./2.+(m_m2[0]-m_m2[1])/(2.*m_isrspkey[3]);
    double E1 = x*Eprime;
    double E2 = (1.-x)*Eprime;
    p[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-sqr(m_m[0])));
    p[1] = Vec4D(E2,(-1.)*Vec3D(p[0]));
    m_flux = 1./(2.*sqrt(sqr(m_isrspkey[3]-m_m2[0]-m_m2[1])-4.*m_m2[0]*m_m2[1]));
    if (p_beamhandler->On()==0 && p_isrhandler->On()==0) {
      double eb1=p_beamhandler->GetBeam(0)->Energy();
      double eb2=p_beamhandler->GetBeam(1)->Energy();
      p[0] = Vec4D(eb1,0.,0.,sqrt(sqr(eb1)-sqr(m_m[0])));
      p[1] = Vec4D(eb2,0.0,0.0,-sqrt(sqr(eb2)-sqr(m_m[1])));
      if (!p_massboost) p_massboost = new ATOOLS::Poincare(p[0]+p[1]);
      else *p_massboost=ATOOLS::Poincare(p[0]+p[1]);
      for (int i=0;i<m_nin;++i) p_massboost->Boost(p[i]);
    }
    return 1;
  }
  return 0;
} 

double Phase_Space_Handler::Differential() 
{ 
  return Differential(p_process);
}

double Phase_Space_Handler::Differential(Integrable_Base *const process,
					 const psm::code mode) 
{ 
  PROFILE_HERE;
  if (process->Name()=="BFKL_ME")
    return p_process->Differential(p_cms);
  if (!(mode&psm::pi_call) && p_pi!=NULL) {
    p_active=process;
    p_pi->GeneratePoint(mode);
#ifdef ANALYSE__Phase_Space_Handler
    static PS_Histogram_3D *three=NULL;
    if (three==NULL) {
      const char *name=(process->Name()+"_sprime_y_pt").c_str();
      three = new PS_Histogram_3D(name,name,100,0.0,10.0,
				  100,-10.0,10.0,100,0.0,10.0);
      MYROOT::myroot->AddObject(three,name);
    }
    three->Fill(-log10(m_isrspkey[3]/m_isrspkey[2]),
		m_isrykey[2],-log10(p_lab[2].PPerp2()/m_isrspkey[2]),
		m_psweight==0.0?0.0:m_flux*(m_result_1+m_result_2)*
		p_pi->GenerateWeight());
#endif
    return p_pi->Value()*p_pi->GenerateWeight();
  }
  p_info->ResetAll();
  if (m_nin>1) {
    if (!(mode&psm::no_lim_isr)) p_isrhandler->Reset();
    if (p_beamhandler->On()>0) { 
      p_beamhandler->SetSprimeMin(m_smin);
      p_beamhandler->SetLimits();
      p_beamchannels->GeneratePoint(m_beamspkey,m_beamykey,
				    p_beamhandler->On()); 
      if (!p_beamhandler->MakeBeams(p_lab)) return 0.;
      if (!(mode&psm::no_lim_isr)) 
	p_isrhandler->SetSprimeMax(m_beamspkey[3]*
				   p_isrhandler->Upper1()*
				   p_isrhandler->Upper2());
      p_isrhandler->SetPole(m_beamspkey[3]);
    }
    if (!(mode&psm::no_lim_isr)) p_isrhandler->SetSprimeMin(m_smin);
    if (!(mode&psm::no_dice_isr)) {
      p_isrhandler->SetLimits();
      p_isrhandler->SetMasses(p_process->Selected()->Flavours(),m_nout);
      if (p_isrhandler->On()>0) { 
	if (p_isrhandler->KMROn()) {
	  if (mode&psm::pi_kp) p_kpchannels->
	    GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn(),p_pi);
	  else p_kpchannels->
	    GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn());
	  if (mode&psm::pi_z) p_zchannels->
	    GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn(),p_pi);
	  else p_zchannels->
	    GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn());
	  p_isrhandler->
	    SetSprimeMax(p_isrhandler->SprimeMax()+
			 (m_isrkpkey[0](0)+m_isrkpkey[1](0)).Abs2());
	  if (p_isrhandler->SprimeMax()<p_isrhandler->SprimeMin()) return 0.0;
	}
	if (mode&psm::pi_isr) 
	  p_isrchannels->GeneratePoint(m_isrspkey,m_isrykey, 
				       p_isrhandler->On(),p_pi);
	else p_isrchannels->
	  GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->On());
      }
    }
    if (!p_isrhandler->MakeISR(p_lab,m_nvec,
			       p_process->Selected()->Flavours(),m_nin+m_nout)) {
      if (p_beamchannels) p_beamchannels->NoDice();    
      if (p_isrchannels)  p_isrchannels->NoDice();    
      if (p_zchannels)  p_zchannels->NoDice();    
      if (p_kpchannels)  p_kpchannels->NoDice();    
      p_fsrchannels->NoDice();
      return 0.;
    }
    if (p_beamhandler->On()>0 || p_isrhandler->On()>0) {
      if (p_isrhandler->On()==0) m_isrspkey[3]=m_beamspkey[3];
      process->Selector()->UpdateCuts(m_isrspkey[3],m_beamykey[2]+m_isrykey[2],p_cuts);
    }
    else {
      MakeIncoming(p_lab);
    }
  }
  if (mode&psm::pi_fsr) p_fsrchannels->GeneratePoint(p_lab,p_cuts,p_pi);
  else p_fsrchannels->GeneratePoint(p_lab,p_cuts);
  if (!Check4Momentum(p_lab)) {
    msg_Out()<<"WARNING in Phase_Space_Handler::Differential : Check4Momentum(p) failed"<<std::endl;
    for (int i=0;i<m_nin+m_nout;++i) msg_Events()<<i<<":"<<p_lab[i]
 						 <<" ("<<p_lab[i].Abs2()<<")"<<std::endl;
    return 0.;
  }
  double Q2 = -1.;
  m_result_1 = m_result_2 = 0.;
  for (int i=0;i<m_nvec;++i) p_cms[i]=p_lab[i];
  if (m_nin>1) {
    if (p_isrhandler->On()>0) p_isrhandler->BoostInLab(p_lab,m_nvec);
    if (p_beamhandler->On()>0) p_beamhandler->BoostInLab(p_lab,m_nvec);
    if (p_massboost) for (int i=0;i<m_nvec;++i) 
      p_massboost->BoostBack(p_lab[i]);
  }
  if (p_process->NAddOut()>0) 
    p_process->SetAddMomenta(p_isrhandler->KMRMomenta());
  // First part : flin[0] coming from Beam[0] and flin[1] coming from Beam[1]
  if (process->Trigger(p_lab)) {
    m_result_1 = 1.;
    if (m_nin>1) {
      Q2 = process->CalculateScale(p_lab);
      if (p_isrhandler->KMROn()>0) {
	m_mu2key[0][0] = process->Scale(stp::kp21);
	m_mu2key[1][0] = process->Scale(stp::kp22);
      }
      if (p_isrhandler->On()>0 && !(mode&psm::no_dice_isr)) {
	if (!p_isrhandler->CalculateWeight(Q2)) return 0.0;
	p_isrchannels->GenerateWeight(p_isrhandler->On());
 	m_result_1 *= p_isrchannels->Weight();
	if (p_isrhandler->KMROn()) {
	  p_zchannels->GenerateWeight(p_isrhandler->KMROn());
	  m_result_1 *= p_zchannels->Weight();
 	  p_kpchannels->GenerateWeight(p_isrhandler->KMROn());
  	  m_result_1 *= p_kpchannels->Weight();
	}
      }
      if (p_beamhandler->On()>0) {
	p_beamhandler->CalculateWeight(Q2);
	p_beamchannels->GenerateWeight(p_beamhandler->On());
	m_result_1 *= p_beamchannels->Weight() * p_beamhandler->Weight();
      }
    }
    p_fsrchannels->GenerateWeight(p_cms,p_cuts);
    m_psweight = m_result_1 *= p_fsrchannels->Weight();
    if (m_nin>1) {
      if (p_isrhandler->On()==3) m_result_2 = m_result_1;
      if (p_isrhandler->KMROn()==0) m_result_1 *= process->Differential(p_cms);
      else m_result_1 *= process->Differential(p_lab);
    }
    else m_result_1 *= process->Differential(p_lab);
    if (m_nin>1 && p_isrhandler->On()==3) {
      Rotate(p_cms);
      if (!p_isrhandler->CalculateWeight2(Q2)) return 0.0;
      if (m_result_2 != 0.) m_result_2 *= process->Differential2();
      else m_result_2 = 0.;
    }
  }
  if (m_nin>1 && (p_isrhandler->On()>0 || p_beamhandler->On()>0)) {
    m_psweight*=m_flux=p_isrhandler->Flux();
  }
#ifdef ANALYSE__Phase_Space_Handler
  Analysis_Map::const_iterator ait=s_psspy.find(process);
  if (ait==s_psspy.end()) {
    const char *name=(process->Name()+"_sprime_y").c_str();
    PS_Histogram_2D *psh = 
      new PS_Histogram_2D(name,name,100,0.0,10.0,100,-10.0,10.0);
    ait=s_psspy.insert(std::pair<Integrable_Base*,
		       PS_Histogram_2D*>(process,psh)).first;
    MYROOT::myroot->AddObject(psh,name);
    MYROOT::myroot->SetDrawOption("lego2");
  }
  if (p_pi==NULL) ait->second->Fill(-log10(m_isrspkey[3]/m_isrspkey[2]),
				    m_isrykey[2],m_psweight==0.0?0.0:
				    m_flux*(m_result_1+m_result_2));
  else ait->second->Fill(-log10(m_isrspkey[3]/m_isrspkey[2]),
			 m_isrykey[2],m_psweight==0.0?0.0:
			 m_flux*(m_result_1+m_result_2)*
			 p_pi->GenerateWeight());
  if (p_pi==NULL) {
    static PS_Histogram_3D *three=NULL;
    if (three==NULL) {
      const char *name=(process->Name()+"_sprime_y_pt").c_str();
      three = new PS_Histogram_3D(name,name,100,0.0,10.0,
				  100,-10.0,10.0,100,0.0,10.0);
      MYROOT::myroot->AddObject(three,name);
    }
    three->Fill(-log10(m_isrspkey[3]/m_isrspkey[2]),m_isrykey[2],
		-log10(p_lab[2].PPerp2()/m_isrspkey[2]),
		m_psweight==0.0?0.0:m_flux*(m_result_1+m_result_2));
  }
  if (p_isrhandler->KMROn()>0) {
    Analysis_Map::const_iterator ait=s_psz.find(process);
    if (ait==s_psz.end()) {
      const char *name=(process->Name()+"_z1_z2").c_str();
      PS_Histogram_2D *psh = 
	new PS_Histogram_2D(name,name,100,0.0,10.0,100,10.0,10.0);
      ait=s_psz.insert(std::pair<Integrable_Base*,
		       PS_Histogram_2D*>(process,psh)).first;
      MYROOT::myroot->AddObject(psh,name);
      MYROOT::myroot->SetDrawOption("lego2");
    }
    ait->second->Fill(-log10(m_isrzkey[0][2]*m_isrzkey[1][2]),
	  	      log10(m_isrzkey[0][2]/m_isrzkey[1][2]),
		      m_psweight==0.0?0.0:
		      m_flux*(m_result_1+m_result_2));
    ait=s_pskp.find(process);
    if (ait==s_pskp.end()) {
      const char *name=(process->Name()+"_kp1_kp2").c_str();
      PS_Histogram_2D *psh = 
	new PS_Histogram_2D(name,name,100,0.0,10.0,100,-10.0,10.0);
      ait=s_pskp.insert(std::pair<Integrable_Base*,
			PS_Histogram_2D*>(process,psh)).first;
      MYROOT::myroot->AddObject(psh,name);
      MYROOT::myroot->SetDrawOption("lego2");
    }
    ait->second->Fill(-log10(m_isrkpkey[0][3]/m_isrkpkey[0][2]),
		      log10(m_isrkpkey[1][3]/m_isrkpkey[1][2]),
		      m_psweight==0.0?0.0:
		      m_flux*(m_result_1+m_result_2));
  }
#endif
  return m_flux*(m_result_1+m_result_2);
}

bool Phase_Space_Handler::Check4Momentum(const ATOOLS::Vec4D *p) 
{
  Vec4D pin,pout;
  pin = pout = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<m_nin;i++) pin += p[i];
  for (int i=m_nin;i<m_nin+m_nout;i++) pout += p[i];
  double sin = pin.Abs2(), sout = pout.Abs2();
  if (pin!=pout || !ATOOLS::IsEqual(sin,sout)) {
    msg_Error()<<"Phase_Space_Handler::Check4Momentum(..): "
		       <<"Difference: "<<pin-pout<<std::endl;
    return false;
  }
  return true;
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::SameEvent() 
{
  return OneEvent(-1.,1);
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::SameWeightedEvent() 
{
  return WeightedEvent(1);
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::OneEvent(const double mass,const int mode)
{
  PROFILE_HERE;
  bool use_overflow=true;
  if (m_nin==1) use_overflow=false;
  if ((mass<0) && (!m_initialized)) InitIncoming();
  if ((mass>0) && (m_nin==1)) InitIncoming(mass);
  m_weight=1.;
  double value;
  bool rot(false);
  for (int j=1, i=1;i<m_maxtrials+1;i++) {
    if (mode==0) {
      p_process->DeSelect();
      p_process->SelectOne();
    }
    else {
      if (!(p_process->Selected())) {
	msg_Error()<<" ERROR: in Phase_Space_Handler::OneEvent() "<<std::endl;
	return NULL;
      }
    }
    double max = p_process->Selected()->Max();
    if (p_process->Selected()->Overflow()==0.) {
      p_process->Selected()->RestoreInOrder(); // use last order for overflow events
      p_isrhandler->SetRunMode(1);
      value = Differential(p_process->Selected());
      rot=false;
    }
    else {
      value = p_process->Selected()->Overflow();
      std::vector<double> & xinfo =p_process->Selected()->XInfo();
      m_beamspkey[3] = xinfo[0];
      m_beamykey[2]  = xinfo[1];
      m_isrspkey[3]  = xinfo[2];
      m_isrykey[2]   = xinfo[3];
      p_isrhandler->MakeISR(p_lab,m_nvec,p_process->Selected()->Flavours(),m_nin+m_nout);
      for (int i=0;i<m_nvec;++i) p_lab[i]=p_process->Selected()->Momenta()[i];
      if (value > max) {
	p_process->Selected()->SetOverflow(value-max);
	value=max;
      }
      else {
	p_process->Selected()->SetOverflow(0.);
      }
      m_result_1=value;
      m_result_2=0.;
    }
    if (value == 0.) {
      ++j;
    }
    else {
      double disc = 0.;
      if (value > max) {
	if (!use_overflow) {
	  // don't use overflow
	  msg_Out()<<"WARNING in Phase_Space_Handler::OneEvent :"<<std::endl
		   <<"   Shifted maximum in "<<p_process->Selected()->Name()<<" : "
		    <<p_process->Selected()->Max()<<" -> "<<value<<std::endl;
	  p_process->Selected()->SetMax(value*1.001);
	  p_process->SetMax(0.);
	}
	else {
	  // use overflow
	  p_process->Selected()->SetOverflow(value-max);
	  std::vector<double> & xinfo =p_process->Selected()->XInfo();
	  xinfo[0]=m_beamspkey[3];
	  xinfo[1]=m_beamykey[2];
	  xinfo[2]=m_isrspkey[3];
	  xinfo[3]=m_isrykey[2];
	  value=max;	
	}
      }
      else disc  = max*ATOOLS::ran.Get();
      if (value >= disc) {
	m_sumtrials += i;m_events ++;
	if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	  Rotate(p_lab);
	  if (p_process->NAddOut()>0) {
	    ATOOLS::Vec4D *addvecs=(ATOOLS::Vec4D*)p_process->AddMomenta();
	    Rotate(addvecs,p_process->NAddOut());
	  }
	  p_process->Selected()->SetMomenta(p_lab);
	  p_process->Selected()->SwapInOrder();
	  rot=true;
	}
	else {
	  p_process->Selected()->SetMomenta(p_lab);
	}
	return new Blob_Data<Weight_Info>
	  (Weight_Info(1.0,p_process->EnhanceFactor(),
		       p_process->Selected()->TotalXS(),1,1));
      }
      else j=1;
    }
  }
  m_sumtrials += m_maxtrials;
  msg_Out()<<"WARNING in Phase_Space_Handler::OneEvent() : "
	   <<" too many trials for "<<p_process->Selected()->Name()<<std::endl
	   <<"   Efficiency = "<<double(m_events)/double(m_sumtrials)*100.<<" %."<<std::endl;
  return NULL;
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::WeightedEvent(int mode)
{
  PROFILE_HERE;
  if (!m_initialized) InitIncoming();
  
  double value;
  for (int i=1;i<m_maxtrials+1;i++) {
    if (mode==0) {
      if (!p_process->ReSelect(i)) return NULL;
    }
    else {
      if (!(p_process->Selected())) {
	msg_Error()<<"Phase_Space_Handler::WeightedEvent(): "
		   <<"No process selected."<<std::endl;
	return 0;
      }
    }
    Integrable_Base *selected=p_process->Selected();
    selected->RestoreInOrder();
    p_isrhandler->SetRunMode(1);
    value = Differential(selected,(psm::code)mode);
    if (value > 0.) {
      m_sumtrials+=i;
      ++m_events;
      if (p_process->Selected()->Name()!="BFKL_ME") {
	if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	  Rotate(p_lab);
	  if (p_process->NAddOut()>0) {
	    ATOOLS::Vec4D* addvecs=(ATOOLS::Vec4D*)p_process->AddMomenta();
	    Rotate(addvecs,p_process->NAddOut());
	  }
	  selected->SetMomenta(p_lab);
	  selected->SwapInOrder();
	}
	else {
	  selected->SetMomenta(p_lab);
	}
      }
      m_weight=value;
      m_trials=i;
      return new Blob_Data<Weight_Info>
	(Weight_Info(m_weight,p_process->GMin()>=0.0?
		     selected->ProcWeight()/p_process->Parent()->ProcWeight():1.0,
		     m_weight,m_trials,m_trials));
    }
    // call from amisic
    if ((psm::code)mode&psm::no_lim_isr ||
	(psm::code)mode&psm::no_dice_isr) return NULL;
  }
  msg_Out()<<"WARNING in Phase_Space_Handler::WeightedEvent() : "
	   <<" too many trials for "<<p_process->Selected()->Name()<<std::endl;
  m_weight=0.;
  return NULL;
} 

void Phase_Space_Handler::AddPoint(const double value) 
{ 
  p_process->AddPoint(value); 
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p)
{
  if (m_nin==2) m_isrspkey[3] = sqr(ATOOLS::rpa.gen.Ecms());
  Single_Channel * TestCh = new Rambo(m_nin,m_nout,p_flavours);
  MakeIncoming(p);
  TestCh->GeneratePoint(p,p_cuts);
  delete TestCh;
}

void Phase_Space_Handler::WriteOut(const std::string &pID,const bool force) 
{
  Data_Reader read(" ",";","!","=");
  int ovf(0);
  if (!read.ReadFromFile(ovf,"GENERATE_RESULT_DIRECTORY")) ovf=0;
  msg_Tracking()<<"Write out channels into directory : "<<pID<<std::endl;
  int  mode_dir = 448;
  ATOOLS::MakeDir(pID.c_str(),mode_dir,force|ovf); 
  if (p_beamchannels != 0) p_beamchannels->WriteOut(pID+"/MC_Beam");
  if (p_isrchannels  != 0) p_isrchannels->WriteOut(pID+"/MC_ISR");
  if (p_zchannels != 0) p_zchannels->WriteOut(pID+"/MC_KMR_Z");
  if (p_kpchannels!= 0) p_kpchannels->WriteOut(pID+"/MC_KMR_KP");
  if (p_fsrchannels  != 0) p_fsrchannels->WriteOut(pID+"/MC_FSR");
  std::string help     = (pID+"/Random").c_str();
  if (p_pi!=NULL) {
    ATOOLS::MakeDir((pID+"/PI/").c_str(),mode_dir,force); 
    std::ofstream piinfo((pID+"/PI/Integrators").c_str());
    if (piinfo.good()) {
      piinfo<<p_pi->Key()<<" "<<p_pi->Point().size()<<"\n";
      p_pi->WriteOut(pID+"/PI/");
    }
  }
  ran.WriteOutStatus(help.c_str());
}

bool Phase_Space_Handler::ReadIn(const std::string &pID,const size_t exclude) 
{
  msg_Info()<<"Read in channels from directory : "<<pID<<std::endl;
  bool okay = 1;
  if (p_beamchannels!=NULL && !(exclude&1)) okay = okay && p_beamchannels->ReadIn(pID+"/MC_Beam");
  if (p_isrchannels!=NULL && !(exclude&2)) okay = okay && p_isrchannels->ReadIn(pID+"/MC_ISR");
  if (p_zchannels!=NULL && !(exclude&4)) okay = okay && p_zchannels->ReadIn(pID+"/MC_KMR_Z");
  if (p_kpchannels!=NULL && !(exclude&8)) okay = okay && p_kpchannels->ReadIn(pID+"/MC_KMR_KP");
  if (p_fsrchannels!=NULL && !(exclude&16)) okay = okay && p_fsrchannels->ReadIn(pID+"/MC_FSR");
  std::ifstream piinfo((pID+"/PI/Integrators").c_str());
  if (piinfo.good()) {
    size_t dim;
    std::string key;
    piinfo>>key>>dim;
    p_pi = new PI_Interface(this,key,dim);
    p_pi->SetMode((psm::code)m_use_pi);
    msg_Tracking()<<"Phase_Space_Handler::ReadIn("<<pID
		  <<"/PI/): Read in PI."<<std::endl;
    p_pi->ReadIn(pID+"/PI/");
  }
  if (rpa.gen.RandomSeed()==1234 && !(exclude&32)) {
    std::string filename     = (pID+"/Random").c_str();
    ran.ReadInStatus(filename.c_str(),0);
  }
  return okay;
}

void Phase_Space_Handler::Rotate(ATOOLS::Vec4D *const p,const size_t n)
{
  int nvec=n==0?m_nin+m_nout:n;
  for (int i=0;i<nvec;i++) p[i] = Vec4D(p[i][0],(-1.)*Vec3D(p[i]));
}

bool Phase_Space_Handler::LoadChannelLibraries() 
{
  if (p_channellibnames->size()==0) return 1;
  InitCuts();
  for (std::list<std::string>::iterator it=p_channellibnames->begin();it!=p_channellibnames->end();++it) {
    Single_Channel * sc = SetChannel(m_nin,m_nout,p_flavours,(*it),GetInfo());
    if (sc==0) {
      msg_Error()<<"Phase_Space_Handler:"
			 <<"Channels are not compiled and linked yet."<<std::endl
			 <<"Type 'make install' and run again."<<std::endl;
      return 0;
    }
    else {
      sc->SetName((*it));
      p_fsrchannels->Add(sc);
    }
  }
  return 1;
}

bool Phase_Space_Handler::CreateIntegrators()
{
  /*if (p_cuts == 0) {
    p_cuts = new Cut_Data();
    p_cuts->Init(nin+nout,psflavs);
    if (proc->Selector()) (proc->Selector())->BuildCuts(p_cuts);
    } */
  if (m_nin==1) {
    if (m_nout==2) m_inttype = 0;
    if (m_inttype<4)  m_inttype = 0;
    else m_inttype = 4;
  }
  if(!LoadChannelLibraries()) return 0;
  if (m_nin==2) {
    if (p_beamhandler && p_beamhandler->On()>0) {
      if (!(MakeBeamChannels())) {
	msg_Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<std::endl
		   <<"   did not construct any isr channels !"<<std::endl;
      }
    }
    if (p_isrhandler) {
      if (p_isrhandler->On()>0) {
	if (!(MakeISRChannels())) {
	  msg_Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<std::endl
		     <<"   did not construct any isr channels !"<<std::endl;
	}
      }
      if (p_isrhandler->KMROn()) {
	if (!(MakeKMRChannels())) {
	  msg_Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<std::endl
		     <<"   did not construct any kmr channels !"<<std::endl;
	}
      }
    }
  }
  if (m_nin==2) { 
    if (m_nout==2&&m_inttype==2) m_inttype=6;
    if ((m_inttype<3||m_inttype>20) && (p_fsrchannels!=0)) p_fsrchannels->DropAllChannels();
  }
  if (p_process->Name()!="BFKL_ME") {
    switch (m_inttype) {
    case 0:
      {
	bool kk_fs=false;
	for (int i=0;i<m_nout;i++){
	  if (p_flavours[i+m_nin].IsKK()) kk_fs=true;
	}
	if (kk_fs) {
	  p_fsrchannels->Add(new RamboKK(m_nin,m_nout,p_flavours));
	  break;
	}
      }
      
      if (m_nin==1 && m_nout==2) p_fsrchannels->Add(new Decay2Channel(m_nin,m_nout,p_flavours));
      else p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
      break;
    case 1: 
      p_fsrchannels->Add(new Sarge(m_nin,m_nout));
      break;
    case 2: 
      {
	Permutation pp(m_nin+m_nout-1);
	for (int j=0;j<pp.MaxNumber();j++) {
	  int* pm = pp.Get(j);
	  if (pm[1]==0||pm[m_nin+m_nout-3]==0) 
	    p_fsrchannels->Add(new VHAAG(m_nin,m_nout,j));
 	}
      }
      break;
    case 3: 
      p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
      DropRedundantChannels();
      break;
    case 4:case 5:case 6: 
      DropRedundantChannels();
      break;
    default:
      msg_Error()<<"Wrong phasespace integration switch ! Using RAMBO as default."<<std::endl;
      p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
    }  
  }
  msg_Tracking()<<"Initialized Phase_Space_Integrator (\n\t";
  if (p_beamchannels) msg_Tracking()<<p_beamchannels->Name()<<","<<p_beamchannels->Number()<<";\n\t";
  if (p_isrchannels) msg_Tracking()<<p_isrchannels->Name()<<","<<p_isrchannels->Number()<<";\n\t";
  if (p_zchannels) msg_Tracking()<<p_zchannels->Name()<<","<<p_zchannels->Number()<<";\n\t";
  if (p_kpchannels) msg_Tracking()<<p_kpchannels->Name()<<","<<p_kpchannels->Number()<<";\n\t";
  if (p_fsrchannels) msg_Tracking()<<p_fsrchannels->Name()<<","<<p_fsrchannels->Number()<<")"<<std::endl;
  return 1;
}

void Phase_Space_Handler::DropRedundantChannels()
{
  p_fsrchannels->Reset();
  int number = p_fsrchannels->Number();
  if (number<2) return;
  int *marker = new int[number];  
  for (short int i=0;i<number;i++) marker[i] = 0; 
  /*Vec4D** perm_vec = new Vec4D*[number]; 
  for (short int i=0;i<number;i++) perm_vec[i] = new Vec4D[m_nin+m_nout+1];
  // Create Momenta
  int rannum   = 1 + 2 + 3*(m_nout-2);
  double * rans = new double[rannum];
  for (short int i=0;i<rannum;i++) rans[i] = ran.Get();  
  // Init markers for deletion and results to compare.
  double * res    = new double[number];
  for (short int i=0;i<number;i++) { marker[i] = 0;res[i] = 0.; }
  for (short int i=0;i<number;i++) {
    perm_vec[i][0] = Vec4D(rpa.gen.Ecms()/2.,0.,0.,rpa.gen.Ecms()/2.);
    perm_vec[i][1] = Vec4D(rpa.gen.Ecms()/2.,0.,0.,-rpa.gen.Ecms()/2.); 
    p_fsrchannels->GeneratePoint(i,perm_vec[i],p_process->Cuts(),rans);
    p_fsrchannels->GenerateWeight(i,perm_vec[i],p_process->Cuts());
    res[i] = p_fsrchannels->Weight();
    if (res[i]==0.) marker[i] = 1;
  }
  delete[] rans;*/
  // kick identicals & permutations
  for (short int i=0;i<number;i++) {
    if (marker[i]==0) {
      for (short int j=i+1;j<number;j++) {
	if (marker[j]==0) {
	  //if ( (Compare(perm_vec[i],perm_vec[j])) && 
	  //     (ATOOLS::IsEqual(res[i],res[j])) ) {
	  if (CompareCh(p_fsrchannels->ChID(i),p_fsrchannels->ChID(j))) {
	    marker[j] = 1; 
	  }
	}
      }
    }
  }
  // kick non-resonants
  /*
    int max_reson    = 0;
    Flavour * fl_res = 0;
    int * reson      = new int[number];
    for (short int i=0;i<number;i++) {
    if (marker[i]==0) {
    reson[i]     = fsrchannels->CountResonances(i,fl_res);
    if (reson[i]!=0) {
    //shorten
    int hit    = 0;
    for (short int j=0;j<reson[i];j++) {
    if (sqr(fl_res[j].Mass())>ycut*sqr(rpa.gen.Ecms()) &&
    sqr(fl_res[j].Mass())<sqr(rpa.gen.Ecms())) 
    hit++;
    }
    reson[i] = hit;
    if (reson[i]>max_reson) max_reson = reson[i];
    }
    else reson[i] = -1;
    }
    else reson[i] = -1;
    }
    //Drop them
    for (short int i=0;i<number;i++) {
    if (reson[i]<max_reson && reson[i]!=-1) marker[i] = 1;
    }
    delete [] reson;
  */
  int count = 0;
  for (short int i=0;i<number;i++) {
    if (marker[i]) {
      p_fsrchannels->DropChannel(i-count);
      count++;
    }
  }
  //delete [] res;
  delete [] marker; 
  //for (short int i=0;i<number;i++) delete [] perm_vec[i];
  //delete [] perm_vec; 
}
  
bool Phase_Space_Handler::CompareCh(std::string C1,std::string C2)
{
  int l=Min(C1.length(),C1.length());
  for (int i=0;i<l;i++) {
    if (C1[i]!=C2[i]) return 0;
    if (C1[i]=='Z') return 1;
  }
  return 1;
}


bool Phase_Space_Handler::Compare(const Vec4D *p1,const Vec4D *p2)
{
  if (m_nout==2) {
    for (short int i=0;i<m_nout;i++) { 
      if (p1[m_nin+i] != p2[m_nin+i]) return 0;
    }
    return 1;
  }
  else {
    //Identicals
    for (short int i=0;i<m_nout;i++) {
      if (p1[i+m_nin] != p2[i+m_nin]) return 0;
    }
    return 1;
    //Permutations - not reached right now.
    int * perm = new int[m_nout];
    for (short int i=0;i<m_nout;i++) perm[i] = 0; 
    
    int over = 0;
    int hit,sw1;
    for(;;) {
      sw1 = 1;
      for(short int i=0;i<m_nout;i++) {
	for (short int j=i+1;j<m_nout;j++) 
	  if (perm[i]==perm[j]) {sw1 = 0; break;}
      }    
      if (sw1) {
	hit = 1;
	for (short int i=0;i<m_nout;i++) {
	  if (p1[i+m_nin] != p2[perm[i]+m_nin]) {
	    hit = 0;
	    break;
	  }
	}
	if (hit) return 1;
      }
      for (short int j=m_nout-1;j>=0;j--) {
	if ((perm[j]+1)<m_nout) {
	  perm[j]++;            
	  break;
	}
	else {
	  perm[j] = 0;
	  if (j==0) over = 1;
	}
      }
      if (over) break;
    }
    delete[] perm;
    return 0;
  }
}

bool Phase_Space_Handler::MakeBeamChannels()
{
  if (m_beamparams.size()>0) return CreateBeamChannels();
  Channel_Info ci;
  // default : Beamstrahlung
  if ((p_flavours[0].IsLepton()) && (p_flavours[1].IsLepton())) {
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(p_beamhandler->Exponent(1));
    m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  else {
    // Laser Backscattering spectrum
    //if ((p_flavours[0].IsPhoton()) || (p_flavours[1].IsPhoton())) {
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(1.);
    m_beamparams.push_back(ci);
    ci.parameters.clear();

    ci.type = 3;
    (ci.parameters).push_back(p_beamhandler->Peak());
    (ci.parameters).push_back(p_beamhandler->Exponent(1));
    (ci.parameters).push_back(0.7);
    m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  int    type;
  double mass,width;
  double thmin=0.,thmax=0.;
  for (size_t i=0;i<p_fsrchannels->Number();i++) {
    type=0; 
    mass=width=0.;
    if (p_process) p_fsrchannels->ISRInfo(i,type,mass,width);
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    if (type==2) {
      if (thmax==0.) { thmax=mass; thmin=mass; }
      thmin = ATOOLS::Min(thmin,mass);
      thmax = ATOOLS::Max(thmax,mass);
      continue;
    }
    ci.type = type;
    (ci.parameters).push_back(mass);
    if (type==1) (ci.parameters).push_back(width);
    if (type==2) (ci.parameters).push_back(1.5);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    bool add=true;
    for (size_t j=0;j<m_beamparams.size();j++) if (m_beamparams[j]==ci) add=false;
    if (add) m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  if (thmax>0.) {
    ci.type = 2;
    (ci.parameters).push_back(thmax);
    (ci.parameters).push_back(1.5);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    m_beamparams.push_back(ci);
    if (thmin<thmax) {
      (ci.parameters)[0]=thmin;
      m_beamparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  return CreateBeamChannels();
}

bool Phase_Space_Handler::MakeISRChannels()
{
  if (m_isrparams.size()>0) return CreateISRChannels();
  Channel_Info ci;
  int    type;
  double mass,width;
  double thmin=0.,thmax=0.;
  for (size_t i=0;i<p_fsrchannels->Number();i++) {
    type=0; 
    mass=width=0.;
    p_fsrchannels->ISRInfo(i,type,mass,width);
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    if (type==2) {
      if (thmax==0.) { thmax=mass; thmin=mass; }
      thmin = ATOOLS::Min(thmin,mass);
      thmax = ATOOLS::Max(thmax,mass);
      continue;
    }
    ci.type = type;
    (ci.parameters).push_back(mass);
    if (type==1) (ci.parameters).push_back(width);
    if (type==2) (ci.parameters).push_back(2.);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    bool add=true;
    for (size_t j=0;j<m_isrparams.size();j++) if (m_isrparams[j]==ci) add=false; 
    if (add) m_isrparams.push_back(ci);
    ci.parameters.clear();
  }
  if (thmax>0.) {
    ci.type = 2;
    (ci.parameters).push_back(thmax);
    (ci.parameters).push_back(2.);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    m_isrparams.push_back(ci);
    if (thmin<thmax) {
      (ci.parameters)[0]=thmin;
      m_isrparams.push_back(ci);
    }
    ci.parameters.clear();
  }
  
  if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) {
    if ((p_flavours[0].IsLepton() && p_flavours[1].Strong()) ||
	(p_flavours[1].IsLepton() && p_flavours[0].Strong())) {
      //The DIS case
      ci.type = 0; 
      (ci.parameters).push_back(1.);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
    }
    else {
      // leptons : 1/s'^2 and 1/(s-s')^beta, sharp FW-BW peak
      //     ci.type = 0;
      //     (ci.parameters).push_back(.5);
      //     (ci.parameters).push_back(1.);
      //     m_isrparams.push_back(ci);
      //     ci.parameters.clear();
      //     ci.type = 0;
      //     (ci.parameters).push_back(2.);
      //     (ci.parameters).push_back(1.);
      //     m_isrparams.push_back(ci);
      //     ci.parameters.clear();
      ci.type = 3;
      (ci.parameters).push_back(p_isrhandler->Exponent(1));
      (ci.parameters).push_back(1.00000001);
      (ci.parameters).push_back(1.);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  else {
    // default : 1/s'
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(0.5);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    if (thmax==0.) {
      double yexp=.2;
      ci.type = 0;
      (ci.parameters).push_back(.99);
      (ci.parameters).push_back(yexp);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
      ci.type = 0;
      (ci.parameters).push_back(1.5);
      (ci.parameters).push_back(yexp);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
      ci.type = 0;
      (ci.parameters).push_back(2.);
      (ci.parameters).push_back(yexp);
      m_isrparams.push_back(ci);
      ci.parameters.clear();
    }
  }
  return CreateISRChannels();
}

void Phase_Space_Handler::MakeZChannels(const int type)
{
  Channel_Info ci;
  ci.type=type;
  ci.parameters.push_back(.0625);
  ci.parameters.push_back(.99);
  m_zparams.push_back(ci);
  ci.parameters[0]=.125;
  m_zparams.push_back(ci);
  ci.parameters[0]=.25;
  m_zparams.push_back(ci);
  ci.parameters[0]=.5;
  m_zparams.push_back(ci);
  ci.parameters[0]=.75;
  m_zparams.push_back(ci);
  ci.parameters[0]=.875;
  m_zparams.push_back(ci);
  ci.parameters.clear();
}

bool Phase_Space_Handler::MakeKMRChannels()
{
  if (m_zparams.size()>0) return CreateKMRChannels();
  Channel_Info ci;
  if (!p_flavours[0].IsLepton() && !p_flavours[1].IsLepton()) {
    // z channels
    if (p_flavours[0].IsQuark() && p_flavours[1].IsQuark()) {
      MakeZChannels(1);
    }
    else if (p_flavours[0].IsQuark() && p_flavours[1].IsGluon()) {
      MakeZChannels(1);
      MakeZChannels(2);
    }
    else if (p_flavours[0].IsGluon() && p_flavours[1].IsQuark()) {
      MakeZChannels(1);
      MakeZChannels(3);
    }
    else if ((p_flavours[0].IsGluon() && p_flavours[1].IsGluon()) ||
	     (p_flavours[0].IsJet() && p_flavours[1].IsJet())) {
      MakeZChannels(1);
      MakeZChannels(2);
      MakeZChannels(3);
      MakeZChannels(4);
    }
    // k_\perp channels
    ci.type=0;
    ci.parameters.push_back(1.001);
    m_kpparams.push_back(ci);
    ci.parameters[0]=1.167;
    m_kpparams.push_back(ci);
    ci.parameters[0]=1.333;
    m_kpparams.push_back(ci);
    ci.parameters[0]=1.5;
    m_kpparams.push_back(ci);
    ci.parameters[0]=1.667;
    m_kpparams.push_back(ci);
    ci.parameters[0]=1.833;
    m_kpparams.push_back(ci);
    ci.parameters[0]=2.0;
    m_kpparams.push_back(ci);
    ci.parameters.clear();
  }
  return CreateKMRChannels();
}

bool Phase_Space_Handler::CreateBeamChannels()
{
  if (m_beamparams.size() < 1) return 0;
  int beam = p_beamhandler->On();
  Single_Channel * channel;   
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case 0:
      channel = new Simple_Pole_Central_V(m_beamparams[i].parameters[0]," beam",p_info,beam);
      p_beamchannels->Add(channel); 	
      if (beam==3) {
	channel = new Simple_Pole_Forward_V(m_beamparams[i].parameters[0],
					    m_beamparams[i].parameters[1]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new Simple_Pole_Backward_V(m_beamparams[i].parameters[0],
					     m_beamparams[i].parameters[1]," beam",p_info);
	p_beamchannels->Add(channel);
      }
      break;
    case 1:
      channel = new Resonance_Central_V(m_beamparams[i].parameters[0],
					m_beamparams[i].parameters[1]," beam",p_info,beam);
      p_beamchannels->Add(channel);
      if (beam==3) {
	channel = new Resonance_Uniform_V(m_beamparams[i].parameters[0],
					  m_beamparams[i].parameters[1]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new Resonance_Forward_V(m_beamparams[i].parameters[0],
					  m_beamparams[i].parameters[1],
					  m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new Resonance_Backward_V(m_beamparams[i].parameters[0],
					   m_beamparams[i].parameters[1],
					   m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
      }
      break;
    case 2:
      channel = new Threshold_Central_V(m_beamparams[i].parameters[0],
					m_beamparams[i].parameters[1]," beam",p_info,beam);
      p_beamchannels->Add(channel);
      if (beam==3) {
	channel = new Threshold_Forward_V(m_beamparams[i].parameters[0],
					  m_beamparams[i].parameters[1],
					  m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new Threshold_Backward_V(m_beamparams[i].parameters[0],
					   m_beamparams[i].parameters[1],
					   m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
      }
      break;
    case 3:
      if ((p_flavours[0].IsPhoton()) || (p_flavours[1].IsPhoton())) {
	channel = new LBS_Compton_Peak_Central_V(m_beamparams[i].parameters[1],
						 m_beamparams[i].parameters[0]," beam",p_info,beam);
	p_beamchannels->Add(channel);
	if (beam==3) {
	  channel = new LBS_Compton_Peak_Forward_V(m_beamparams[i].parameters[1],
						   m_beamparams[i].parameters[0],
						   m_beamparams[i].parameters[2]," beam",p_info);
	  p_beamchannels->Add(channel);
	  channel = new LBS_Compton_Peak_Backward_V(m_beamparams[i].parameters[1],
						    m_beamparams[i].parameters[0],
						    m_beamparams[i].parameters[2]," beam",p_info);
	  p_beamchannels->Add(channel);
	}
      }
      break;
    }
  }
  return 1;
}

bool Phase_Space_Handler::CreateISRChannels()
{
  if (m_isrparams.size() < 1) return 0;
  int isr = p_isrhandler->On();
  Single_Channel * channel;   
  if (m_use_pi&psm::pi_isr && isr==3) {
    channel = new Flat_ISR_V(1.0," isr",p_info);
    p_isrchannels->Add(channel);
    return true;
  }
  int length = m_isrparams.size();
  for (int i=0;i<length;i++) {
    switch (m_isrparams[i].type) {
    case 0:
      if (isr==3) {
#ifdef PDF_Channels
	// these channels are slower than the usual isr channels
	// however they perform better for hadron pdf's due to 
	// a more realistic x1 / x2 - mapping
 	channel = new Simple_Pole_PDF_Uniform_V(m_isrparams[i].parameters[0]," isr",p_info);
 	p_isrchannels->Add(channel);
	channel = new Simple_Pole_PDF_Forward_V(m_isrparams[i].parameters[0],
						m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Simple_Pole_PDF_Backward_V(m_isrparams[i].parameters[0],
						 m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
#else
 	channel = new Simple_Pole_Uniform_V(m_isrparams[i].parameters[0]," isr",p_info);
 	p_isrchannels->Add(channel);
	channel = new Simple_Pole_Forward_V(m_isrparams[i].parameters[0],
					    m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Simple_Pole_Backward_V(m_isrparams[i].parameters[0],
					     m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
#endif
      }
      else {
	//The channels used for DIS
	channel = new Simple_Pole_Central_V(m_isrparams[i].parameters[0]," isr",p_info,isr);
	p_isrchannels->Add(channel);
	m_use_pi=m_use_pi|psm::pi_isr;
      }
      break;
    case 1:
      if (isr==3) {
	channel = new Resonance_Uniform_V(m_isrparams[i].parameters[0],
					  m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Resonance_Forward_V(m_isrparams[i].parameters[0],
					  m_isrparams[i].parameters[1],
					  m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Resonance_Backward_V(m_isrparams[i].parameters[0],
					   m_isrparams[i].parameters[1],
					   m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
      }
      else {
	channel = new Resonance_Central_V(m_isrparams[i].parameters[0],
					  m_isrparams[i].parameters[1]," isr",p_info,isr);
      }
      break;
    case 2:
      channel = new Threshold_Central_V(m_isrparams[i].parameters[0],
					m_isrparams[i].parameters[1]," isr",p_info,isr);
      p_isrchannels->Add(channel);
      if (isr==3) {
	channel = new Threshold_Forward_V(m_isrparams[i].parameters[0],
					  m_isrparams[i].parameters[1],
					  m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Threshold_Backward_V(m_isrparams[i].parameters[0],
					   m_isrparams[i].parameters[1],
					   m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
      }
      break;
    case 3:
      channel = new Leading_Log_Central_V(m_isrparams[i].parameters[0],
					  m_isrparams[i].parameters[1]," isr",p_info,isr);
      p_isrchannels->Add(channel);
      if (isr==3) {
	channel = new Leading_Log_Forward_V(m_isrparams[i].parameters[0],
					    m_isrparams[i].parameters[1],
					    m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Leading_Log_Backward_V(m_isrparams[i].parameters[0],
					     m_isrparams[i].parameters[1],
					     m_isrparams[i].parameters[2]," isr",p_info);
	p_isrchannels->Add(channel);
      }
      break;
    }
  }
  return 1;
}

bool Phase_Space_Handler::CreateKMRChannels()
{
  if (m_zparams.size() < 1) return 0;
  Single_Channel *channel=NULL;   
  for (size_t i=0;i<m_zparams.size();++i) {
    switch (m_zparams[i].type) {
    case 0:
    case 1: channel = new Leading_Log_Z_QQ(m_zparams[i].parameters[0],
					   m_zparams[i].parameters[1],"",p_info);
      p_zchannels->Add(channel);
      if (m_zparams[i].type!=0) break;
    case 2: channel = new Leading_Log_Z_QG(m_zparams[i].parameters[0],
					   m_zparams[i].parameters[1],"",p_info);
      p_zchannels->Add(channel);
      if (m_zparams[i].type!=0) break;
    case 3: channel = new Leading_Log_Z_GQ(m_zparams[i].parameters[0],
					   m_zparams[i].parameters[1],"",p_info);
      p_zchannels->Add(channel);
      if (m_zparams[i].type!=0) break;
    case 4: channel = new Leading_Log_Z_GG(m_zparams[i].parameters[0],
					   m_zparams[i].parameters[1],"",p_info);
      p_zchannels->Add(channel);
      if (m_zparams[i].type!=0) break;
    default:
      break;
    }
  }
  for (size_t i=0;i<m_kpparams.size();++i) {
    channel = new LL_KPerp(m_kpparams[i].parameters[0],"",p_info);
    p_kpchannels->Add(channel);
  }
  return 1;
}

void Phase_Space_Handler::ISRChannels(const int i,Channel_Info &ci) const 
{
  if (i<(int)m_isrparams.size()) {
    ci.type       = m_isrparams[i].type;
    ci.parameters = m_isrparams[i].parameters;
    return;
  }
  else {
    msg_Error()<<"Error in Phase_Space_Handler::Isrchannels("<<i<<")"<<std::endl
		       <<"  delimiter out of bounds."<<std::endl;
    abort();
  }
}

void Phase_Space_Handler::BeamChannels(const int i,Channel_Info &ci) const 
{
  if (i<(int)m_beamparams.size()) {
    ci.type       = m_beamparams[i].type;
    ci.parameters = m_beamparams[i].parameters;
    return;
  }
  else {
    msg_Error()<<"Error in Phase_Space_Handler::Beamchannels("<<i<<")"<<std::endl
		       <<"  delimiter out of bounds."<<std::endl;
    abort();
  }
}

Integration_Info *const Phase_Space_Handler::GetInfo() 
{
  if (p_info==NULL) return p_info = new Integration_Info();
  return p_info;
}

void Phase_Space_Handler::DeleteInfo() 
{
  delete p_info;
  p_info=NULL;
}

bool Phase_Space_Handler::CreatePI()
{
  size_t dim=0, sdim[4]={0,0,0,0};
  if ((psm::code)m_use_pi&psm::pi_isr) 
    if (p_isrchannels!=NULL && p_isrchannels->Number()!=0) 
      dim+=sdim[0]=(*p_isrchannels)[0]->Dimension()+1;
  if ((psm::code)m_use_pi&psm::pi_fsr) 
    if (p_fsrchannels!=NULL && p_fsrchannels->Number()!=0) 
      dim+=sdim[1]=(*p_fsrchannels)[0]->Dimension()+1;
  if ((psm::code)m_use_pi&psm::pi_z) 
    if (p_zchannels!=NULL && p_zchannels->Number()!=0) 
      dim+=sdim[2]=(*p_zchannels)[0]->Dimension()+1;
  if ((psm::code)m_use_pi&psm::pi_kp) 
    if (p_kpchannels!=NULL && p_kpchannels->Number()!=0) 
      dim+=sdim[3]=(*p_kpchannels)[0]->Dimension()+1;
  if (dim>25) THROW(critical_error,"Dimension too large for PI.");
  if (dim==0) {
    msg_Error()<<"Phase_Space_Handler::CreatePI(): "
		       <<"Zero dimensional phase space. Abort.";
    return false;
  }
  msg_Info()<<"Phase_Space_Handler::CreatePI(..): "
	    <<"Creating "<<dim<<"-dimensional PI\n";
  msg_Tracking()<<"   '"<<p_process->Name()<<"'\n";
  p_pi = new PI_Interface(this,p_process->Name(),dim);
  p_pi->SetMode((psm::code)m_use_pi);
  (*p_pi)->Initialize();
  if (sdim[0]>0) {
    double sum=0.0;
    std::vector<double> alpha(p_isrchannels->Number()-1);
    for (size_t j=0, i=0;i<p_isrchannels->Number()-1;++i) 
      if ((*p_isrchannels)[i]->Alpha()>0.0)
	alpha[j++]=sum+=(*p_isrchannels)[i]->Alpha();
      else alpha.pop_back();
    (*p_pi)->Reserve(p_isrchannels->Name(),sdim[0],1);
    (*p_pi)->Split(p_isrchannels->Name(),0,alpha);
    p_isrchannels->FixAlpha();
  }
  if (sdim[1]>0) {
    double sum=0.0;
    std::vector<double> alpha(p_fsrchannels->Number()-1);
    for (size_t j=0, i=0;i<p_fsrchannels->Number()-1;++i) 
      if ((*p_fsrchannels)[i]->Alpha()>0.0)
	alpha[j++]=sum+=(*p_fsrchannels)[i]->Alpha();
      else alpha.pop_back();
    (*p_pi)->Reserve(p_fsrchannels->Name(),sdim[1],1);
    (*p_pi)->Split(p_fsrchannels->Name(),0,alpha);
    p_fsrchannels->FixAlpha();
  }
  if (sdim[2]>0) {
    double sum=0.0;
    std::vector<double> alpha(p_zchannels->Number()-1);
    for (size_t j=0, i=0;i<p_zchannels->Number()-1;++i)
      if ((*p_zchannels)[i]->Alpha()>0.0)
	alpha[j++]=sum+=(*p_zchannels)[i]->Alpha();
      else alpha.pop_back();
    (*p_pi)->Reserve(p_zchannels->Name(),sdim[2],1);
    (*p_pi)->Split(p_zchannels->Name(),0,alpha);
    p_zchannels->FixAlpha();
  }
  if (sdim[3]>0) {
    double sum=0.0;
    std::vector<double> alpha(p_kpchannels->Number()-1);
    for (size_t j=0, i=0;i<p_kpchannels->Number()-1;++i)
      if ((*p_kpchannels)[i]->Alpha()>0.0)
	alpha[j++]=sum+=(*p_kpchannels)[i]->Alpha();
      else alpha.pop_back();
    (*p_pi)->Reserve(p_kpchannels->Name(),sdim[3],1);
    (*p_pi)->Split(p_kpchannels->Name(),0,alpha);
    p_kpchannels->FixAlpha();
  }
  return true;
}

typedef Single_Channel * (*Lib_Getter_Function)(int nin,int nout,ATOOLS::Flavour* fl
					    , ATOOLS::Integration_Info * const info,Phase_Space_Handler *psh);

Single_Channel * Phase_Space_Handler::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						 std::string& pID, ATOOLS::Integration_Info * const info)
{
  int pos=pID.find("/");
  std::string libname=ATOOLS::rpa.gen.Variable("SHERPA_LIB_PATH")+
    "/libProc_"+pID.substr(0,pos)+LIB_SUFFIX;
  std::string gettername="Getter_"+pID.substr(pos+1); 

  char * error;
  void * module;
  Lib_Getter_Function GetterFunction;

  // try loading library
  module = dlopen(libname.c_str(),RTLD_LAZY);
  error  = dlerror();
  if (module==NULL) return 0;

  GetterFunction = (Lib_Getter_Function)dlsym(module,gettername.c_str());
  error  = dlerror();
  if (error!=NULL) return 0;

  return GetterFunction(nin,nout,fl,info,this);
}

template Weight_Info &ATOOLS::Blob_Data_Base::Get<Weight_Info>();

namespace ATOOLS {
  std::ostream & operator<<(std::ostream & s, const PHASIC::Weight_Info & wi)
  {
    s<<" weight="<<wi.weight<<"   ntrial="<<wi.ntrial<<std::endl;
    return s;
  }
}

namespace ATOOLS {
template <> Blob_Data<Weight_Info>::~Blob_Data() {}

template class Blob_Data<Weight_Info>;
}
