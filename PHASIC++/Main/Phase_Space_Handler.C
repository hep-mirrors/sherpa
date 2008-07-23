#include "Phase_Space_Handler.H"

#include "Phase_Space_Integrator.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Rambo.H"
#include "RamboKK.H"
#include "Sarge.H"
#include "VHAAG.H"
#include "FSR_Channel.H"
#include "ISR_Vegas.H"
#include "Running_AlphaS.H"
#include "Permutation.H"
#include "Color_Integrator.H"
#include "Helicity_Integrator.H"

#include "Library_Loader.H"
#include "Run_Parameter.H"
#include "Blob.H"
#include "Message.H"  
#include "Random.H"
#include "Shell_Tools.H"
#include "MyStrStream.H"
#include "Data_Reader.H"
#include "Data_Writer.H"
#include "Model_Base.H"

#include <dlfcn.h>

#define PTS long unsigned int
#define PT(ARG) (PTS)(ARG)

using namespace PHASIC;
using namespace ATOOLS;
using namespace BEAM;
using namespace PDF;

Integration_Info *PHASIC::Phase_Space_Handler::p_info=NULL;

Phase_Space_Handler::Phase_Space_Handler(Integrable_Base *proc,
					 ISR_Handler *ih,Beam_Spectra_Handler *bh, double error): 
  m_name(proc->Name()), p_process(proc), p_active(proc), p_integrator(NULL), p_cuts(NULL),
  p_beamhandler(bh), p_isrhandler(ih), p_fsrchannels(NULL),
  p_isrchannels(NULL), p_beamchannels(NULL), p_flavours(NULL), p_cms(NULL), p_lab(NULL), p_massboost(NULL),
  m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(0), m_initialized(0), m_sintegrator(0),
  m_maxtrials(1000000), m_sumtrials(0), m_events(0), m_E(ATOOLS::rpa.gen.Ecms()), m_s(m_E*m_E), 
  m_weight(1.), p_colint(NULL), p_helint(NULL)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(rpa.GetPath());
  dr.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  m_error    = dr.GetValue<double>("ERROR",0.01);
  m_inttype  = dr.GetValue<int>("INTEGRATOR",6);
  m_fin_opt  = dr.GetValue<std::string>("FINISH_OPTIMIZATION","On")=="On"?1:0;
  if (error>0.) {
    m_error   = error;
    m_fin_opt = 0;
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
    if (p_isrhandler && p_isrhandler->On()>0) {
      p_isrchannels = new Multi_Channel("isr_"+proc->Name());
    }
  }
  if (m_nin==2) {
    m_isrspkey.Assign("s' isr",4,0,p_info);
    m_isrykey.Assign("y isr",3,0,p_info);
    m_beamspkey.Assign("s' beam",4,0,p_info);
    m_beamykey.Assign("y beam",3,0,p_info);
    m_mu2key[0].Assign("mu2_1",1,0,p_info);
    m_mu2key[1].Assign("mu2_2",1,0,p_info);
    p_beamhandler->AssignKeys(p_info);
    p_isrhandler->AssignKeys(p_info);
  }
#ifdef USING__Threading
  m_uset=0;
#endif
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p_helint!=NULL) delete p_helint;
  if (p_colint!=NULL) delete p_colint;
  if (p_fsrchannels)  { delete p_fsrchannels;  p_fsrchannels = 0;   }
  if (p_isrchannels)  { delete p_isrchannels;  p_isrchannels = 0;   }
  if (p_beamchannels) { delete p_beamchannels; p_beamchannels  = 0; }
  if (p_cuts)         { delete p_cuts;         p_cuts = 0;          }
  if (p_massboost)    { delete p_massboost;    p_massboost = 0;     }
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
  p_process->FillOnshellConditions();
  if (p_process->Selector()) (p_process->Selector())->BuildCuts(p_cuts);
}

bool Phase_Space_Handler::InitIntegrators()
{
  if (p_process->ColorScheme()==cls::sample && p_colint==NULL) {
    p_colint = new Color_Integrator();
  }
  if (p_process->HelicityScheme()==hls::sample && p_helint==NULL) {
    p_helint = new Helicity_Integrator();
  }
  return true;
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
  if (MODEL::s_model->Name()==std::string("ADD") && p_isrhandler->On()==0 && p_beamhandler->On()==0) {
    if (rpa.gen.Ecms()>MODEL::s_model->ScalarConstant(std::string("M_cut"))) {
      msg_Error()<<"Warning in Phase_Space_Handler::Integrate() :"<<std::endl
		 <<"   Use of model ADD at a c.m. energy of "<<rpa.gen.Ecms()<<" GeV,"<<std::endl
		 <<"   but internal string/cut-off scale of model is "
		 <<MODEL::s_model->ScalarConstant(std::string("M_cut"))<<" GeV."<<std::endl
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
  }
  msg_Debugging()<<"  FSR    : "<<p_fsrchannels->Name()<<" ("<<p_fsrchannels<<") "
		 <<"  ("<<p_fsrchannels->Number()<<","<<p_fsrchannels->N()<<")"<<std::endl;
#ifdef USING__Threading
  pthread_cond_init(&m_sme_cnd,NULL);
  pthread_cond_init(&m_tme_cnd,NULL);
  pthread_mutex_init(&m_sme_mtx,NULL);
  pthread_mutex_init(&m_tme_mtx,NULL);
  pthread_mutex_lock(&m_sme_mtx);
  pthread_mutex_lock(&m_tme_mtx);
  pthread_cond_init(&m_sps_cnd,NULL);
  pthread_cond_init(&m_tps_cnd,NULL);
  pthread_mutex_init(&m_sps_mtx,NULL);
  pthread_mutex_init(&m_tps_mtx,NULL);
  pthread_mutex_lock(&m_sps_mtx);
  pthread_mutex_lock(&m_tps_mtx);
  m_uset=1;
  m_sig=1;
  int tec(0);
  if ((tec=pthread_create(&m_met,NULL,&CalculateME,(void*)this))) {
    THROW(fatal_error,"Cannot create matrix element thread");
  }
  if ((tec=pthread_create(&m_pst,NULL,&CalculatePS,(void*)this)))
    THROW(fatal_error,"Cannot create phase space thread");
#endif
  double res(0.0);
  if (m_nin==2) res=p_integrator->Calculate(this,m_error,m_fin_opt);
  if (m_nin==1) res=p_integrator->CalculateDecay(this,sqrt(p_lab[0].Abs2()),m_error);
#ifdef USING__Threading
  m_uset=0;
  m_sig=0;
  // terminate ps calc thread
  pthread_cond_wait(&m_sps_cnd,&m_sps_mtx);
  if ((tec=pthread_join(m_pst,NULL)))
    THROW(fatal_error,"Cannot join phase space thread");
  pthread_mutex_unlock(&m_tps_mtx);
  pthread_mutex_unlock(&m_sps_mtx);
  pthread_mutex_destroy(&m_tps_mtx);
  pthread_mutex_destroy(&m_sps_mtx);
  pthread_cond_destroy(&m_tps_cnd);
  pthread_cond_destroy(&m_sps_cnd);
  // terminate me calc thread
  pthread_cond_wait(&m_sme_cnd,&m_sme_mtx);
  if ((tec=pthread_join(m_met,NULL)))
    THROW(fatal_error,"Cannot join matrix element thread");
  pthread_mutex_unlock(&m_tme_mtx);
  pthread_mutex_unlock(&m_sme_mtx);
  pthread_mutex_destroy(&m_tme_mtx);
  pthread_mutex_destroy(&m_sme_mtx);
  pthread_cond_destroy(&m_tme_cnd);
  pthread_cond_destroy(&m_sme_cnd);
#endif
  return res;
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

void Phase_Space_Handler::CalculateME()
{
  if (m_nin==1) {
    m_result_1=p_active->Differential(p_lab);
  }
  else {
    double Q2(p_active->CalculateScale(p_lab));
    if (Q2<0.0) {
      m_mu2key[0][0]=p_active->Scale(stp::kp21);
      m_mu2key[1][0]=p_active->Scale(stp::kp22);
    }
    if (p_isrhandler->On()>0 && 
	!(m_cmode&psm::no_dice_isr))
      if (!p_isrhandler->CalculateWeight(Q2)) return; 
    m_result_1=p_active->Differential(p_cms);
    if (p_isrhandler->On()!=3) {
      m_result_2=0.0;
    }
    else {
      if (!p_isrhandler->CalculateWeight2(Q2)) return;
      m_result_2=p_active->Differential2();
    }
    if (p_beamhandler->On()>0) {
      p_beamhandler->CalculateWeight(Q2);
      m_result_1*=p_beamhandler->Weight();
      m_result_2*=p_beamhandler->Weight();
    }
  }
}

void Phase_Space_Handler::CalculatePS()
{
  m_psweight=1.0;
  if (m_nin>1) {
    if (p_isrhandler->On()>0 && 
	!(m_cmode&psm::no_dice_isr)) {
      p_isrchannels->GenerateWeight(p_isrhandler->On());
      m_psweight*=p_isrchannels->Weight();
    }
    if (p_beamhandler->On()>0) {
      p_beamchannels->GenerateWeight(p_beamhandler->On());
      m_psweight*=p_beamchannels->Weight();
    }
  }
  p_fsrchannels->GenerateWeight(p_cms,p_cuts);
  m_psweight*=p_fsrchannels->Weight();
}

#ifdef USING__Threading
void *Phase_Space_Handler::CalculateME(void *arg)
{
  Phase_Space_Handler *psh((Phase_Space_Handler*)arg);
  while (true) {
    // wait for psh to signal
    pthread_mutex_lock(&psh->m_sme_mtx);
    pthread_mutex_unlock(&psh->m_sme_mtx);
    pthread_cond_signal(&psh->m_sme_cnd);
    if (psh->m_sig==0) return NULL;
    psh->CalculateME();
    // signal psh to continue
    pthread_cond_wait(&psh->m_tme_cnd,&psh->m_tme_mtx);
  }
  return NULL;
}

void *Phase_Space_Handler::CalculatePS(void *arg)
{
  Phase_Space_Handler *psh((Phase_Space_Handler*)arg);
  while (true) {
    // wait for psh to signal
    pthread_mutex_lock(&psh->m_sps_mtx);
    pthread_mutex_unlock(&psh->m_sps_mtx);
    pthread_cond_signal(&psh->m_sps_cnd);
    if (psh->m_sig==0) return NULL;
    psh->CalculatePS();
    // signal psh to continue
    pthread_cond_wait(&psh->m_tps_cnd,&psh->m_tps_mtx);
  }
  return NULL;
}
#endif

double Phase_Space_Handler::Differential(Integrable_Base *const process,
					 const psm::code mode) 
{ 
  m_cmode=mode;
  p_active=process;
  if (process->Name().find("BFKL")==0)
    return p_process->Differential(p_cms);
  p_info->ResetAll();
  if (process->ColorScheme()==cls::sample &&
      !p_colint->GeneratePoint(true)) return 0.0;
  if (process->HelicityScheme()==hls::sample &&
      !p_helint->GeneratePoint()) return 0.0;
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
	p_isrchannels->GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->On());
      }
    }
    if (!p_isrhandler->MakeISR(p_lab,m_nvec,
			       p_process->Selected()->Flavours(),m_nin+m_nout)) {
      if (p_beamchannels) p_beamchannels->NoDice();    
      if (p_isrchannels)  p_isrchannels->NoDice();    
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
  p_fsrchannels->GeneratePoint(p_lab,p_cuts);
  if (!Check4Momentum(p_lab)) {
    msg_Out()<<"WARNING in Phase_Space_Handler::Differential : Check4Momentum(p) failed"<<std::endl;
    for (int i=0;i<m_nin+m_nout;++i) msg_Events()<<i<<":"<<p_lab[i]
 						 <<" ("<<p_lab[i].Abs2()<<")"<<std::endl;
    return 0.;
  }
  m_result_1 = m_result_2 = 0.;
  for (int i=0;i<m_nvec;++i) p_cms[i]=p_lab[i];
  if (m_nin>1) {
    if (p_isrhandler->On()>0) p_isrhandler->BoostInLab(p_lab,m_nvec);
    if (p_beamhandler->On()>0) p_beamhandler->BoostInLab(p_lab,m_nvec);
    if (p_massboost) for (int i=0;i<m_nvec;++i) 
      p_massboost->BoostBack(p_lab[i]);
  }
  m_result_2=m_result_1=0.0;
  if (process->Trigger(p_lab)) {
#ifdef USING__Threading
    if (m_uset) {
      // start me calc
      pthread_cond_wait(&m_sme_cnd,&m_sme_mtx);
      // start ps calc
      pthread_cond_wait(&m_sps_cnd,&m_sps_mtx);
      // wait for ps calc to finish
      pthread_mutex_lock(&m_tps_mtx);
      pthread_mutex_unlock(&m_tps_mtx);
      pthread_cond_signal(&m_tps_cnd);
      // wait for me calc to finish
      pthread_mutex_lock(&m_tme_mtx);
      pthread_mutex_unlock(&m_tme_mtx);
      pthread_cond_signal(&m_tme_cnd);
    }
    else {
      CalculateME();
      CalculatePS();
    }
#else
    CalculateME();
    CalculatePS();
#endif
    m_result_1*=m_psweight;
    m_result_2*=m_psweight;
    if (m_nin>1 && p_isrhandler->On()==3) Rotate(p_cms);
  }
  if (m_nin>1 && (p_isrhandler->On()>0 || p_beamhandler->On()>0)) {
    m_psweight*=m_flux=p_isrhandler->Flux();
  }
  return m_flux*(m_result_1+m_result_2);
}

bool Phase_Space_Handler::Check4Momentum(const ATOOLS::Vec4D *p) 
{
  Vec4D pin,pout;
  pin = pout = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<m_nin;i++) pin += p[i];
  for (int i=m_nin;i<m_nin+m_nout;i++) pout += p[i];
  double sin = pin.Abs2(), sout = pout.Abs2();
  static double accu(sqrt(Accu()));
  if (!IsEqual(pin,pout,accu) || !IsEqual(sin,sout,accu)) {
    msg_Error()<<"Phase_Space_Handler::Check4Momentum(..): "
		       <<"Difference: "<<pin-pout<<" "<<sin-sout<<std::endl;
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
        double xf1(0.0), xf2(0.0);
	if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	  Rotate(p_lab);
	  p_process->Selected()->SetMomenta(p_lab);
	  p_process->Selected()->SwapInOrder();
	  rot=true;
          xf1=p_isrhandler->XF1(1);
          xf2=p_isrhandler->XF2(1);
	}
	else {
	  p_process->Selected()->SetMomenta(p_lab);
          xf1=p_isrhandler->XF1(0);
          xf2=p_isrhandler->XF2(0);
	}
	return new Blob_Data<Weight_Info>
	  (Weight_Info(1.0,p_process->EnhanceFactor(),
		       p_process->Selected()->TotalXS(),xf1,xf2,1,1));
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
      double xf1(0.0), xf2(0.0);
      if (p_process->Selected()->Name().find("BFKL")!=0) {
	if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	  Rotate(p_lab);
	  selected->SetMomenta(p_lab);
	  selected->SwapInOrder();
          xf1=p_isrhandler->XF1(1);
          xf2=p_isrhandler->XF2(1);
	}
	else {
	  selected->SetMomenta(p_lab);
          xf1=p_isrhandler->XF1(0);
          xf2=p_isrhandler->XF2(0);
	}
      }
      m_weight=value;
      m_trials=i;
      return new Blob_Data<Weight_Info>
	(Weight_Info(m_weight,p_process->GMin()>=0.0?
		     selected->ProcWeight()/p_process->Parent()->ProcWeight():1.0,
		     m_weight,xf1,xf2,m_trials,m_trials));
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
  if (p_process->ColorScheme()==cls::sample)
    while (!p_colint->GeneratePoint()); 
  if (p_process->HelicityScheme()==hls::sample)
    while (!p_helint->GeneratePoint()); 
  delete TestCh;
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,int nin,int nout, Flavour* flav)
{
  if (nin==2) {
    m_isrspkey[3] = sqr(ATOOLS::rpa.gen.Ecms());
    MakeIncoming(p);
  }
  Single_Channel * TestCh = new Rambo(nin,nout,flav);
  TestCh->GeneratePoint(p,p_cuts);
  delete TestCh;
}

void Phase_Space_Handler::WriteOut(const std::string &pID,const bool force) 
{
  Data_Reader read(" ",";","!","=");
  int ovf(0);
  if (!read.ReadFromFile(ovf,"GENERATE_RESULT_DIRECTORY")) ovf=0;
  msg_Tracking()<<"Write out channels into directory : "<<pID<<std::endl;
  ATOOLS::MakeDir(pID,force|ovf); 
  if (p_beamchannels != 0) p_beamchannels->WriteOut(pID+"/MC_Beam");
  if (p_isrchannels  != 0) p_isrchannels->WriteOut(pID+"/MC_ISR");
  if (p_fsrchannels  != 0) p_fsrchannels->WriteOut(pID+"/MC_FSR");
  std::string help     = (pID+"/Random").c_str();
  if (p_helint!=NULL) p_helint->WriteOut(pID);
  ran.WriteOutStatus(help.c_str());
  Data_Writer writer;
  writer.SetOutputPath(pID+"/");
  writer.SetOutputFile("Statistics.dat");
  writer.MatrixToFile(m_stats);
}

bool Phase_Space_Handler::ReadIn(const std::string &pID,const size_t exclude) 
{
  msg_Info()<<"Read in channels from directory : "<<pID<<std::endl;
  bool okay = 1;
  if (p_beamchannels!=NULL && !(exclude&1)) okay = okay && p_beamchannels->ReadIn(pID+"/MC_Beam");
  if (p_isrchannels!=NULL && !(exclude&2)) okay = okay && p_isrchannels->ReadIn(pID+"/MC_ISR");
  if (p_fsrchannels!=NULL && !(exclude&16)) okay = okay && p_fsrchannels->ReadIn(pID+"/MC_FSR");
  if (p_helint!=NULL) p_helint->ReadIn(pID);
  if (rpa.gen.RandomSeed()==1234 && !(exclude&32)) {
    std::string filename     = (pID+"/Random").c_str();
    ran.ReadInStatus(filename.c_str(),0);
  }
  Data_Reader reader;
  reader.SetAddCommandLine(false);
  reader.SetInputPath(pID+"/");
  reader.SetInputFile("Statistics.dat");
  std::vector<std::vector<double> > stats;
  if (reader.MatrixFromFile(stats,"")) m_stats=stats;
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
    }
  }
  if (m_nin==2) { 
    if (m_nout==2&&m_inttype==2) m_inttype=6;
    if ((m_inttype<3||m_inttype>20) && (p_fsrchannels!=0)) p_fsrchannels->DropAllChannels();
  }
  if (p_process->Name().find("BFKL")!=0) {
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
	VHAAG *firsthaag=NULL,*hlp=NULL;
	Permutation pp(m_nin+m_nout-1);
	for (int j=0;j<pp.MaxNumber();j++) {
	  int* pm = pp.Get(j);
	  if (pm[1]==0||pm[m_nin+m_nout-3]==0) 
	    p_fsrchannels->Add(hlp=new VHAAG(m_nin,m_nout,j,firsthaag));
	  if (!firsthaag) firsthaag=hlp;
 	}
      }
      break;
    case 3: 
      p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
      DropRedundantChannels();
      break;
    case 4:case 5:case 6: 
      DropRedundantChannels();
      m_sintegrator=p_process->FillSIntegrator(p_fsrchannels);
      break;
    default:
      msg_Error()<<"Wrong phasespace integration switch ! Using RAMBO as default."<<std::endl;
      p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
    }  
  }
  msg_Tracking()<<"Initialized Phase_Space_Integrator (\n\t";
  if (p_beamchannels) msg_Tracking()<<p_beamchannels->Name()<<","<<p_beamchannels->Number()<<";\n\t";
  if (p_isrchannels) msg_Tracking()<<p_isrchannels->Name()<<","<<p_isrchannels->Number()<<";\n\t";
  if (p_fsrchannels) msg_Tracking()<<p_fsrchannels->Name()<<","<<p_fsrchannels->Number()<<")"<<std::endl;
  return 1;
}

bool Phase_Space_Handler::UpdateIntegrators()
{
  if (!m_sintegrator) return false;
  double error=Process()->TotalVar()/Process()->TotalResult();
  msg_Info()<<om::blue
	    <<Process()->TotalResult()*rpa.Picobarn()
	    <<" pb"<<om::reset<<" +- ( "<<om::red
	    <<Process()->TotalVar()*rpa.Picobarn()
	    <<" pb = "<<error*100<<" %"<<om::reset<<" ) "
	    <<FSRIntegrator()->ValidN()<<" ( "
	    <<(FSRIntegrator()->ValidN()*1000/FSRIntegrator()->N())/10.0<<" % ) "<<std::endl;
  p_process->UpdateIntegrator(p_fsrchannels);
  return true;
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
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(0.99);
    m_beamparams.push_back(ci);
    ci.parameters.clear();
    // Laser Backscattering spectrum
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
  int length = m_isrparams.size();
  for (int i=0;i<length;i++) {
    switch (m_isrparams[i].type) {
    case 0:
      if (isr==3) {
 	channel = new Simple_Pole_Uniform_V(m_isrparams[i].parameters[0]," isr",p_info);
 	p_isrchannels->Add(channel);
	channel = new Simple_Pole_Forward_V(m_isrparams[i].parameters[0],
					    m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
	channel = new Simple_Pole_Backward_V(m_isrparams[i].parameters[0],
					     m_isrparams[i].parameters[1]," isr",p_info);
	p_isrchannels->Add(channel);
      }
      else {
	//The channels used for DIS
	channel = new Simple_Pole_Central_V(m_isrparams[i].parameters[0]," isr",p_info,isr);
	p_isrchannels->Add(channel);
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

typedef Single_Channel * (*Lib_Getter_Function)(int nin,int nout,ATOOLS::Flavour* fl
					    , ATOOLS::Integration_Info * const info,Phase_Space_Handler *psh);

Single_Channel * Phase_Space_Handler::SetChannel(int nin,int nout,ATOOLS::Flavour* fl,
						 std::string& pID, ATOOLS::Integration_Info * const info)
{
  size_t pos(pID.find("/"));
  s_loader->AddPath(rpa.gen.Variable("SHERPA_LIB_PATH"));
  Lib_Getter_Function gf = (Lib_Getter_Function)
    PT(s_loader->GetLibraryFunction("Proc_"+pID.substr(0,pos),
				    "Getter_"+pID.substr(pos+1)));
  if (gf==NULL) return NULL;
  return gf(nin,nout,fl,info,this);
}

void Phase_Space_Handler::AddStats(const std::vector<double> &stats)
{ 
  std::vector<double> nstats(1,m_stats.size()+1);
  nstats.insert(nstats.end(),stats.begin(),stats.end());
  m_stats.push_back(nstats); 
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
