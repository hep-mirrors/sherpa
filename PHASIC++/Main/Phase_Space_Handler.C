#include "PHASIC++/Main/Phase_Space_Handler.H"

#include "PHASIC++/Main/Phase_Space_Integrator.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Channels/FSR_Channels.H"
#include "PHASIC++/Channels/ISR_Channels.H"
#include "PHASIC++/Channels/Beam_Channels.H"
#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Process/Process_Info.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Message.H"  
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Data_Writer.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Smart_Pointer.C"

using namespace PHASIC;
using namespace ATOOLS;
using namespace BEAM;
using namespace PDF;

Integration_Info *PHASIC::Phase_Space_Handler::p_info=NULL;

namespace ATOOLS { template class SP(Phase_Space_Handler); }

Phase_Space_Handler::Phase_Space_Handler(Process_Integrator *proc,
					 ISR_Handler *ih,Beam_Spectra_Handler *bh, double error): 
  m_name(proc->Process()->Name()), p_process(proc), p_active(proc), p_integrator(NULL), p_cuts(NULL),
  p_eint(NULL), p_beamhandler(bh), p_isrhandler(ih), p_fsrchannels(NULL),
  p_isrchannels(NULL), p_beamchannels(NULL), p_massboost(NULL),
  m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(0), m_dmode(1), m_initialized(0), m_sintegrator(0),
  m_maxtrials(1000000), m_sumtrials(0), m_events(0), m_E(ATOOLS::rpa.gen.Ecms()), m_s(m_E*m_E), 
  m_weight(1.)
{
  Data_Reader dr(" ",";","!","=");
  dr.AddComment("#");
  dr.AddWordSeparator("\t");
  dr.SetInputPath(rpa.GetPath());
  dr.SetInputFile(rpa.gen.Variable("INTEGRATION_DATA_FILE"));
  m_error    = dr.GetValue<double>("ERROR",0.01);
  m_maxtrials = dr.GetValue<int>("MAX_TRIALS",1000000);
  m_fin_opt  = dr.GetValue<std::string>("FINISH_OPTIMIZATION","On")=="On"?1:0;
  if (error>0.) {
    m_error   = error;
  }
  p_flavours=proc->Process()->Flavours();
  p_fsrchannels = new FSR_Channels(this,"fsr_"+proc->Process()->Name());
  m_m[0] = p_flavours[0].Mass(); m_m2[0] = m_m[0]*m_m[0];
  if (m_nin==2) {
    m_m[1] = p_flavours[1].Mass(); m_m2[1] = m_m[1]*m_m[1]; 
    if (p_beamhandler) {
      if (p_beamhandler->On()>0) p_beamchannels = new Beam_Channels(this,"beam_"+proc->Process()->Name());
    }
    if (p_isrhandler && p_isrhandler->On()>0) {
      p_isrchannels = new ISR_Channels(this,"isr_"+proc->Process()->Name());
    }
  }
  if (m_nin==2) {
    m_isrspkey.Assign("s' isr",4,0,p_info);
    m_isrykey.Assign("y isr",3,0,p_info);
    m_isrxkey.Assign("x isr",5,0,p_info);
    m_beamspkey.Assign("s' beam",4,0,p_info);
    m_beamykey.Assign("y beam",3,0,p_info);
    p_beamhandler->AssignKeys(p_info);
  }
#ifdef USING__Threading
  m_uset=0;
#endif
  m_nvec=m_nin+m_nout;
  p_lab.resize(m_nvec);  
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p_fsrchannels) delete p_fsrchannels;
  if (p_isrchannels) delete p_isrchannels;
  if (p_beamchannels) delete p_beamchannels;
  if (p_cuts) delete p_cuts;
  if (p_eint) delete p_eint;
  if (p_massboost) delete p_massboost;
  delete p_integrator;
}

void Phase_Space_Handler::InitCuts() 
{
  if (p_cuts!=NULL) delete p_cuts;
  p_cuts = new Cut_Data();
  p_process->Process()->InitCuts(p_cuts);
  p_process->Process()->FillOnshellConditions();
  p_process->Process()->BuildCuts(p_cuts);
}

bool Phase_Space_Handler::InitIncoming() 
{
  if (!(MakeIncoming(&p_lab.front())) ) {
    msg_Error()<<"Phase_Space_Handler::Integrate : Error !"<<std::endl
	       <<"  Either too little energy for initial state"
	       <<"  ("<<m_E<<" vs "<<m_m[0]+m_m[1]<<") or "<<std::endl
	       <<"  bad number of incoming particles ("<<m_nin<<")."<<std::endl;
    return 0;
  } 
  if (m_nin>1) {
    m_smin=ATOOLS::Max(sqr(p_process->ISRThreshold()),p_cuts->Smin());
  }
  m_initialized=1;
  return 1;
}

double Phase_Space_Handler::Integrate() 
{
  if (p_process->Points()>0 && p_process->TotalError()<dabs(m_error*p_process->TotalXS())) 
    return p_process->TotalXS()*rpa.Picobarn();
  p_integrator = new Phase_Space_Integrator();
  if (!InitIncoming()) return 0;
  if (MODEL::s_model->Name()==std::string("ADD") && p_isrhandler->On()==0 && p_beamhandler->On()==0) {
    if (rpa.gen.Ecms()>MODEL::s_model->ScalarConstant(std::string("M_cut"))) {
      msg_Error()<<"Warning in Phase_Space_Handler::Integrate() :"<<std::endl
		 <<"   Use of model ADD at a c.m. energy of "<<rpa.gen.Ecms()<<" GeV,"<<std::endl
		 <<"   but internal string/cut-off scale of model is "
		 <<MODEL::s_model->ScalarConstant(std::string("M_cut"))<<" GeV."<<std::endl
		 <<"   Return 0 pb as cross section for process "<<p_process->Process()->Name()<<std::endl;
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
  if (m_nout>3 && (p_process->Process()->ThreadInfo()&1)) {
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
  }
#endif
  if (p_beamchannels) p_beamchannels->Print();
  if (p_isrchannels) p_isrchannels->Print();
  p_fsrchannels->Print();
  SetUpEnhance();
  m_dmode=0;
  double res(0.0);
  if (m_nin==2) res=p_integrator->Calculate(this,m_error,m_fin_opt);
  if (m_nin==1) res=p_integrator->CalculateDecay(this,m_error);
  m_dmode=1;
#ifdef USING__Threading
  if (m_uset) {
  m_uset=0;
  m_sig=0;
  int tec(0);
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
  }
#endif
  return res;
}

bool Phase_Space_Handler::MakeIncoming(ATOOLS::Vec4D *const p) 
{
  if (m_nin == 1) {
    m_E = m_m[0];
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

void Phase_Space_Handler::SetUpEnhance()
{
  std::string func(p_process->EnhanceFunction());
  if (func!="1" && p_eint==NULL) {
    p_eint = new Algebra_Interpreter();
    p_eint->SetTagReplacer(this);
    for (int i(0);i<m_nin+m_nout;++i) 
      p_eint->AddTag("p["+ToString(i)+"]",ToString(p_lab[i]));
    p_eint->Interprete(func);
  }
}

std::string Phase_Space_Handler::ReplaceTags(std::string &expr) const
{
  return p_eint->ReplaceTags(expr);
}

Term *Phase_Space_Handler::ReplaceTags(Term *term) const
{
  term->Set(p_lab[term->Id()]);
  return term;
}

void Phase_Space_Handler::AssignId(Term *term)
{
  term->SetId(ToType<int>
	      (term->Tag().substr
	       (2,term->Tag().length()-3)));
}

double Phase_Space_Handler::EnhanceFactor()
{
  if (p_eint==NULL) return 1.0;
  return p_eint->Calculate()->Get<double>();
}

double Phase_Space_Handler::Differential() 
{ 
  return Differential(p_process);
}

void Phase_Space_Handler::CalculateME()
{
  if (m_nin==1) {
    m_result_1=p_active->Process()->Differential(p_lab);
  }
  else {
    m_result_1=p_active->Process()->Differential(p_lab);
    if (p_isrhandler->On()!=3) {
      m_result_2=0.0;
    }
    else {
      m_result_2=p_active->Process()->Differential2();
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
  p_fsrchannels->GenerateWeight(&p_lab.front(),p_cuts);
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

double Phase_Space_Handler::Differential(Process_Integrator *const process,
					 const psm::code mode) 
{ 
  m_cmode=mode;
  p_active=process;
  if (!process->Process()->GeneratePoint()) return 0.0;
  p_info->ResetAll();
  if (m_nin>1) {
    if (!(mode&psm::no_lim_isr)) p_isrhandler->Reset();
    if (p_beamhandler->On()>0) { 
      p_beamhandler->SetSprimeMin(m_smin);
      p_beamhandler->SetLimits();
      p_beamchannels->GeneratePoint(m_beamspkey,m_beamykey,
				    p_beamhandler->On()); 
      if (!p_beamhandler->MakeBeams(&p_lab.front())) return 0.;
      if (!(mode&psm::no_lim_isr)) 
	p_isrhandler->SetSprimeMax(m_beamspkey[3]*
				   p_isrhandler->Upper1()*
				   p_isrhandler->Upper2());
      p_isrhandler->SetPole(m_beamspkey[3]);
    }
    if (!(mode&psm::no_lim_isr)) p_isrhandler->SetSprimeMin(m_smin);
    if (!(mode&psm::no_dice_isr)) {
      p_isrhandler->SetLimits(m_isrspkey.Doubles(),m_isrykey.Doubles(),
			      m_isrxkey.Doubles());
      p_isrhandler->SetMasses(process->Process()->Selected()->Flavours());
      if (p_isrhandler->On()>0) { 
	p_isrchannels->GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->On());
      }
    }
    if (!p_isrhandler->MakeISR(m_isrspkey[3],m_isrykey[2],
			     p_lab,process->Process()->
			     Selected()->Flavours())) {
      if (p_beamchannels) p_beamchannels->NoDice();    
      if (p_isrchannels)  p_isrchannels->NoDice();    
      p_fsrchannels->NoDice();
      return 0.;
    }
    if (p_beamhandler->On()>0 || p_isrhandler->On()>0) {
      if (p_isrhandler->On()==0) m_isrspkey[3]=m_beamspkey[3];
      p_cuts->Update(m_isrspkey[3],m_beamykey[2]+m_isrykey[2]);
    }
    else {
      MakeIncoming(&p_lab.front());
    }
  }
  if (m_nin>1) {
    if (p_isrhandler->On()>0) p_isrhandler->BoostInLab(&p_lab.front(),m_nvec);
    if (p_beamhandler->On()>0) p_beamhandler->BoostInLab(&p_lab.front(),m_nvec);
    if (p_massboost) for (int i=0;i<m_nvec;++i) 
      p_massboost->BoostBack(p_lab[i]);
  }
  p_fsrchannels->GeneratePoint(&p_lab.front(),p_cuts);
  if (!Check4Momentum(p_lab)) {
    msg_Out()<<"WARNING in Phase_Space_Handler::Differential : Check4Momentum(p) failed"<<std::endl;
    for (int i=0;i<m_nin+m_nout;++i) msg_Events()<<i<<":"<<p_lab[i]
 						 <<" ("<<p_lab[i].Abs2()<<")"<<std::endl;
    return 0.;
  }
  m_result_2=m_result_1=0.0;
  if (process->Process()->Trigger(p_lab)) {
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
      CalculatePS();
      CalculateME();
      if (m_result_1+m_result_2==0.) { return 0.;}
    }
#else
    CalculatePS();
    CalculateME();
    if (m_result_1+m_result_2==0.) { return 0.;}
#endif
    msg_Debugging()<<"csum: me = "<<m_result_1<<" / "<<m_result_2
		   <<", ps = "<<m_psweight<<", p[2] = "<<p_lab[2]
		   <<" "<<p_active->Process()->Name()<<"\n";
    m_result_1*=m_psweight;
    m_result_2*=m_psweight;
  }
  NLO_subevtlist* nlos=p_active->Process()->GetSubevtList();
  if (nlos) {
    for (size_t i=0;i<nlos->size();i++) {
      (*nlos)[i]->m_flip=0;
    }
    (*nlos)*=m_psweight;
    (*nlos).MultMEwgt(m_psweight);
  }
  return (m_result_1+m_result_2);
}

bool Phase_Space_Handler::Check4Momentum(const ATOOLS::Vec4D_Vector &p) 
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

Weight_Info *Phase_Space_Handler::OneEvent(Process_Base *const proc)
{
  if (!m_initialized) InitIncoming();
  if (proc==NULL) THROW(fatal_error,"No process.");
  Process_Integrator *cur(proc->Integrator());
  m_weight=1.;
  for (int j=1, i=1;i<m_maxtrials+1;i++) {
    cur->RestoreInOrder();
    p_isrhandler->SetRunMode(1);
    double value = Differential(cur), max = cur->Max();
    if (value == 0.) {
      ++j;
    }
    else {
      double disc = 0.;
      if (value > max) {
	msg_Info()<<METHOD<<"(): Point for '"<<cur->Process()->Name()
		  <<"' exceeds maximum by "<<value/cur->Max()-1.0<<"."<<std::endl;
      }
      else disc  = max*ATOOLS::ran.Get();
      if (value >= disc) {
	m_sumtrials += i;m_events ++;
        double xf1(0.0), xf2(0.0), mu12(0.0), mu22(0.0), dxs(0.0);
	if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	  cur->SetMomenta(p_lab);
	  cur->SwapInOrder();
	  dxs=m_result_2/m_psweight;
          xf1=p_isrhandler->XF1(1);
          xf2=p_isrhandler->XF2(1);
          mu12=p_isrhandler->MuF2(1);
          mu22=p_isrhandler->MuF2(0);
	}
	else {
	  cur->SetMomenta(p_lab);
	  dxs=m_result_1/m_psweight;
          xf1=p_isrhandler->XF1(0);
          xf2=p_isrhandler->XF2(0);
          mu12=p_isrhandler->MuF2(0);
          mu22=p_isrhandler->MuF2(1);
	}
	return new Weight_Info(1,value,dxs,1.0,xf1,xf2,mu12,mu22);
      }
      else j=1;
    }
  }
  m_sumtrials += m_maxtrials;
  msg_Out()<<"WARNING in Phase_Space_Handler::OneEvent() : "
	   <<" too many trials for "<<proc->Selected()->Name()<<std::endl
	   <<"   Efficiency = "<<double(m_events)/double(m_sumtrials)*100.<<" %."<<std::endl;
  return NULL;
}

Weight_Info *Phase_Space_Handler::WeightedEvent(Process_Base *const proc,int mode)
{
  if (!m_initialized) InitIncoming();
  if (proc==NULL) THROW(fatal_error,"No process.");
  Process_Integrator *selected(proc->Integrator());
  double value;
  selected->RestoreInOrder();
  p_isrhandler->SetRunMode(1);
  
  value = Differential(selected,(psm::code)mode);
  if (value != 0.) {
    ++m_events;
    double xf1(0.0), xf2(0.0), mu12(0.0), mu22(0.0), dxs(0.0);
    ME_wgtinfo* wgtinfo=p_active->Process()->GetMEwgtinfo();
    double pnf=dabs(m_result_1)/dabs(m_result_1+m_result_2);
    if (pnf < ATOOLS::ran.Get()) {
      selected->SetMomenta(p_lab);
      selected->SwapInOrder();
      dxs=m_result_2/m_psweight;
      xf1=p_isrhandler->XF1(1);
      xf2=p_isrhandler->XF2(1);
      mu12=p_isrhandler->MuF2(1);
      mu22=p_isrhandler->MuF2(0);
      NLO_subevtlist* nlos=p_active->Process()->GetSubevtList();
      if (nlos) {
	(*nlos).MultMEwgt(1./(1.-pnf));
	for (size_t l=0;l<nlos->size();l++) 
	  (*nlos)[l]->m_flip=1;
      }	
      if (wgtinfo) {
	(*wgtinfo)*=m_psweight/(1.-pnf);
	wgtinfo->m_x1=p_isrhandler->X1();
	wgtinfo->m_x2=p_isrhandler->X2();
	wgtinfo->Flip();
      }
    }
    else {
      selected->SetMomenta(p_lab);
      dxs=m_result_1/m_psweight;
      xf1=p_isrhandler->XF1(0);
      xf2=p_isrhandler->XF2(0);
      mu12=p_isrhandler->MuF2(0);
      mu22=p_isrhandler->MuF2(1);
      NLO_subevtlist* nlos=p_active->Process()->GetSubevtList();
      if (nlos) (*nlos).MultMEwgt(1./pnf);
      if (wgtinfo) {
	(*wgtinfo)*=m_psweight/pnf;
	wgtinfo->m_x1=p_isrhandler->X1();
	wgtinfo->m_x2=p_isrhandler->X2();
      }
    }
    m_weight=value;
    return new Weight_Info(2,value,dxs,1.0,xf1,xf2,mu12,mu22);
  }

  return NULL;
} 

void Phase_Space_Handler::AddPoint(const double value) 
{ 
  p_process->AddPoint(value); 
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    ATOOLS::Vec4D_Vector cp,ATOOLS::Flavour_Vector fl,
				    const Subprocess_Info *info,size_t &n)
{
  size_t nin(fl.size());
  for (size_t i(0);i<nin;++i) msg_Debugging()<<fl[i]<<" ";
  msg_Debugging()<<"->";
  fl.resize(nin+info->m_ps.size());
  cp.resize(nin+info->m_ps.size());
  for (size_t i(0);i<info->m_ps.size();++i) {
    fl[nin+i]=info->m_ps[i].m_fl;
    msg_Debugging()<<" "<<fl[nin+i];
  }
  msg_Debugging()<<" {\n";
  if (info->m_ps.size()==1) {
    for (size_t i(0);i<nin;++i) cp.back()+=cp[i];
  }
  else {
    Single_Channel * TestCh = new Rambo(nin,info->m_ps.size(),&fl.front());
    TestCh->GeneratePoint(&cp.front(),(Cut_Data*)(NULL));
    delete TestCh;
    if (nin==1) {
      Poincare cms(cp.front());
      for (size_t i(1);i<cp.size();++i) cms.BoostBack(cp[i]);
    }
  }
  for (size_t i(0);i<info->m_ps.size();++i) {
    msg_Indent();
    if (info->m_ps[i].m_ps.empty()) {
      msg_Debugging()<<"p["<<n<<"] = "<<cp[nin+i]<<", m = "
		     <<sqrt(dabs(cp[nin+i].Abs2()))<<" ("<<fl[nin+i]<<")\n";
      p[n++]=cp[nin+i];
    }
    else {
      msg_Debugging()<<"P["<<nin+i<<"] = "<<cp[nin+i]<<", m = "
		     <<sqrt(dabs(cp[nin+i].Abs2()))<<" ("<<fl[nin+i]<<")\n";
      Vec4D_Vector ncp(1,cp[nin+i]);
      Flavour_Vector nfl(1,info->m_ps[i].m_fl);
      TestPoint(p,ncp,nfl,&info->m_ps[i],n);
    }
  }
  msg_Debugging()<<"}\n";
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    const Process_Info *info)
{
  DEBUG_FUNC("");
  Flavour_Vector fl_i(info->m_ii.GetExternal());
  Vec4D_Vector cp(fl_i.size());
  if (fl_i.size()==1) {
    p[0]=cp[0]=Vec4D(fl_i[0].Mass(),0.0,0.0,0.0);
    msg_Debugging()<<"p[0] = "<<p[0]<<"\n";
  }
  else {
    double m[2]={fl_i[0].Mass(),fl_i[1].Mass()};
    double E=rpa.gen.Ecms();
    if (info->m_fi.m_ps.size()==1)
      E=info->m_fi.m_ps.front().m_fl.Mass();
    if (E<m[0]+m[1]) return;
    double x=1.0/2.0+(m[0]*m[0]-m[1]*m[1])/(2.0*E*E);
    p[0]=cp[0]=Vec4D(x*E,0.0,0.0,sqrt(sqr(x*E)-m[0]*m[0]));
    p[1]=cp[1]=Vec4D((1.0-x)*E,Vec3D(-p[0]));
    msg_Debugging()<<"p[0] = "<<p[0]<<"\np[1] = "<<p[1]<<"\n";
  }
  
  unsigned int osd_counter=0;
  for (size_t i=0;i<info->m_fi.GetDecayInfos().size();i++)
    if (info->m_fi.GetDecayInfos()[i].m_osd) osd_counter++;  
    
  if (osd_counter==info->m_fi.GetDecayInfos().size()) {
    size_t n(fl_i.size());
    TestPoint(p,cp,fl_i,&info->m_fi,n);
  }
  else {
    Flavour_Vector fl_f(info->m_fi.GetExternal());
    Flavour_Vector fl_tot(fl_i);
    fl_tot.insert(fl_tot.end(),fl_f.begin(),fl_f.end());
    //
    Single_Channel * TestCh = new Rambo(fl_i.size(),fl_f.size(),&fl_tot.front());
    TestCh->GeneratePoint(p,(Cut_Data*)(NULL));
    //
    delete TestCh;
  }
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p,
				    const size_t &nin,const size_t &nout,
				    const Flavour_Vector &flavs)
{
  if (nin==1) {
    p[0]=Vec4D(flavs[0].Mass(),0.0,0.0,0.0);
  }
  else {
    double m[2]={flavs[0].Mass(),flavs[1].Mass()};
    double E=0.5*rpa.gen.Ecms();
    if (E<m[0]+m[1]) return;
    double x=1.0/2.0+(m[0]*m[0]-m[1]*m[1])/(2.0*E*E);
    p[0]=Vec4D(x*E,0.0,0.0,sqrt(sqr(x*E)-m[0]*m[0]));
    p[1]=Vec4D((1.0-x)*E,Vec3D(-p[0]));
  }
  Single_Channel * TestCh = new Rambo(nin,nout,&flavs.front());
  TestCh->GeneratePoint(p,(Cut_Data*)(NULL));
  delete TestCh;
}

void Phase_Space_Handler::WriteOut(const std::string &pID) 
{
  if (p_beamchannels != 0) p_beamchannels->WriteOut(pID+"/MC_Beam");
  if (p_isrchannels  != 0) p_isrchannels->WriteOut(pID+"/MC_ISR");
  if (p_fsrchannels  != 0) p_fsrchannels->WriteOut(pID+"/MC_FSR");
  std::string help     = (pID+"/Random").c_str();
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
  if (rpa.gen.RandomSeed()==1234 && !(exclude&32)) {
    std::string filename     = (pID+"/Random").c_str();
    ran.ReadInStatus(filename.c_str());
  }
  Data_Reader reader;
  reader.SetAddCommandLine(false);
  reader.SetInputPath(pID+"/");
  reader.SetInputFile("Statistics.dat");
  std::vector<std::vector<double> > stats;
  if (reader.MatrixFromFile(stats,"")) m_stats=stats;
  return okay;
}

bool Phase_Space_Handler::CreateIntegrators()
{
  m_sintegrator=p_fsrchannels->Initialize();
  if (m_nin==2) {
    if (p_beamhandler && p_beamhandler->On()>0) {
      if (!p_beamchannels->Initialize()) return false;
    }
    if (p_isrhandler && p_isrhandler->On()>0) {
      if (!p_isrchannels->Initialize()) return false;
    }
  }
  msg_Tracking()<<"Initialized Phase_Space_Integrator (\n\t";
  if (p_beamchannels) msg_Tracking()<<p_beamchannels->Name()<<","<<p_beamchannels->Number()<<";\n\t";
  if (p_isrchannels) msg_Tracking()<<p_isrchannels->Name()<<","<<p_isrchannels->Number()<<";\n\t";
  if (p_fsrchannels) msg_Tracking()<<p_fsrchannels->Name()<<","<<p_fsrchannels->Number()<<")"<<std::endl;
  return true;
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
  p_process->Process()->UpdateIntegrator(this);
  return true;
}

Integration_Info* Phase_Space_Handler::GetInfo() 
{
  if (p_info==NULL) return p_info = new Integration_Info();
  return p_info;
}

void Phase_Space_Handler::DeleteInfo() 
{
  delete p_info;
  p_info=NULL;
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
    s<<" w = "<<wi.m_weight<<std::endl;
    return s;
  }
}

namespace ATOOLS {
template <> Blob_Data<Weight_Info>::~Blob_Data() {}

template class Blob_Data<Weight_Info>;
}
