#include "Phase_Space_Handler.H"

#include "Phase_Space_Integrator.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Rambo.H"
#include "RamboKK.H"
#include "Sarge.H"
#include "Simple_Pole.H"
#include "Resonance.H"
#include "Leading_Log.H"
#include "Threshold.H"
#include "LBS_Compton_Peak.H"
#include "Leading_Log_Z.H"
#include "LL_KPerp.H"
#include "LDL_KPerp.H"
#include "FSR_Channel.H"
#include "Running_AlphaS.H"

#include "Run_Parameter.H"
#include "Message.H"  
#include "Random.H"

#include <sys/types.h>
#include <sys/stat.h>

#ifdef PROFILE__All
#include "prof.hh"
#else
#include "Blob.H"
#ifdef PROFILE__Phase_Space_Handler
#include "prof.hh"
#else
#define PROFILE_HERE ;
#endif
#endif

using namespace PHASIC;
using namespace ATOOLS;
using namespace BEAM;
using namespace PDF;
using namespace std;

Integration_Info *PHASIC::Phase_Space_Handler::p_info=NULL;

Phase_Space_Handler::Phase_Space_Handler(Integrable_Base *proc,
					 ISR_Handler *ih,Beam_Spectra_Handler *bh): 
  m_name(proc->Name()), p_process(proc), p_integrator(NULL), p_cuts(NULL),
  p_beamhandler(bh), p_isrhandler(ih), p_fsrchannels(NULL), p_zchannels(NULL), p_kpchannels(NULL), 
  p_isrchannels(NULL), p_beamchannels(NULL), p_flavours(NULL), p_cms(NULL), p_lab(NULL), 
  m_nin(proc->NIn()), m_nout(proc->NOut()), m_nvec(0), m_initialized(0),
  m_maxtrials(100000), m_sumtrials(0), m_events(0), m_E(ATOOLS::rpa.gen.Ecms()), m_s(m_E*m_E), 
  m_weight(1.)
{
  Data_Read dr(rpa.GetPath()+string("/Integration.dat"));
  m_error = dr.GetValue<double>("ERROR",0.01);
  m_inttype = dr.GetValue<int>("INTEGRATOR",3);
  p_flavours = new Flavour[m_nin+m_nout];
  for (int i=0;i<m_nin+m_nout;i++) p_flavours[i] = proc->Flavours()[i];
  p_fsrchannels = new Multi_Channel(string("fsr_")+proc->Name());
  m_m[0] = p_flavours[0].Mass(); m_m2[0] = m_m[0]*m_m[0];
  if (m_nin==2) {
    m_m[1] = p_flavours[1].Mass(); m_m2[1] = m_m[1]*m_m[1]; 
    if (p_beamhandler) {
      if (p_beamhandler->On()>0) p_beamchannels = new Multi_Channel(string("beam_")+proc->Name());
    }
    if (p_isrhandler) {
      if (p_isrhandler->On()>0) {
	p_isrchannels = new Multi_Channel(string("isr_")+proc->Name());
 	if (p_isrhandler->KMROn()>0) {
 	  p_zchannels = new Multi_Channel(string("kmr_z_")+proc->Name());
 	  p_kpchannels = new Multi_Channel(string("kmr_kp_")+proc->Name());
	}
      }
    }
  }
  ATOOLS::msg.Tracking()<<"Phase_Space_Handler::Phase_Space_Handler(..): "
			<<"Creating initial mapping keys ...\n";
  m_isrspkey.Assign("s' isr",4,0,p_info);
  m_isrykey.Assign("y isr",3,0,p_info);
  m_beamspkey.Assign("s' beam",4,0,p_info);
  m_beamykey.Assign("y beam",3,0,p_info);
  ATOOLS::msg.Tracking()<<"... done."<<std::endl;
  p_beamhandler->AssignKeys(p_info);
  p_isrhandler->AssignKeys(p_info);
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p_fsrchannels)  { delete p_fsrchannels;  p_fsrchannels = 0;   }
  if (p_kpchannels)   { delete p_kpchannels;   p_kpchannels = 0;    }
  if (p_zchannels)    { delete p_zchannels;    p_zchannels = 0;     }
  if (p_isrchannels)  { delete p_isrchannels;  p_isrchannels = 0;   }
  if (p_beamchannels) { delete p_beamchannels; p_beamchannels  = 0; }
  if (p_cuts)         { delete p_cuts;         p_cuts = 0;          }
  delete [] p_cms;
  delete [] p_lab;
  delete [] p_flavours;
  delete p_integrator;
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
  msg.Debugging()<<"Phase_Space_Handler::Integrate with : "<<endl;
  if (m_nin>1) {
    if (p_beamchannels) 
      msg.Debugging()<<"  Beam   : "<<p_beamchannels->Name()<<" ("<<p_beamchannels<<") "
		     <<"  ("<<p_beamchannels->Number()<<","<<p_beamchannels->N()<<")"<<endl;
    if (p_isrchannels) 
      msg.Debugging()<<"  ISR    : "<<p_isrchannels->Name()<<" ("<<p_isrchannels<<") "
		     <<"  ("<<p_isrchannels->Number()<<","<<p_isrchannels->N()<<")"<<endl;
    if (p_zchannels) 
      msg.Debugging()<<"  KMR Z  : "<<p_zchannels->Name()<<" ("<<p_zchannels<<") "
 		     <<"  ("<<p_zchannels->Number()<<","<<p_zchannels->N()<<")"<<endl;
    if (p_kpchannels) 
      msg.Debugging()<<"  KMR kp : "<<p_kpchannels->Name()<<" ("<<p_kpchannels<<") "
 		     <<"  ("<<p_kpchannels->Number()<<","<<p_kpchannels->N()<<")"<<endl;
  }
  msg.Debugging()<<"  FSR    : "<<p_fsrchannels->Name()<<" ("<<p_fsrchannels<<") "
		 <<"  ("<<p_fsrchannels->Number()<<","<<p_fsrchannels->N()<<")"<<endl;
  if (!(MakeIncoming(p_lab)) ) {
    msg.Error()<<"Phase_Space_Handler::Integrate : Error !"<<endl
	       <<"  Either too little energy for initial state"
	       <<"  ("<<m_E<<" vs "<<m_m[0]+m_m[1]<<") or "<<endl
	       <<"  bad number of incoming particles ("<<m_nin<<")."<<endl;
    return 0;
  } 
  if (m_nin>1) {
    InitCuts();
    if (p_beamhandler->On()>0) {
      p_beamhandler->SetSprimeMin(ATOOLS::Max(sqr(p_process->ISRThreshold()),p_cuts->Smin()));
      p_beamchannels->SetRange(p_beamhandler->SprimeRange(),p_beamhandler->YRange());
      p_beamchannels->GetRange();
    }
    if (p_isrhandler->On()>0) {
      p_isrhandler->SetSprimeMin(ATOOLS::Max(ATOOLS::sqr(p_process->ISRThreshold()),
					     p_cuts->Smin()));
      msg.Debugging()<<"In Phase_Space_Handler::Integrate : "<<p_beamhandler->On()<<":"
		     <<p_isrhandler->On()<<endl
		     <<"   "<<p_isrhandler->SprimeMin()<<" ... "<<p_isrhandler->SprimeMax()
		     <<" ... "<<p_isrhandler->Pole()<<endl
		     <<"  for Threshold = "<<p_process->ISRThreshold()<<"  "
		     <<p_process->Name()<<endl;
      p_isrchannels->SetRange(p_isrhandler->SprimeRange(),p_isrhandler->YRange());
      p_isrchannels->GetRange();
    }
    msg.SetPrecision(6);
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
      msg.Error()<<"Warning in Phase_Space_Handler::Integrate() :"<<endl
		 <<"   Use of model ADD at a c.m. energy of "<<rpa.gen.Ecms()<<" GeV,"<<endl
		 <<"   but internal string/cut-off scale of model is "
		 <<rpa.gen.ScalarConstant(std::string("M_cut"))<<" GeV."<<endl
		 <<"   Return 0 pb as cross section for process "<<p_process->Name()<<endl;
      return 0.;
    }
  }
  if (m_nin==2) return p_integrator->Calculate(this,m_error);
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
    double Eprime = sqrt(m_isrspkey[3]);
    if ((m_E<m_m[0]+m_m[1])) return 0;
    double x = 1./2.+(m_m2[0]-m_m2[1])/(2.*m_isrspkey[3]);
    double E1 = x*Eprime;
    double E2 = (1.-x)*Eprime;
    p[0] = Vec4D(E1,0.,0.,sqrt(sqr(E1)-sqr(m_m[0])));
    p[1] = Vec4D(E2,(-1.)*Vec3D(p[0]));
    m_flux = 1./(2.*sqrt(sqr(m_isrspkey[3]-m_m2[0]-m_m2[1])-4.*m_m2[0]*m_m2[1]));
    return 1;
  }
  return 0;
} 

double Phase_Space_Handler::Differential() 
{ 
  return Differential(p_process);
}

double Phase_Space_Handler::Differential(Integrable_Base *const process) 
{ 
  PROFILE_HERE;
  p_info->ResetAll();
  p_isrhandler->Reset();
  if (m_nin>1) {
    if (p_beamhandler->On()>0) { 
      p_beamhandler->SetLimits();
      p_beamchannels->GeneratePoint(m_beamspkey,m_beamykey,p_beamhandler->On()); 
      if (!p_beamhandler->MakeBeams(p_lab)) return 0.;
      p_isrhandler->SetSprimeMax(m_beamspkey[3]*
				 p_isrhandler->Upper1()*p_isrhandler->Upper2());
      p_isrhandler->SetPole(m_beamspkey[3]);
    }
    p_isrhandler->SetLimits();
    if (p_isrhandler->On()>0) { 
      p_isrchannels->GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->On());
      if (p_isrhandler->KMROn()) {
	p_kpchannels->GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn());
	p_zchannels->GeneratePoint(m_isrspkey,m_isrykey,p_isrhandler->KMROn());
      }
    }
    if (!p_isrhandler->MakeISR(p_lab,m_nvec)) return 0.;
    if (p_beamhandler->On()>0 || p_isrhandler->On()>0) {
      process->Selector()->UpdateCuts(m_isrspkey[3],m_beamykey[2]+m_isrykey[2],p_cuts);
    }
  }
  p_fsrchannels->GeneratePoint(p_lab,p_cuts);
  if (!Check4Momentum(p_lab)) {
    msg.Events()<<"Phase_Space_Handler Check4Momentum(p) failed"<<endl;
    for (int i=0;i<m_nin+m_nout;++i) msg.Events()<<i<<":"<<p_lab[i]
 						 <<" ("<<p_lab[i].Abs2()<<")"<<endl;
    return 0.;
  }
  double KFactor = 1., Q2 = -1.;
  m_result_1 = m_result_2 = 0.;
  for (int i=0;i<m_nvec;++i) p_cms[i]=p_lab[i];
  if (p_isrhandler->On()>0) p_isrhandler->BoostInLab(p_lab,m_nvec);
  if (p_beamhandler->On()>0) p_beamhandler->BoostInLab(p_lab,m_nvec);
  // First part : flin[0] coming from Beam[0] and flin[1] coming from Beam[1]
  bool trigger = 0;
  if ((process->Selector())->Trigger(p_lab)) {
    m_result_1 = 1.;
    if (m_nin>1) {
      trigger = 1;
      Q2 = process->Scale(p_lab);
      if (p_isrhandler->On()>0) {
	p_isrhandler->CalculateWeight(Q2);
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
      KFactor *= process->KFactor(Q2);
    }
    p_fsrchannels->GenerateWeight(p_cms,p_cuts);
    m_result_1 *= KFactor * p_fsrchannels->Weight();
    if (p_isrhandler->On()==3) m_result_2 = m_result_1;
    if (p_isrhandler->KMROn()==0) m_result_1 *= process->Differential(p_cms);
    else m_result_1 *= process->Differential(p_lab);
  }
  if (p_isrhandler->On()==3 && trigger==1) {
    Rotate(p_cms);
    p_isrhandler->CalculateWeight2(Q2);
    if (m_result_2 > 0.) m_result_2 *= process->Differential2();
    else m_result_2 = 0.;
  }
  if (m_nin>1 && (p_isrhandler->On()>0 || p_beamhandler->On()>0)) {
    m_flux=p_isrhandler->Flux();
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
  if (!(ATOOLS::IsZero((sin-sout)/(sin+sout)))) {
    ATOOLS::msg.Error()<<"Phase_Space_Handler::Check4Momentum(..): "
		       <<"Difference: "<<pin-pout<<std::endl;
    return false;
  }
  return true;
}

bool Phase_Space_Handler::SameEvent() 
{
  return OneEvent(-1.,1);
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::SameWeightedEvent() 
{
  return WeightedEvent(1);
}

bool Phase_Space_Handler::OneEvent(double mass,int mode)
{
  PROFILE_HERE;
  const bool use_overflow=true;

  if ((mass<0) && (!m_initialized)) InitIncoming();
  if ((mass>0) && (m_nin==1)) InitIncoming(mass);
  m_weight=1.;
  double value;
  for (int i=1;i<m_maxtrials+1;i++) {
    if (mode==0) {
      p_process->DeSelect();
      p_process->SelectOne();
    }
    else {
      if (!(p_process->Selected())) {
	msg.Error()<<" ERROR: in Phase_Space_Handler::OneEvent() "<<endl;
	return false;
      }
    }
    double max = p_process->Selected()->Max();

    if (p_process->Selected()->Overflow()==0.) {
      p_process->Selected()->RestoreInOrder(); // use last order for overflow events
      // *SH*    if (p_isrchannels) p_isrchannels->SetRange(p_isrhandler->SprimeRange(),p_isrhandler->YRange());
      value = Differential(p_process->Selected());
    }
    else {
      value = p_process->Selected()->Overflow();
      msg.Debugging()<<om::green<<" old Overflow event: "<<value<<" ("<<max<<") "<<om::bold<<value/max<<om::reset<<std::endl;
      msg.Debugging()<<"     proc "<<p_process->Selected()->Name()<<std::endl;

      std::vector<double> & xinfo =p_process->Selected()->XInfo();
      m_beamspkey[3] = xinfo[0];
      m_beamykey[2]  = xinfo[1];
      m_isrspkey[3]  = xinfo[2];
      m_isrykey[2]   = xinfo[3];
      p_isrhandler->MakeISR(p_lab,m_nvec);
      for (int i=0;i<m_nvec;++i)  p_lab[i]=p_process->Selected()->Momenta()[i];

      if (value > max) {
	p_process->Selected()->SetOverflow(value-max);
	value=max;
      }
      else {
	p_process->Selected()->SetOverflow(0.);
	msg.Debugging()<<om::blue<<" Overflow event finished: "<<value<<" ("<<max<<")"<<om::reset<<std::endl;
	msg.Debugging()<<"     proc "<<p_process->Selected()->Name()<<std::endl;
      }
      m_result_1=value;
      m_result_2=0.;
    }
    if (value > 0.) {
      double disc = 0.;
      if (value > max) {
	if (use_overflow) {
	  // don't use overflow
	  msg.Events()<<"Shifted maximum in "<<p_process->Selected()->Name()<<" : "
		      <<p_process->Selected()->Max()<<" -> "<<value<<endl;
	  p_process->Selected()->SetMax(value*1.001);
	  p_process->SetMax(0.);
	}
	else {
	  // use overflow
	  msg.Debugging()<<om::red<<" Overflow event: "<<value<<" ("<<max<<") "<<om::bold<<value/max<<om::reset<<std::endl;
	  msg.Debugging()<<"     proc "<<p_process->Selected()->Name()<<std::endl;
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
	  p_process->Selected()->SetMomenta(p_lab);
	  p_process->Selected()->SwapInOrder();
	}
	else {
	  p_process->Selected()->SetMomenta(p_lab);
	}
	return true;
      }
      else {
	msg.Tracking()<<" rejected "<<std::endl;
      }
    }
  }
  m_sumtrials += m_maxtrials;
  msg.Error()<<"Phase_Space_Handler::OneEvent() : "
	     <<" too many trials for "<<p_process->Selected()->Name()<<endl
	     <<"   Efficiency = "<<double(m_events)/double(m_sumtrials)*100.<<" %."<<endl;
  return false;
}

ATOOLS::Blob_Data_Base *Phase_Space_Handler::WeightedEvent(int mode)
{
  if (!m_initialized) InitIncoming();
  
  double value;
  for (int i=1;i<m_maxtrials+1;i++) {
    if (mode==0) {
      p_process->DeSelect();
      p_process->SelectOne();
    }
    else {
      if (!(p_process->Selected())) {
	msg.Error()<<" ERROR: in Phase_Space_Handler::WeightedEvent() "<<endl;
	return 0;
      }
    }
    p_process->Selected()->RestoreInOrder();
    value = Differential(p_process->Selected());
    if (value > 0.) {
      m_sumtrials+=i;
      ++m_events;
      if (m_result_1 < (m_result_1+m_result_2)*ATOOLS::ran.Get()) {
	Rotate(p_lab);
	p_process->Selected()->SetMomenta(p_lab);
	p_process->Selected()->SwapInOrder();
      }
      else {
	p_process->Selected()->SetMomenta(p_lab);
      }
      m_weight=value;
      m_trials=i;
      return new Blob_Data<Weight_Info>(Weight_Info(m_weight,m_trials));
    }
  }
  msg.Error()<<"Phase_Space_Handler::WeightedEvent() : "
	     <<" too many trials for "<<p_process->Selected()->Name()<<endl;
  m_weight=0.;
  return 0;
} 

void Phase_Space_Handler::AddPoint(const double value) 
{ 
  p_process->AddPoint(value); 
}

void Phase_Space_Handler::TestPoint(ATOOLS::Vec4D *const p)
{
  m_isrspkey[3] = sqr(ATOOLS::rpa.gen.Ecms());
  Single_Channel * TestCh = new Rambo(m_nin,m_nout,p_flavours);
  MakeIncoming(p);
  TestCh->GeneratePoint(p,p_cuts);
  delete TestCh;
}

void Phase_Space_Handler::WriteOut(const std::string &pID) 
{
  msg.Debugging()<<"Write out channels into directory : "<<pID<<endl;
  int  mode_dir = 448;
  mkdir(pID.c_str(),mode_dir); 
  if (p_beamchannels != 0) p_beamchannels->WriteOut(pID+string("/MC_Beam"));
  if (p_isrchannels  != 0) p_isrchannels->WriteOut(pID+string("/MC_ISR"));
  if (p_zchannels != 0) p_zchannels->WriteOut(pID+string("/MC_KMR_Z"));
  if (p_kpchannels!= 0) p_kpchannels->WriteOut(pID+string("/MC_KMR_KP"));
  if (p_fsrchannels  != 0) p_fsrchannels->WriteOut(pID+string("/MC_FSR"));
  string help     = (pID+string("/Random")).c_str();
  ran.WriteOutStatus(help.c_str());
}

bool Phase_Space_Handler::ReadIn(const std::string &pID) 
{
  msg.Debugging()<<"Read in channels from directory : "<<pID<<endl;
  bool okay = 1;
  if (p_beamchannels != 0) okay = okay && p_beamchannels->ReadIn(pID+string("/MC_Beam"));
  if (p_isrchannels  != 0) okay = okay && p_isrchannels->ReadIn(pID+string("/MC_ISR"));
  if (p_zchannels != 0) okay = okay && p_zchannels->ReadIn(pID+string("/MC_KMR_Z"));
  if (p_kpchannels!= 0) okay = okay && p_kpchannels->ReadIn(pID+string("/MC_KMR_KP"));
  if (p_fsrchannels  != 0) okay = okay && p_fsrchannels->ReadIn(pID+string("/MC_FSR"));
  if (rpa.gen.RandomSeed()==1234) {
    string filename     = (pID+string("/Random")).c_str();
    ran.ReadInStatus(filename.c_str(),0);
  }
  return okay;
}

void Phase_Space_Handler::Rotate(Vec4D *const p) 
{
  for (int i=0;i<m_nin+m_nout;i++) p[i] = Vec4D(p[i][0],(-1.)*Vec3D(p[i]));
}

bool Phase_Space_Handler::CreateIntegrators()
{
  /*if (p_cuts == 0) {
    p_cuts = new Cut_Data();
    p_cuts->Init(nin+nout,psflavs);
    if (proc->Selector()) (proc->Selector())->BuildCuts(p_cuts);
    } */
  if (m_nin==1) m_inttype = 0;
  if (m_nin==2) {
    if (p_beamhandler && p_beamhandler->On()>0) {
      if (!(MakeBeamChannels())) {
	msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		   <<"   did not construct any isr channels !"<<endl;
      }
    }
    if (p_isrhandler) {
      if (p_isrhandler->On()>0) {
	if (!(MakeISRChannels())) {
	  msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		     <<"   did not construct any isr channels !"<<endl;
	}
      }
      if (p_isrhandler->KMROn()) {
	if (!(MakeKMRChannels())) {
	  msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		     <<"   did not construct any kmr channels !"<<endl;
	}
      }
    }
  }
  if (m_nin==2) { 
    if (m_inttype < 3 || m_inttype == 5 && (p_fsrchannels!=0)) p_fsrchannels->DropAllChannels();
  }
  switch (m_inttype) {
  case 0: 
    if (m_nin==1 && m_nout==2) p_fsrchannels->Add(new Decay2Channel(m_nin,m_nout,p_flavours));
    else p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
    break;
  case 1: 
    p_fsrchannels->Add(new Sarge(m_nin,m_nout));
    break;
  case 2: 
    p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
    p_fsrchannels->Add(new Sarge(m_nin,m_nout));
    DropRedundantChannels();
    break;
  case 3: 
    p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
    DropRedundantChannels();
    break;
  case 4: 
    DropRedundantChannels();
    break;
  case 5:
    p_fsrchannels->Add(new RamboKK(m_nin,m_nout,p_flavours));
    break;    
  default:
    msg.Error()<<"Wrong phasespace integration switch ! Using RAMBO as default."<<endl;
    p_fsrchannels->Add(new Rambo(m_nin,m_nout,p_flavours));
  }  
  msg.Debugging()<<"Initialized Phase_Space_Integrator (\n\t";
  if (p_beamchannels) msg.Debugging()<<p_beamchannels->Name()<<","<<p_beamchannels->Number()<<";\n\t";
  if (p_isrchannels) msg.Debugging()<<p_isrchannels->Name()<<","<<p_isrchannels->Number()<<";\n\t";
  if (p_zchannels) msg.Debugging()<<p_zchannels->Name()<<","<<p_zchannels->Number()<<";\n\t";
  if (p_kpchannels) msg.Debugging()<<p_kpchannels->Name()<<","<<p_kpchannels->Number()<<";\n\t";
  if (p_fsrchannels) msg.Debugging()<<p_fsrchannels->Name()<<","<<p_fsrchannels->Number()<<")"<<endl;
  return 1;
}

void Phase_Space_Handler::DropRedundantChannels()
{
  p_fsrchannels->Reset();
  int number = p_fsrchannels->Number();
  msg.Debugging()<<"In Phase_Space_Handler::DropRedundantChannels"
		 <<"("<<p_fsrchannels->Name()<<")."<<endl
		 <<"    Start with "<<number<<" added channel(s)."<<endl;
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
	    msg.Debugging()<<"  "<<(p_fsrchannels->Channel(i))->Name()
			   <<" and "<<(p_fsrchannels->Channel(j))->Name()
			   <<" are identical."<<endl;
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
  msg.Debugging()<<"    "<<count<<" channel(s) were deleted."<<endl;
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
  for (int i=0;i<p_fsrchannels->Number();i++) {
    type=0; 
    mass=width=0.;
    if (p_process) p_fsrchannels->ISRInfo(i,type,mass,width);
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    bool add=true;
    for (size_t j=0;j<m_beamparams.size();j++) if (m_beamparams[j]==ci) add=false;
    if (add) m_beamparams.push_back(ci);
    ci.parameters.clear();
  }
  return CreateBeamChannels();
}

bool Phase_Space_Handler::MakeISRChannels()
{
  if (m_isrparams.size()>0) return CreateISRChannels();
  Channel_Info ci;
  if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) {
    // leptons : 1/s'^2 and 1/(s-s')^beta, sharp FW-BW peak
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(1.);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    ci.type = 0;
    (ci.parameters).push_back(2.);
    (ci.parameters).push_back(1.);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    ci.type = 3;
    (ci.parameters).push_back(p_isrhandler->Exponent(1));
    (ci.parameters).push_back(1.00000001);
    (ci.parameters).push_back(1.);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
  }
  else {
    // default : 1/s'
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(0.5);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    ci.type = 0;
    (ci.parameters).push_back(0.99);
    (ci.parameters).push_back(0.5);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
    ci.type = 0;
    (ci.parameters).push_back(2.);
    (ci.parameters).push_back(0.5);
    m_isrparams.push_back(ci);
    ci.parameters.clear();
  }
  int    type;
  double mass,width;
  for (int i=0;i<p_fsrchannels->Number();i++) {
    type=0; 
    mass=width=0.;
    p_fsrchannels->ISRInfo(i,type,mass,width);
    if (type==0 || type==3 ||
	(type==1 && (ATOOLS::IsZero(mass) || ATOOLS::IsZero(width))) ||
	(type==2 && ATOOLS::IsZero(mass))) continue;
    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((p_flavours[0].IsLepton()) || (p_flavours[1].IsLepton())) (ci.parameters).push_back(1.);
    else (ci.parameters).push_back(.5);
    bool add=true;
    for (size_t j=0;j<m_isrparams.size();j++) if (m_isrparams[j]==ci) add=false; 
    if (add) m_isrparams.push_back(ci);
    ci.parameters.clear();
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
      //MakeZChannels(1);
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
  Single_Channel * channel;   
  for (size_t i=0;i<m_beamparams.size();i++) {
    switch (m_beamparams[i].type) {
    case 0:
      // channel = new Simple_Pole_Uniform(m_beamparams[i].parameters[0],"",p_info);
      // p_beamchannels->Add(channel);
      channel = new Simple_Pole_Central(m_beamparams[i].parameters[0]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Simple_Pole_Forward(m_beamparams[i].parameters[0],
					m_beamparams[i].parameters[1]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Simple_Pole_Backward(m_beamparams[i].parameters[0],
					 m_beamparams[i].parameters[1]," beam",p_info);
      p_beamchannels->Add(channel);
      break;
    case 1:
      channel = new Resonance_Uniform(m_beamparams[i].parameters[0],
				      m_beamparams[i].parameters[1]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Resonance_Central(m_beamparams[i].parameters[0],
				      m_beamparams[i].parameters[1]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Resonance_Forward(m_beamparams[i].parameters[0],
				      m_beamparams[i].parameters[1],
				      m_beamparams[i].parameters[2]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Resonance_Backward(m_beamparams[i].parameters[0],
				       m_beamparams[i].parameters[1],
				       m_beamparams[i].parameters[2]," beam",p_info);
      p_beamchannels->Add(channel);
      break;
    case 2:
      channel = new Threshold_Central(m_beamparams[i].parameters[0]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Threshold_Forward(m_beamparams[i].parameters[0],
				      m_beamparams[i].parameters[2]," beam",p_info);
      p_beamchannels->Add(channel);
      channel = new Threshold_Backward(m_beamparams[i].parameters[0],
				       m_beamparams[i].parameters[2]," beam",p_info);
      p_beamchannels->Add(channel);
      break;
    case 3:
      if ((p_flavours[0].IsPhoton()) || (p_flavours[1].IsPhoton())) {
	channel = new LBS_Compton_Peak_Central(m_beamparams[i].parameters[1],
					       m_beamparams[i].parameters[0]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new LBS_Compton_Peak_Forward(m_beamparams[i].parameters[1],
					       m_beamparams[i].parameters[0],
					       m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
	channel = new LBS_Compton_Peak_Backward(m_beamparams[i].parameters[1],
						m_beamparams[i].parameters[0],
						m_beamparams[i].parameters[2]," beam",p_info);
	p_beamchannels->Add(channel);
      }
      break;
    }
  }
  return 1;
}

bool Phase_Space_Handler::CreateISRChannels()
{
  if (m_isrparams.size() < 1) return 0;
  Single_Channel * channel;   
  int length = m_isrparams.size();
  for (int i=0;i<length;i++) {
    switch (m_isrparams[i].type) {
    case 0:
      channel = new Simple_Pole_Uniform(m_isrparams[i].parameters[0]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Simple_Pole_Forward(m_isrparams[i].parameters[0],
					m_isrparams[i].parameters[1]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Simple_Pole_Backward(m_isrparams[i].parameters[0],
					 m_isrparams[i].parameters[1]," isr",p_info);
      p_isrchannels->Add(channel);
      break;
    case 1:
      channel = new Resonance_Uniform(m_isrparams[i].parameters[0],
				      m_isrparams[i].parameters[1]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Resonance_Forward(m_isrparams[i].parameters[0],
				      m_isrparams[i].parameters[1],
				      m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Resonance_Backward(m_isrparams[i].parameters[0],
				       m_isrparams[i].parameters[1],
				       m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
      break;
    case 2:
      channel = new Threshold_Central(m_isrparams[i].parameters[0]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Threshold_Forward(m_isrparams[i].parameters[0],
				      m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Threshold_Backward(m_isrparams[i].parameters[0],
				       m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
      break;
    case 3:
      channel = new Leading_Log_Central(m_isrparams[i].parameters[0],
					m_isrparams[i].parameters[1]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Leading_Log_Forward(m_isrparams[i].parameters[0],
					m_isrparams[i].parameters[1],
					m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
      channel = new Leading_Log_Backward(m_isrparams[i].parameters[0],
					 m_isrparams[i].parameters[1],
					 m_isrparams[i].parameters[2]," isr",p_info);
      p_isrchannels->Add(channel);
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
  double lambda2=ATOOLS::sqr(ATOOLS::Flavour(ATOOLS::kf::Z).Mass());
  lambda2*=exp(-M_PI/(MODEL::as->Beta0(lambda2)*(*MODEL::as)(lambda2)));
  channel = new LDL_KPerp(lambda2,"",p_info);
  p_kpchannels->Add(channel);
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
    ATOOLS::msg.Error()<<"Error in Phase_Space_Handler::Isrchannels("<<i<<")"<<std::endl
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
    ATOOLS::msg.Error()<<"Error in Phase_Space_Handler::Beamchannels("<<i<<")"<<std::endl
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

template Weight_Info ATOOLS::Blob_Data_Base::Get<Weight_Info>();

namespace ATOOLS {
  std::ostream & operator<<(std::ostream & s, const PHASIC::Weight_Info & wi)
  {
    s<<" weight="<<wi.weight<<"   ntrial="<<wi.ntrial<<std::endl;
    return s;
  }
}

template class Blob_Data<Weight_Info>;

