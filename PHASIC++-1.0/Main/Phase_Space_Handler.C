#include "Phase_Space_Handler.H"

#include "Phase_Space_Generator.H"
#include "Phase_Space_Integrator.H"
#include "Beam_Handler.H"
#include "ISR_Handler.H"
#include "Process_Base.H"
#include "XS_Base.H"


#include "Run_Parameter.H"
#include "Message.H"  
#include "Random.H"
#include "Rambo.H"
#include "Sarge.H"

using namespace PHASIC;
using namespace AMEGIC;
using namespace EXTRAXS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;

Phase_Space_Handler::Phase_Space_Handler(Process_Base * _proc,
					 ISR_Handler * _ih,Beam_Handler * _bh) 
  : proc(_proc), ih(_ih), bh(_bh)
{
  xs   = 0;
  nin  = proc->Nin();
  nout = proc->Nout();
  nvec = proc->Nvec()+1;
  name = proc->Name();
  Init(proc->Flavs());
}

Phase_Space_Handler::Phase_Space_Handler(XS_Base * _xs,
					 ISR_Handler * _ih,Beam_Handler * _bh) 
  : xs(_xs), ih(_ih), bh(_bh)
{
  proc = 0;
  nin  = xs->Nin();
  nout = xs->Nout();
  nvec = nin+nout;
  name = xs->Name();
  Init(xs->Flavs());
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p)            { delete [] p;         p            = 0; }
  if (psi)          { delete psi;          psi          = 0; }
  if (psgen)        { delete psgen;        psgen        = 0; }
  if (psflavs)      { delete [] psflavs;   psflavs      = 0; }
  if (isrchannels)  { delete isrchannels;  isrchannels  = 0; }
  if (fsrchannels)  { delete fsrchannels;  fsrchannels  = 0; }
  if (beamchannels) { delete beamchannels; beamchannels = 0; }
  // if (kin)          { delete kin;          kin          = 0; }
  if (proc) msg.Debugging()<<"Deleted Phase_Space_Handler for "<<proc->Name()<<endl;
  if (xs)   msg.Debugging()<<"Deleted Phase_Space_Handler for "<<xs->Name()<<endl;
}

void Phase_Space_Handler::Init(Flavour * _fl) {
  Data_Read dr(rpa.GetPath()+string("/Integration.dat"));

  error      = dr.GetValue<double>("ERROR");
  int_type   = dr.GetValue<int>("INTEGRATOR");
  
  psflavs    = new Flavour[nin+nout];
  for (int i=0;i<nin+nout;i++) psflavs[i] = _fl[i];
  p          = new vec4d[nvec];  
  msg.Debugging()<<"Initialize new vectors : "<<nvec<<endl;

  m1 = _fl[0].mass(); m12 = m1*m1;
  if (nin==2) {
    m2   = _fl[1].mass(); m22 = m2*m2;
  }

  E          = AORGTOOLS::rpa.gen.Ecms();
  s = sprime = E*E;

  maxtrials  = 100000;
  sumtrials  = 0;
  events     = 0;
  
  psi = 0; psgen = 0; beamchannels = 0; isrchannels = 0; fsrchannels = 0; 
  
  
  //  msg.Tracking()<<"Initialized a new Phase_Space_Handler(proc) for "<<proc->Name()<<endl;
  msg.Tracking()<<"   ("<<ih->Type()<<", "<<nin<<"  ->  "<<nout<<" process)"<<endl;
}

/* ----------------------------------------------------------------------

   Channel creation

   ---------------------------------------------------------------------- */

bool Phase_Space_Handler::CreateChannelLibrary(string ptype,string pID)
{
  if (nin ==1 ) return 1;  // Take RAMBO as default : no problems 

  msg.Tracking()<<"Creating Multichannel for phasespace integration for "
		<<"  "<<ptype<<"/"<<pID<<endl;   

  int ngraph = proc->NumberOfDiagrams();
    
  psgen       = new Phase_Space_Generator(nin,nout);
  fsrchannels = new Multi_Channel(string("fsr_")+proc->Name());

  bool newch  = psgen->Construct(fsrchannels,ptype,pID,psflavs,proc); 

  if (newch) {
    msg.Error()<<fsrchannels->Number()<<" new Channels produced for "<<pID<<" ! "<<endl
	       <<"After program termination please enter \"make install\" and rerun !"<<endl;
    return 0;
  }
  else {
    msg.Tracking()<<"No new Channels produced for "<<pID<<" ! "<<endl
		  <<"  added the following channels to the fs multi-channel : "<<endl;
    for (short int i=0;i<fsrchannels->Number();i++)
      msg.Tracking()<<"     "<<(fsrchannels->Channel(i))->Name()<<endl;
    msg.Debugging()<<"Program continues."<<endl;
    return 1;
  }
}

/* ----------------------------------------------------------------------

   Setting up the integrator

   ---------------------------------------------------------------------- */


bool Phase_Space_Handler::CreateIntegrators()
{
  msg.Debugging()<<"In Phase_Space_Handler::CreateIntegrators"<<endl;

  if (xs) psgen = new Phase_Space_Generator(xs->Nin(),xs->Nout());

  if (nin==1) int_type = 0;
  if (bh) {
    if ((nin==2) && bh && (bh->On()>0) ) {
      if (xs)   beamchannels = new Multi_Channel(string("beam_")+xs->Name());
      if (proc) beamchannels = new Multi_Channel(string("beam_")+proc->Name());

      if (!(MakeBeamChannels())) {
	msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		 <<"   did not construct any isr channels !"<<endl;
      }
      if (beamchannels) 
	msg.Debugging()<<"  ("<<beamchannels->Name()<<","<<beamchannels->Number()<<";";
    }
    else {
      msg.Debugging()<<" no Beam-Handling needed : "<<bh->Name()
		     <<" for "<<nin<<" incoming particles."<<endl;
    }
  }

  if ((nin==2) && (ih->On()>0)) {
    if (xs)   isrchannels = new Multi_Channel(string("isr_")+xs->Name());
    if (proc) isrchannels = new Multi_Channel(string("isr_")+proc->Name());

    if (!(MakeISRChannels())) {
      msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		 <<"   did not construct any isr channels !"<<endl;
    }
    if (isrchannels) 
      msg.Debugging()<<"  ("<<isrchannels->Name()<<","<<isrchannels->Number()<<";";
  }
  else {
    msg.Debugging()<<" no ISR needed : "<<ih->Name()
		   <<" for "<<nin<<" incoming particles."<<endl;
  }
  
  if (xs) {
    fsrchannels = new Multi_Channel(string("fsr_")+xs->Name());
    MakeFSRChannels();

    msg.Debugging()<<"Initialized Phase_Space_Integrator "<<endl<<"   (";
    if (isrchannels) msg.Debugging()<<" "<<isrchannels->Name()<<","<<isrchannels->Number()<<";";
    msg.Debugging()<<" "<<fsrchannels->Name()<<","<<fsrchannels->Number()<<")"<<endl;

    return 1;
  }


  if (proc) {
    msg.Debugging()<<" "<<fsrchannels->Name()<<","<<fsrchannels->Number()<<")"<<endl
		   <<"    integration mode = "<<int_type<<endl;
    
    if (int_type < 3) fsrchannels->DropAllChannels();
    
    switch (int_type) {
    case 0: 
      fsrchannels->Add(new Rambo(nin,nout,psflavs));
      break;
    case 1: 
      fsrchannels->Add(new Sarge(nin,nout));
      break;
    case 2: 
      fsrchannels->Add(new Rambo(nin,nout,psflavs));
      fsrchannels->Add(new Sarge(nin,nout));
      DropRedundantChannels();
      break;
    case 3: 
      fsrchannels->Add(new Rambo(nin,nout,psflavs));
      DropRedundantChannels();
      break;
    case 4: 
      DropRedundantChannels();
      break;
    default:
      msg.Error()<<"Wrong phasespace integration switch ! ";
      msg.Error()<<"Using as default : RAMBO ."<<endl;
      fsrchannels->Add(new Rambo(nin,nout,psflavs));
    }  
    
    msg.Debugging()<<"Initialized Phase_Space_Integrator (";
    if (beamchannels) msg.Debugging()<<beamchannels->Name()<<","<<beamchannels->Number()<<";";
    if (isrchannels)  msg.Debugging()<<isrchannels->Name()<<","<<isrchannels->Number()<<";";
    msg.Debugging()<<" "<<fsrchannels->Name()<<","<<fsrchannels->Number()<<")"<<endl;

    return 1;
  }
  return 0;
}


void Phase_Space_Handler::CollectChannels() {
  msg.Debugging()<<"Phase_Space_Handler::CollectChannels("<<proc->Name()<<") : "<<endl;
  fsrchannels = new Multi_Channel(string("fsr_")+proc->Name());
  msg.Debugging()<<"New multichannel init : "<<proc->Name()<<"  "<<fsrchannels<<endl;
  AddChannels(proc,fsrchannels,beam_params,isr_params);

  msg.Debugging()<<"    in total : "<<fsrchannels->Number()<<" FSR, ";
  msg.Debugging()<<3*isr_params.size()<<" ISR channels."<<endl;
}

void Phase_Space_Handler::AddChannels(Process_Base * _proc,Multi_Channel * _fsr,
				      vector<Channel_Info> & _beamparams,
				      vector<Channel_Info> & _isrparams) {
  bool         addit;
  Channel_Info ci;

  msg.Debugging()<<"In AddChannels("<<_proc->Name()<<", "<<_proc->Size()<<")"<<endl;

  for (int i=0;i<_proc->Size();i++) {
    if ((*_proc)[i]->Partner() == NULL) AddChannels((*_proc)[i],_fsr,_beamparams,_isrparams);
    else {
      msg.Debugging()<<"Test : "<<(*_proc)[i]->Name()<<" "<<(*_proc)[i]->Partner()->Name()<<endl;
      if ((*_proc)[i]->Partner() == (*_proc)[i]) {
	Single_Channel * sc;
	int next; string chname;
	msg.Debugging()<<"Phase_Space_Handler::Add "<<(*_proc)[i]->NumberOfFSRIntegrators()
		       <<" Channels of "<<(*_proc)[i]->Name()<<endl;
	for (int j=0;j<(*_proc)[i]->NumberOfFSRIntegrators();j++) { 
	  chname = ((*_proc)[i]->FSRIntegrator(j))->Name();
	  if ( (chname!=string("Rambo")) && (chname!=string("Sarge")) ) { 
	    next   = chname.find(string("--"));
	    chname = chname.substr(0,next);
	    sc   = psgen->SetChannel(nin,nout,psflavs,
				     ((*_proc)[i]->FSRIntegrator(j))->ChNumber(),chname);
	    sc->SetName(((*_proc)[i]->FSRIntegrator(j))->Name());
	    // sc = new Single_Channel((*_proc)[i]->FSRIntegrator(j));
	    _fsr->Add( sc );
	  }
	}

	if (bh->On()>0) {
	  msg.Debugging()<<"Phase_Space_Handler::Add "
			 <<(*_proc)[i]->NumberOfBeamIntegrators()<<" "
			 <<"Beam-Channels of "<<(*_proc)[i]->Name()<<endl;
	  for (int j=0;j<(*_proc)[i]->NumberOfBeamIntegrators()/3;j++) {
	    (*_proc)[i]->BeamChannels(j,ci);
	    addit = 1;
	    for (int k=0;k<beam_params.size();k++) {
	      if (beam_params[k]==ci) { addit = 0; break; }
	    }
	    if (addit) beam_params.push_back(ci);
	  }
	}

	if (ih->On()>0) {
	  msg.Debugging()<<"Phase_Space_Handler::Add "
			 <<(*_proc)[i]->NumberOfISRIntegrators()<<" "
			 <<"ISRChannels of "<<(*_proc)[i]->Name()<<endl;
	  for (int j=0;j<(*_proc)[i]->NumberOfISRIntegrators()/3;j++) {
	    (*_proc)[i]->ISRChannels(j,ci);
	    addit = 1;
	    msg.Debugging()<<j<<" th channel for "<<(*_proc)[i]->Name()<<" : "
			   <<" check for "<<ci.type<<" in "<<isr_params.size()<<endl; 
	    for (int k=0;k<isr_params.size();k++) {
	      if (isr_params[k]==ci) { addit = 0; break; }
	    }
	    if (addit) isr_params.push_back(ci);
	  }
	}
      }
      else {
	msg.Debugging()<<"Omit channels of "<<(*_proc)[i]->Name()<<", "
		       <<" partner = "<<(*_proc)[i]->Partner()->Name()<<endl;
      }
    }
  }
  msg.Debugging()<<"Leave AddChannels for "<<_proc->Name()<<endl;
};


bool Phase_Space_Handler::MakeFSRChannels()
{
  if (!xs) return 0;
  return psgen->CreateFSRChannels(fsrchannels,psflavs);
}


bool Phase_Space_Handler::MakeISRChannels()
{
  // For process groups : Harvest the members.
  if (!xs) {
    if ((proc) && (isr_params.size() > 0)) {
      return psgen->CreateISRChannels(isrchannels,psflavs,isr_params);
    }
  }

  Channel_Info ci;

  double deltay[2];
  
  deltay[0] = log(ih->Upper1());
  deltay[1] = log(ih->Upper2());

  msg.Out()<<"*** DeltaY1 / 2 = "<<deltay[0]<<" / "<<deltay[1]<<endl;


  if ((psflavs[0].islepton()) || (psflavs[1].islepton())) {
    // leptons : 1/s'^2 and 1/(s-s')^beta, sharp FW-BW peak
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(1.);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();

    ci.type = 0;
    (ci.parameters).push_back(2.);
    (ci.parameters).push_back(1.);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();

    ci.type = 3;
    (ci.parameters).push_back(ih->Exponent(1));
    (ci.parameters).push_back(1.);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();
  }
  else {
    // default : 1/s'
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();

    ci.type = 0;
    (ci.parameters).push_back(0.99);
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();
  }

  bool   addit;
  int    type,maxnumber;
  double mass,width;
    
  if (proc) { maxnumber = fsrchannels->Number(); }
  if (xs)   { maxnumber = xs->ISRNumber(); }
    
  for (int i=0;i<maxnumber;i++) {
    type = 0; mass = width = 0.;
    if (proc) fsrchannels->ISRInfo(i,type,mass,width);
    if (xs)   xs->ISRInfo(i,type,mass,width);
    
    msg.Debugging()<<i<<" : "<<type<<"/"<<mass<<"/"<<width<<endl;
    if (AMATOOLS::IsZero(mass) || AMATOOLS::IsZero(width)) continue;
    if ((type == 0) || (type == 3))                        continue;

    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((psflavs[0].islepton()) || (psflavs[1].islepton())) (ci.parameters).push_back(1.);
                                                       else (ci.parameters).push_back(.5);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);

    addit = 1;
    for (int j=0;j<isr_params.size();j++) {
      if (isr_params[j]==ci) { addit = 0; break; }
    }
    if (addit) isr_params.push_back(ci);
    ci.parameters.clear();
  }

  msg.Debugging()<<" create "<<3*isr_params.size()<<" ISR channels."<<endl;
  return psgen->CreateISRChannels(isrchannels,psflavs,isr_params);
}




bool Phase_Space_Handler::MakeBeamChannels()
{
  if (xs)  { msg.Error()<<"Beam Handling not implemented for EXTRA_XS !"<<endl; abort(); }

  // For process groups : Harvest the members.
  if (!xs) {
    if ((proc) && (beam_params.size() > 0)) {
      return psgen->CreateBeamChannels(beamchannels,psflavs,beam_params);
    }
  }

  double deltay[2];
  deltay[0] = log(bh->Upper1());
  deltay[1] = log(bh->Upper2());

  Channel_Info ci;

  // default : Beamstrahlung
  if ((psflavs[0].islepton()) && (psflavs[1].islepton())) {
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(bh->Exponent(1));
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    beam_params.push_back(ci);
    ci.parameters.clear();
  }

  // Laser Backscattering spectrum
  if ((psflavs[0].isphoton()) || (psflavs[1].isphoton())) {
    ci.type = 0;
    (ci.parameters).push_back(.5);
    (ci.parameters).push_back(1.);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    beam_params.push_back(ci);
    ci.parameters.clear();

    ci.type = 3;
    (ci.parameters).push_back(bh->Peak());
    (ci.parameters).push_back(bh->Exponent(1));
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    (ci.parameters).push_back(0.7);
    beam_params.push_back(ci);
    ci.parameters.clear();
  }

  bool   addit;
  int    type,maxnumber;
  double mass,width;
    
  if (proc) { maxnumber = fsrchannels->Number(); }
  //if (xs)   { maxnumber = xs->BeamNumber(); }
    
  for (int i=0;i<maxnumber;i++) {
    type = 0; mass = width = 0.;
    if (proc) fsrchannels->ISRInfo(i,type,mass,width);
    //if (xs)   xs->ISRInfo(i,type,mass,width);
    
    msg.Debugging()<<i<<" : "<<type<<"/"<<mass<<"/"<<width<<endl;
    if (AMATOOLS::IsZero(mass) || AMATOOLS::IsZero(width)) continue;
    if ((type == 0) || (type == 3))                        continue;

    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((psflavs[0].islepton()) || (psflavs[1].islepton())) (ci.parameters).push_back(1.);
                                                       else (ci.parameters).push_back(.5);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);

    addit = 1;
    for (int j=0;j<beam_params.size();j++) {
      if (beam_params[j]==ci) { addit = 0; break; }
    }
    if (addit) beam_params.push_back(ci);
    ci.parameters.clear();
  }

  msg.Debugging()<<" create "<<3*beam_params.size()<<" Beam channels."<<endl;
  return psgen->CreateBeamChannels(beamchannels,psflavs,beam_params);
}



/* ----------------------------------------------------------------------

   Integration

   ---------------------------------------------------------------------- */

double Phase_Space_Handler::Integrate() 
{
  psi        = new Phase_Space_Integrator();
  
  msg.Debugging()<<"Phase_Space_Handler::Integrate with : "<<endl;
  if (beamchannels) 
    msg.Debugging()<<"  Beam : "<<beamchannels->Name()<<" ("<<beamchannels<<") "
		   <<"  ("<<beamchannels->Number()<<","<<beamchannels->N()<<")"<<endl;
  if (isrchannels) 
    msg.Debugging()<<"  ISR  : "<<isrchannels->Name()<<" ("<<isrchannels<<") "
		   <<"  ("<<isrchannels->Number()<<","<<isrchannels->N()<<")"<<endl;

  msg.Debugging()<<"  FSR  : "<<fsrchannels->Name()<<" ("<<fsrchannels<<") "
		 <<"  ("<<fsrchannels->Number()<<","<<fsrchannels->N()<<")"<<endl;
  

  sprime     = sqr(AORGTOOLS::rpa.gen.Ecms());
  if (!(MakeIncoming(p)) ) {
    msg.Error()<<"Phase_Space_Handler::Integrate : Error !"<<endl
	       <<"  Either too little energy for initial state"
	       <<"  ("<<E<<" vs "<<m1+m2<<") or "<<endl
	       <<"  bad number of incoming particles ("<<nin<<")."<<endl;
    return 0;
  } 

  msg.SetPrecision(12);
  if (bh) 
  if (bh->On()>0) {
    beamchannels->GetRange();
    beamchannels->SetRange(bh->SprimeRange(),bh->YRange());
    beamchannels->GetRange();
  }
  if (ih->On()>0) {
    isrchannels->GetRange();
    isrchannels->SetRange(ih->SprimeRange(),ih->YRange());
    isrchannels->GetRange();
  }
  msg.SetPrecision(6);

  return psi->Calculate(this,error);
}


bool Phase_Space_Handler::MakeIncoming(vec4d * _p) {
  if (nin == 1) {
    E     = m1;
    s     = E*E;
    _p[0] = vec4d(E,0.,0.,0.);
    return 1;

    flux = 1./(2.*m1);
  }
  if (nin == 2) {
    double Eprime = sqrt(sprime);
    if ((E<m1+m2)) return 0;
    double x      = 1./2.+(m12-m22)/(2.*sprime);
    double E1     = x*Eprime;
    double E2     = (1.-x)*Eprime;
    _p[0]         = vec4d(E1,0.,0.,sqrt(sqr(E1)-sqr(m1)));
    _p[1]         = vec4d(E2,(-1.)*vec3d(_p[0]));
    
    flux          = 1./(2.*sqrt(sqr(sprime-m12-m22)-4.*m12*m22));

    return 1;
  }
  return 0;
} 

double Phase_Space_Handler::Differential() { 
  if (proc) return Differential(proc);
  if (xs)   return Differential(xs);
  return 0.;
}




double Phase_Space_Handler::Differential(Process_Base * process)
{
  y = 0;
  if (bh->On()>0) { 
    beamchannels->GeneratePoint(sprimeB,yB,bh->On()); 
    if (!(bh->MakeBeams(p,sprimeB,yB))) return 0.;
    if (ih->On()>0) {
      ih->SetSprimeMax(sprimeB*sqrt(ih->Upper1()*ih->Upper2()));
      ih->SetPole(sprimeB);
      isrchannels->SetRange(ih->SprimeRange(),ih->YRange());
    }
    //msg.Out()<<"Beam : "<<sprimeB<<" / "<<yB<<endl;
    sprime = sprimeB; y += yB;
  }

  if (ih->On()>0) { 
    isrchannels->GeneratePoint(sprimeI,yI,ih->On());
    //msg.Out()<<"ISR : "<<sprimeI<<" / "<<yI<<endl;
    if (!(ih->MakeISR(p,sprimeI,yI))) return 0.;
    sprime = sprimeI; y += yI;
  }

  if ( (bh->On()>0) || (ih->On()>0) ) {
    proc->UpdateCuts(sprime,y);
  }
  
  //msg.Out()<<"Before FSR"<<endl;
  //for (int i=0;i<nin+nout;i++) msg.Out()<<" "<<i<<"th : "<<p[i]<<" "<<p[i].abs2()<<endl;

  fsrchannels->GeneratePoint(p,proc->Cuts());

  if (!Check4Momentum(p)) return 0.;

  double value = 0., KFactor = 0., Q2 = -1.;
  bool take = 1;

  result1 = result2 = 0.;

  if (bh->On()>0) bh->BoostInLab(p,nin+nout);
  if (ih->On()>0) ih->BoostInLab(p,nin+nout);

  // First part : flin[0] coming from Beam[0] and flin[1] coming from Beam[1]

  bool trigger = 0;
  if ( (proc->Selector())->Trigger(p)) {
    trigger = 1;
    result1 = 1.;
    Q2 = proc->Scale(p);
    //msg.Out()<<"Scale = "<<Q2<<endl;
    if (ih->On()>0) {
      ih->CalculateWeight(Q2);
      isrchannels->GenerateWeight(sprimeI,yI,ih->On());
      result1 *= isrchannels->Weight();
      ih->BoostInCMS(p,nin+nout);
    }
    if (bh->On()>0) {
      bh->CalculateWeight(Q2);
      beamchannels->GenerateWeight(sprimeB,yB,bh->On());
      result1 *= beamchannels->Weight() * bh->Weight();
      bh->BoostInCMS(p,nin+nout);
    }

    KFactor = proc->KFactor(Q2);
    fsrchannels->GenerateWeight(p,proc->Cuts());
    result1 *= KFactor = fsrchannels->Weight();
    
    if (ih->On()==3) result2 = result1;
    
    result1 *= process->Differential(p);
  }
  // Second part : flin[0] coming from Beam[1] and flin[1] coming from Beam[0]
  if (ih->On()==3) {
    Rotate(p);
    if ( (proc->Selector())->Trigger(p)) {
      ih->CalculateWeight2(Q2);
      if (result2 > 0.) result2 *= process->Differential2();
      else              result2  = 0.;
    }
    else                result2  = 0.;
    Rotate(p);
  }
  if ( (ih->On()>0) || (bh->On()>0) ) 
    flux = 1./(2.*sqrt(sqr(sprime-m12-m22)-4.*m12*m22));
  //  cout<<"weight"<<x1<<" Diff"<<x2<<endl;

  return flux*(result1+result2);
}


bool Phase_Space_Handler::Check4Momentum(vec4d * _p) {
  vec4d pin,pout;
  pin = pout = vec4d(0.,0.,0.,0.);
  for (int i=0;i<nin;i++)        pin  += _p[i];
  for (int i=nin;i<nin+nout;i++) pout += _p[i];
  double sin = pin.abs2(),sout = pout.abs2();
  if (!(AMATOOLS::IsZero((sin-sout)/(sin+sout)))) return 0;
  return 1;
}


double Phase_Space_Handler::Differential(XS_Base * xsec) {
  y = 0;
  if (ih->On()>0) { 
    isrchannels->GeneratePoint(sprimeI,yI,ih->On());
    //msg.Out()<<"ISR : "<<sprimeI<<" / "<<yI<<endl;
    if (!(ih->MakeISR(p,sprimeI,yI))) return 0.;
    sprime = sprimeI; y += yI;
  }

//   if ((ih->On()>0) ) {
//     xsec->UpdateCuts(sprime,y);
//   }
  
  //msg.Out()<<"Before FSR"<<endl;
  //for (int i=0;i<nin+nout;i++) msg.Out()<<" "<<i<<"th : "<<p[i]<<" "<<p[i].abs2()<<endl;

  //  fsrchannels->GeneratePoint(p,xsec->Cuts());
  fsrchannels->GeneratePoint(p);

  if (!Check4Momentum(p)) {
    msg.Out()<<" WARNING: Check4Momentum(p) failed "<<endl;
    return 0.;
  }

  double value = 0., KFactor = 0., Q2 = -1.;
  bool take = 1;

  result1 = result2 = 0.;

  if (ih->On()>0) ih->BoostInLab(p,nin+nout);

  // First part : flin[0] coming from Beam[0] and flin[1] coming from Beam[1]

  bool trigger = 0;
  if ( (xsec->Selector())->Trigger(p)) {
    trigger = 1;
    result1 = 1.;
    Q2 = xsec->Scale(p);
    //msg.Out()<<"Scale = "<<Q2<<endl;
    if (ih->On()>0) {
      ih->CalculateWeight(Q2);
      isrchannels->GenerateWeight(sprimeI,yI,ih->On());
      result1 *= isrchannels->Weight();
      ih->BoostInCMS(p,nin+nout);
    }

    KFactor = xsec->KFactor(Q2);
    fsrchannels->GenerateWeight(p);
    result1 *= KFactor = fsrchannels->Weight();
    
    if (ih->On()==3) result2 = result1;
    
    result1 *= xsec->Differential(p);
  }
  
  //  cout<<"weight"<<x1<<" Diff"<<x2<<endl;
  if ( (ih->On()>0) ) 
    flux = 1./(2.*sqrt(sqr(sprime-m12-m22)-4.*m12*m22));

  return flux*(result1+result2);
 }


bool Phase_Space_Handler::SameEvent() {
  return OneEvent(1);
}

bool Phase_Space_Handler::OneEvent(int mode)
{
  double value;
  for (int i=1;i<maxtrials+1;i++) {
    if (proc) {
      if (mode==0) {
	proc->DeSelect();
	proc->SelectOne();
      }
      else {
	if (!(proc->Selected())) {
	  msg.Error()<<" ERROR: in Phase_Space_Handler::OneEvent() "<<endl;
	  return 0;
	}
      }
      value = Differential(proc->Selected());
    }
    if (xs) {
      xs->DeSelect();
      xs->SelectOne();
      value = Differential(xs->Selected());
    }

    if (value > 0.) {
      double max;
      double disc = 0.;
      if (proc) max = proc->Selected()->Max();
      if (xs)   max = xs->Selected()->Max();
      if (value > max) {
	if (proc) {
	  msg.Events()<<"Shifted maximum in "<<proc->Selected()->Name()<<" : "
		      <<proc->Selected()->Max()<<" -> "<<value<<endl;
	  proc->Selected()->SetMax(value);
	}
	if (xs) {
	  msg.Events()<<"Shifted maximum in "<<xs->Selected()->Name()<<" : "
		      <<xs->Selected()->Max()<<" -> "<<value<<endl;
	  xs->Selected()->SetMax(value);
	}
      }
      else disc  = max*AMATOOLS::Ran.get();
      if (value >= disc) {
	sumtrials += i;events ++;
	msg.Debugging()<<"Phase_Space_Handler::OneEvent() : "<<i<<" trials for ";
	if (proc) msg.Debugging()<<proc->Selected()->Name()<<endl;
	if (xs)   msg.Debugging()<<xs->Selected()->Name()<<endl;
	msg.Debugging()<<"   Efficiency = "<<100./double(i)<<" %."
		       <<"   in total = "<<double(events)/double(sumtrials)*100.<<" %."<<endl;


	if (result1 < (result1+result2)*AMATOOLS::Ran.get()) Rotate(p);
	if (proc) proc->Selected()->SetMomenta(p);
	if (xs)   xs->Selected()->SetMomenta(p);
	for (int i=0;i<nin+nout;i++) {
	  if (proc) msg.Debugging()<<"  "<<proc->Selected()->Flavs()[i]<<" : "<<p[i]<<endl; 
	  if (xs)   msg.Debugging()<<"  "<<xs->Selected()->Flavs()[i]<<" : "<<p[i]<<endl; 
	} 
	return 1.;
      }
    }
  }
  sumtrials += maxtrials;


  msg.Debugging()<<"Phase_Space_Handler::OneEvent() : "
		 <<" too many trials for ";
  if (proc) msg.Debugging()<<proc->Selected()->Name()<<endl;
  if (xs)   msg.Debugging()<<xs->Selected()->Name()<<endl;
  msg.Debugging()<<"   Efficiency = "<<double(events)/double(sumtrials)*100.<<" %."<<endl;


  return 0;
}


double Phase_Space_Handler::WeightedEvent()
{
  return 0.;
} 





void Phase_Space_Handler::AddPoint(const double value) { 
  if (proc) proc->AddPoint(value); 
  if (xs)   xs->AddPoint(value); 
}












void Phase_Space_Handler::DropRedundantChannels()
{
  fsrchannels->Reset();
  int number = fsrchannels->Number();

  msg.Debugging()<<"In Phase_Space_Handler::DropRedundantChannels";
  msg.Debugging()<<"("<<fsrchannels->Name()<<")."<<endl;
  msg.Debugging()<<"    Start with "<<number<<" added channel(s)."<<endl;

  if (number<2) return;
  vec4d** perm_vec = new vec4d*[number]; 
  
  for (short int i=0;i<number;i++) perm_vec[i] = new vec4d[nin+nout+1];
  
  // Create Momenta
  int rannum   = 1 + 2 + 3*(nout-2);
  double * ran = new double[rannum];
  for (short int i=0;i<rannum;i++) ran[i] = Ran.get();  
  // Init markers for deletion and results to compare.
  int    * marker = new int[number];  
  double * res    = new double[number];
  for (short int i=0;i<number;i++) { marker[i] = 0;res[i] = 0.; }

  for (short int i=0;i<number;i++) {
    perm_vec[i][0] = vec4d(rpa.gen.Ecms()/2.,0.,0.,rpa.gen.Ecms()/2.);
    perm_vec[i][1] = vec4d(rpa.gen.Ecms()/2.,0.,0.,-rpa.gen.Ecms()/2.); 
    msg.Debugging()<<"==== "<<i<<" : ";
    msg.Debugging()<<(fsrchannels->Channel(i))->Name()<<"====================="<<endl;
    fsrchannels->GeneratePoint(i,perm_vec[i],proc->Cuts(),ran);
    // for (short int j=0;j<nin+nout;j++) msg.Debugging()<<j<<"th : "<<perm_vec[i][j]<<endl;
    fsrchannels->GenerateWeight(i,perm_vec[i],proc->Cuts());
    res[i] = fsrchannels->Weight();
    if (res[i]==0.) {
      msg.Debugging()<<"  "<<(fsrchannels->Channel(i))->Name()<<" produced a zero weight."<<endl;
      marker[i] = 1;
    }
  }
  delete[] ran;

  // kick identicals & permutations
  for (short int i=0;i<number;i++) {
    if (marker[i]==0) {
      for (short int j=i+1;j<number;j++) {
	if (marker[j]==0) {
	  if ( (Compare(perm_vec[i],perm_vec[j])) && 
	       (AMATOOLS::IsEqual(res[i],res[j])) ) {
	    msg.Debugging()<<"  "<<(fsrchannels->Channel(i))->Name()
			   <<" and "<<(fsrchannels->Channel(j))->Name()
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
    if (sqr(fl_res[j].mass())>ycut*sqr(rpa.gen.Ecms()) &&
    sqr(fl_res[j].mass())<sqr(rpa.gen.Ecms())) 
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
      fsrchannels->DropChannel(i-count);
      count++;
    }
  }
  msg.Debugging()<<"    "<<count<<" channel(s) were deleted."<<endl;


  delete [] res;
  delete [] marker; 
  for (short int i=0;i<number;i++) delete [] perm_vec[i];
  delete [] perm_vec; 
}



/* ----------------------------------------------------------------------

   One test point to check the amplitudes

   ---------------------------------------------------------------------- */

void Phase_Space_Handler::TestPoint(vec4d * _p)
{
  sprime                  = sqr(AORGTOOLS::rpa.gen.Ecms());
  Single_Channel * TestCh = new Rambo(nin,nout,psflavs);
  MakeIncoming(_p);
  TestCh->GeneratePoint(_p,proc->Cuts());
  delete TestCh;
  msg.Debugging()<<"------------------------------------------------"<<endl;
  for (int i=0;i<nin+nout;i++) msg.Debugging()<<i<<" th mom : "<<_p[i]<<endl;
  msg.Debugging()<<"------------------------------------------------"<<endl;
}

/* ----------------------------------------------------------------------

   Phase space I/O

   ---------------------------------------------------------------------- */

void Phase_Space_Handler::WriteOut(string pID) {
  msg.Debugging()<<"Write out channels into directory : "<<pID<<endl;
  int  mode_dir = 448;
  mkdir(pID.c_str(),mode_dir); 
  if (beamchannels != 0) beamchannels->WriteOut(pID+string("/MC_Beam"));
  if (isrchannels  != 0) isrchannels->WriteOut(pID+string("/MC_ISR"));
  if (fsrchannels  != 0) fsrchannels->WriteOut(pID+string("/MC_FSR"));

  char * filename = new char[100];
  string help     = (pID+string("/Random")).c_str();
  for (int pos=0;pos<help.length();pos++) filename[pos] = help[pos];
  filename[help.length()]=0;
  int nran = Ran.WriteOutStatus(filename);
}

bool Phase_Space_Handler::ReadIn(string pID) {
  msg.Debugging()<<"Read in channels from directory : "<<pID<<endl;
  bool okay = 1;

  if (beamchannels != 0) okay = okay && beamchannels->ReadIn(pID+string("/MC_Beam"));
  if (isrchannels  != 0) okay = okay && isrchannels->ReadIn(pID+string("/MC_ISR"));
  if (fsrchannels  != 0) okay = okay && fsrchannels->ReadIn(pID+string("/MC_FSR"));

  char * filename = new char[100];
  string help     = (pID+string("/Random")).c_str();
  for (int pos=0;pos<help.length();pos++) filename[pos] = help[pos];

  Ran.ReadInStatus(filename,0);

  return okay;
}

void Phase_Space_Handler::Rotate(vec4d * _p) 
{
  for (int i=0;i<nin+nout;i++) _p[i] = vec4d(_p[i][0],(-1.)*vec3d(_p[i]));
}


bool Phase_Space_Handler::Compare(vec4d* _p1,vec4d* _p2)
{
  if (nout==2) {
    for (short int i=0;i<nout;i++) { 
      if (_p1[nin+i] != _p2[nin+i]) return 0;
    }
    return 1;
  }
  else {
    //Identicals
    for (short int i=0;i<nout;i++) {
      if (_p1[i+nin] != _p2[i+nin]) return 0;
    }
    return 1;

    //Permutations - not reached right now.
    int * perm = new int[nout];
    for (short int i=0;i<nout;i++) perm[i] = 0; 
    
    int over = 0;
    int hit,sw1;
    for(;;) {
      sw1 = 1;
      for(short int i=0;i<nout;i++) {
	for (short int j=i+1;j<nout;j++) 
	  if (perm[i]==perm[j]) {sw1 = 0; break;}
      }    
      if (sw1) {
	hit = 1;
	for (short int i=0;i<nout;i++) {
	  if (_p1[i+nin] != _p2[perm[i]+nin]) {
	    hit = 0;
	    break;
	  }
	}
	if (hit) return 1;
      }
      for (short int j=nout-1;j>=0;j--) {
	if ((perm[j]+1)<nout) {
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
