#include "Phase_Space_Handler.H"

#include "Phase_Space_Integrator.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Rambo.H"
#include "RamboKK.H"
#include "Sarge.H"
#include "ISR_Channel.H"
#include "FSR_Channel.H"

#include "Run_Parameter.H"
#include "Message.H"  
#include "Random.H"

using namespace PHASIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;

Phase_Space_Handler::Phase_Space_Handler(Integrable_Base * _proc,
					 ISR_Handler * _ih,Beam_Spectra_Handler * _bh) 
  : proc(_proc), ih(_ih), bh(_bh),
    E(AORGTOOLS::rpa.gen.Ecms()), s(E*E), sprime(s),
    maxtrials(100000), sumtrials(0),
    events(0), psi(NULL), p(NULL),
    beamchannels(NULL), isrchannels(NULL), fsrchannels(NULL), psflavs(NULL),
    nin(proc->Nin()), nout(proc->Nout()), nvec(proc->Nvec()+1), name(proc->Name())
{
  Data_Read dr(rpa.GetPath()+string("/Integration.dat"));
  
  error      = dr.GetValue<double>("ERROR",0.01);
  int_type   = dr.GetValue<int>("INTEGRATOR",3);
  
  psflavs    = new Flavour[nin+nout];
  for (int i=0;i<nin+nout;i++) psflavs[i] = proc->Flavs()[i];
  p          = new Vec4D[nvec];  
  msg.Debugging()<<"Initialize new vectors : "<<nvec<<endl;
  
  m1 = psflavs[0].Mass(); m12 = m1*m1;
  if (nin==2) { m2   = psflavs[1].Mass(); m22 = m2*m2; }

  if (nin==2) {
    if (bh) {
      if (bh->On()>0) beamchannels = new Multi_Channel(string("beam_")+proc->Name());
    }
    if (ih) {
      if (ih->On()>0) isrchannels  = new Multi_Channel(string("isr_")+proc->Name());
    }
  }
  fsrchannels = new Multi_Channel(string("fsr_")+proc->Name());

  msg.Tracking()<<"Initialized new Phase_Space_Handler for "<<proc->Name()<<endl
		<<" ("<<ih->Type()<<", "<<nin<<"  ->  "<<nout<<" process)"<<endl;
}

Phase_Space_Handler::~Phase_Space_Handler()
{
  if (p)            { delete [] p;         p            = 0; }
  if (psi)          { delete psi;          psi          = 0; }
  if (psflavs)      { delete [] psflavs;   psflavs      = 0; }
  if (isrchannels)  { delete isrchannels;  isrchannels  = 0; }
  if (fsrchannels)  { delete fsrchannels;  fsrchannels  = 0; }
  if (beamchannels) { delete beamchannels; beamchannels = 0; }
  if (proc) msg.Debugging()<<"Deleted Phase_Space_Handler for "<<proc->Name()<<endl;
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
  

  if (!(MakeIncoming(p)) ) {
    msg.Error()<<"Phase_Space_Handler::Integrate : Error !"<<endl
	       <<"  Either too little energy for initial state"
	       <<"  ("<<E<<" vs "<<m1+m2<<") or "<<endl
	       <<"  bad number of incoming particles ("<<nin<<")."<<endl;
    return 0;
  } 

  msg.SetPrecision(12);
  if (bh) {
    if (bh->On()>0) {
      beamchannels->GetRange();
      beamchannels->SetRange(bh->SprimeRange(),bh->YRange());
      beamchannels->GetRange();
    }
  }
  if (ih) {
    if (ih->On()>0) {
      ih->SetSprimeMin(sqr(proc->ISRThreshold()));
      isrchannels->GetRange();
      msg.Debugging()<<"In Phase_Space_Handler::Integrate : "<<bh->On()<<":"<<ih->On()<<endl
		     <<"   "<<ih->SprimeMin()<<" ... "<<ih->SprimeMax()<<" ... "<<ih->Pole()<<endl;
      isrchannels->SetRange(ih->SprimeRange(),ih->YRange());
      isrchannels->GetRange();
    }
  }
  msg.SetPrecision(6);

  if (nin==2) return psi->Calculate(this,error);
  if (nin==1) return psi->CalculateDecay(this,sqrt(p[0].Abs2()),error);
  return 0.;
}


bool Phase_Space_Handler::MakeIncoming(Vec4D * _p) {
  if (nin == 1) {
    E     = m1;
    s     = E*E;
    _p[0] = Vec4D(E,0.,0.,0.);

    flux = 1./(2.*m1);
    return 1;
  }
  if (nin == 2) {
    double Eprime = sqrt(sprime);
    if ((E<m1+m2)) return 0;
    double x      = 1./2.+(m12-m22)/(2.*sprime);
    double E1     = x*Eprime;
    double E2     = (1.-x)*Eprime;
    _p[0]         = Vec4D(E1,0.,0.,sqrt(sqr(E1)-sqr(m1)));
    _p[1]         = Vec4D(E2,(-1.)*Vec3D(_p[0]));
    
    flux          = 1./(2.*sqrt(sqr(sprime-m12-m22)-4.*m12*m22));
    return 1;
  }
  return 0;
} 

double Phase_Space_Handler::Differential() { 
  return Differential(proc);
}

double Phase_Space_Handler::Differential(Integrable_Base * process) { 
  y = 0;
  if (bh && bh->On()>0) { 
    beamchannels->GeneratePoint(sprimeB,yB,bh->On()); 
    if (!(bh->MakeBeams(p,sprimeB,yB))) return 0.;
    if (ih && ih->On()>0) {
      ih->SetSprimeMax(sprimeB*sqrt(ih->Upper1()*ih->Upper2()));
      ih->SetPole(sprimeB);
      isrchannels->SetRange(ih->SprimeRange(),ih->YRange());
    }

    //msg.Out()<<"Beam : "<<sprimeB<<" / "<<yB<<endl;

    sprime = sprimeB; y += yB;
  }

  if (ih && ih->On()>0) { 
    isrchannels->GeneratePoint(sprimeI,yI,ih->On());

    //msg.Out()<<"ISR : "<<sprimeI<<" / "<<yI<<endl;

    if (!(ih->MakeISR(p,sprimeI,yI))) return 0.;
    sprime = sprimeI; y += yI;
  }

  if ( (bh && bh->On()>0) || (ih && ih->On()>0) ) {
    proc->UpdateCuts(sprime,y);
  }
  fsrchannels->GeneratePoint(p,proc->Cuts());

  if (!Check4Momentum(p)) {
    msg.Out()<<"Phase_Space_Handler Check4Momentum(p) failed"<<endl;
    return 0.;
  }

  double value = 0., KFactor = 0., Q2 = -1.;
  result1      = result2     = 0.;

  if (bh->On()>0) bh->BoostInLab(p,nin+nout);
  if (ih->On()>0) ih->BoostInLab(p,nin+nout);

  // First part : flin[0] coming from Beam[0] and flin[1] coming from Beam[1]

  bool trigger = 0;
 
  if ((proc->Selector())->Trigger(p)) {
    trigger = 1;
    result1 = 1.;
    Q2 = proc->Scale(p);
    if (ih && ih->On()>0) {
      ih->CalculateWeight(Q2);
      isrchannels->GenerateWeight(sprimeI,yI,ih->On());
      result1 *= isrchannels->Weight();
      ih->BoostInCMS(p,nin+nout);
    }
    if (bh && bh->On()>0) {
      bh->CalculateWeight(Q2);
      beamchannels->GenerateWeight(sprimeB,yB,bh->On());
      result1 *= beamchannels->Weight() * bh->Weight();
      bh->BoostInCMS(p,nin+nout);
    }

    KFactor = proc->KFactor(Q2);
    fsrchannels->GenerateWeight(p,proc->Cuts());
    result1 *= KFactor * fsrchannels->Weight();

    if (ih && ih->On()==3) result2 = result1;

    result1 *= process->Differential(p);
  }

  // Second part : flin[0] coming from Beam[1] and flin[1] coming from Beam[0]
  if (ih && ih->On()==3 && trigger==1) {
    Rotate(p);
    ih->CalculateWeight2(Q2);
    if (result2 > 0.) result2 *= process->Differential2();
                 else result2  = 0.;
  }

  if ( (ih && ih->On()>0) || (bh && bh->On()>0) ) 
    flux = 1./(2.*sqrt(sqr(sprime-m12-m22)-4.*m12*m22));

  return flux*(result1+result2);
}

bool Phase_Space_Handler::Check4Momentum(Vec4D * _p) {
  Vec4D pin,pout;
  pin = pout = Vec4D(0.,0.,0.,0.);
  for (int i=0;i<nin;i++)        pin  += _p[i];
  for (int i=nin;i<nin+nout;i++) pout += _p[i];
  double sin = pin.Abs2(),sout = pout.Abs2();
  if (!(AMATOOLS::IsZero((sin-sout)/(sin+sout)))) return 0;
  return 1;
}

bool Phase_Space_Handler::SameEvent() {
  return OneEvent(1);
}

bool Phase_Space_Handler::OneEvent(int mode)
{
  double value;
  for (int i=1;i<maxtrials+1;i++) {
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
    ih->SetSprimeMin(sqr(proc->ISRThreshold()));
    
    if (isrchannels) isrchannels->SetRange(ih->SprimeRange(),ih->YRange());
    value = Differential(proc->Selected());

    if (value > 0.) {
      double max;
      double disc = 0.;
      max = proc->Selected()->Max();
      if (value > max) {
	  msg.Events()<<"Shifted maximum in "<<proc->Selected()->Name()<<" : "
		      <<proc->Selected()->Max()<<" -> "<<value<<endl;
	  proc->Selected()->SetMax(value);
      }
      else disc  = max*AMATOOLS::ran.Get();
      if (value >= disc) {
	sumtrials += i;events ++;
	msg.Debugging()<<"Phase_Space_Handler::OneEvent() : "<<i<<" trials for "
		       <<proc->Selected()->Name()<<endl
		       <<"   Efficiency = "<<100./double(i)<<" %."
		       <<"   in total = "<<double(events)/double(sumtrials)*100.<<" %."<<endl;


	if (result1 < (result1+result2)*AMATOOLS::ran.Get()) Rotate(p);
	proc->Selected()->SetMomenta(p);
	for (int i=0;i<nin+nout;i++) {
	  msg.Debugging()<<"  "<<proc->Selected()->Flavs()[i]<<" : "<<p[i]<<endl; 
	} 
	return 1;
      }
    }
  }
  sumtrials += maxtrials;


  msg.Debugging()<<"Phase_Space_Handler::OneEvent() : "
		 <<" too many trials for "<<proc->Selected()->Name()<<endl
		 <<"   Efficiency = "<<double(events)/double(sumtrials)*100.<<" %."<<endl;
  return 0;
}

double Phase_Space_Handler::WeightedEvent()
{
  return 0.;
} 

void Phase_Space_Handler::AddPoint(const double value) { 
  proc->AddPoint(value); 
}




/* ----------------------------------------------------------------------

   One test point to check the amplitudes

   ---------------------------------------------------------------------- */

void Phase_Space_Handler::TestPoint(Vec4D * _p)
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
  int nran = ran.WriteOutStatus(filename);
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

  ran.ReadInStatus(filename,0);

  return okay;
}

void Phase_Space_Handler::Rotate(Vec4D * _p) 
{
  for (int i=0;i<nin+nout;i++) _p[i] = Vec4D(_p[i][0],(-1.)*Vec3D(_p[i]));
}




bool Phase_Space_Handler::CreateIntegrators()
{
  msg.Debugging()<<"In Phase_Space_Handler::CreateIntegrators"<<endl;

  if (nin==1) int_type = 0;
  if (bh) {
    if ((nin==2) && (bh && bh->On()>0) ) {
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

  if (ih) {
    if ((nin==2) && (ih && ih->On()>0)) {
      if (!(MakeISRChannels())) {
	msg.Error()<<"Error in Phase_Space_Handler::CreateIntegrators !"<<endl
		   <<"   did not construct any isr channels !"<<endl;
      }
      if (isrchannels) 
	msg.Debugging()<<"  ("<<isrchannels->Name()<<","<<isrchannels->Number()<<";";
    }
    else {
      msg.Debugging()<<" no ISR needed           : "<<ih->Name()
		     <<" for "<<nin<<" incoming particles."<<endl;
    }
  }

  msg.Debugging()<<" "<<fsrchannels->Name()<<","<<fsrchannels->Number()<<")"<<endl
		 <<" integration mode = "<<int_type<<endl;
  
  if (int_type < 3 || int_type == 5 && (fsrchannels!=0)) fsrchannels->DropAllChannels();
  
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
  case 5:
    fsrchannels->Add(new RamboKK(nin,nout,psflavs));
      break;    
  default:
    msg.Error()<<"Wrong phasespace integration switch ! Using RAMBO as default."<<endl;
    fsrchannels->Add(new Rambo(nin,nout,psflavs));
  }  
    
  msg.Debugging()<<"Initialized Phase_Space_Integrator (";
  if (beamchannels) msg.Debugging()<<beamchannels->Name()<<","<<beamchannels->Number()<<";";
  if (isrchannels)  msg.Debugging()<<isrchannels->Name()<<","<<isrchannels->Number()<<";";
  if (fsrchannels)  msg.Debugging()<<fsrchannels->Name()<<","<<fsrchannels->Number()<<")"<<endl;
  
  return 1;
}


void Phase_Space_Handler::DropRedundantChannels()
{
  fsrchannels->Reset();
  int number = fsrchannels->Number();
  
  msg.Debugging()<<"In Phase_Space_Handler::DropRedundantChannels"
		 <<"("<<fsrchannels->Name()<<")."<<endl
		 <<"    Start with "<<number<<" added channel(s)."<<endl;

  if (number<2) return;
  Vec4D** perm_vec = new Vec4D*[number]; 
  
  for (short int i=0;i<number;i++) perm_vec[i] = new Vec4D[nin+nout+1];
  
  // Create Momenta
  int rannum   = 1 + 2 + 3*(nout-2);
  double * rans = new double[rannum];
  for (short int i=0;i<rannum;i++) rans[i] = ran.Get();  
  // Init markers for deletion and results to compare.
  int    * marker = new int[number];  
  double * res    = new double[number];
  for (short int i=0;i<number;i++) { marker[i] = 0;res[i] = 0.; }

  for (short int i=0;i<number;i++) {
    perm_vec[i][0] = Vec4D(rpa.gen.Ecms()/2.,0.,0.,rpa.gen.Ecms()/2.);
    perm_vec[i][1] = Vec4D(rpa.gen.Ecms()/2.,0.,0.,-rpa.gen.Ecms()/2.); 
    msg.Debugging()<<"==== "<<i<<" : ";
    msg.Debugging()<<(fsrchannels->Channel(i))->Name()<<"====================="<<endl;
    fsrchannels->GeneratePoint(i,perm_vec[i],proc->Cuts(),rans);
    // for (short int j=0;j<nin+nout;j++) msg.Debugging()<<j<<"th : "<<perm_vec[i][j]<<endl;
    fsrchannels->GenerateWeight(i,perm_vec[i],proc->Cuts());
    res[i] = fsrchannels->Weight();
    if (res[i]==0.) {
      msg.Debugging()<<"  "<<(fsrchannels->Channel(i))->Name()<<" produced a zero weight."<<endl;
      marker[i] = 1;
    }
  }
  delete[] rans;

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


bool Phase_Space_Handler::Compare(Vec4D* _p1,Vec4D* _p2)
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


bool Phase_Space_Handler::MakeBeamChannels()
{
  if (beam_params.size()>0) return CreateBeamChannels();

  double deltay[2];
  deltay[0] = log(bh->Upper1());
  deltay[1] = log(bh->Upper2());

  Channel_Info ci;
  // default : Beamstrahlung
  if ((psflavs[0].IsLepton()) && (psflavs[1].IsLepton())) {
    ci.type = 0;
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(bh->Exponent(1));
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    beam_params.push_back(ci);
    ci.parameters.clear();
  }
  // Laser Backscattering spectrum
  if ((psflavs[0].IsPhoton()) || (psflavs[1].IsPhoton())) {
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
  int    type;
  double mass,width;

  for (int i=0;i<fsrchannels->Number();i++) {
    type = 0; mass = width = 0.;
    if (proc) fsrchannels->ISRInfo(i,type,mass,width);
    
    msg.Debugging()<<i<<" : "<<type<<"/"<<mass<<"/"<<width<<endl;
    if (AMATOOLS::IsZero(mass) || AMATOOLS::IsZero(width)) continue;
    if ((type == 0) || (type == 3))                        continue;

    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((psflavs[0].IsLepton()) || (psflavs[1].IsLepton())) (ci.parameters).push_back(1.);
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
  return CreateBeamChannels();
}


bool Phase_Space_Handler::MakeISRChannels()
{
  if (isr_params.size()>0) return CreateISRChannels();

  Channel_Info ci;

  double deltay[2];  
  deltay[0] = log(ih->Upper1());
  deltay[1] = log(ih->Upper2());
  msg.Out()<<"*** DeltaY1 / 2 = "<<deltay[0]<<" / "<<deltay[1]<<endl
	   <<"*** Exponents :   "<<ih->Exponent(0)<<" / "<<ih->Exponent(1)<<endl;


  if ((psflavs[0].IsLepton()) || (psflavs[1].IsLepton())) {
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

    ci.type = 0;
    (ci.parameters).push_back(2.);
    (ci.parameters).push_back(0.5);
    (ci.parameters).push_back(deltay[0]);
    (ci.parameters).push_back(deltay[1]);
    isr_params.push_back(ci);
    ci.parameters.clear();
  }

  bool   addit;
  int    type;
  double mass,width;
    
  for (int i=0;i<fsrchannels->Number();i++) {
    type = 0; mass = width = 0.;
    fsrchannels->ISRInfo(i,type,mass,width);
    
    msg.Debugging()<<i<<" : "<<type<<"/"<<mass<<"/"<<width<<endl;
    if (AMATOOLS::IsZero(mass) || AMATOOLS::IsZero(width)) continue;
    if ((type == 0) || (type == 3))                        continue;

    ci.type = type;
    (ci.parameters).push_back(mass);
    (ci.parameters).push_back(width);
    if ((psflavs[0].IsLepton()) || (psflavs[1].IsLepton())) (ci.parameters).push_back(1.);
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
  return CreateISRChannels();
}

bool Phase_Space_Handler::CreateBeamChannels()
{
  if (beam_params.size() < 1) return 0;

  Single_Channel * channel;   

  for (int i=0;i<beam_params.size();i++) {
    if ((beam_params[i]).type==0) {
      channel = new SimplePoleCentral(beam_params[i].parameters[0],
				      beam_params[i].parameters[2],
				      beam_params[i].parameters[3]);
      beamchannels->Add(channel);
      channel = new SimplePoleForward(beam_params[i].parameters[0],
				      beam_params[i].parameters[1],
				      beam_params[i].parameters[2],
				      beam_params[i].parameters[3]);
      beamchannels->Add(channel);
      channel = new SimplePoleBackward(beam_params[i].parameters[0],
				       beam_params[i].parameters[1],
				       beam_params[i].parameters[2],
				       beam_params[i].parameters[3]);
      beamchannels->Add(channel);
    } 
    if ((beam_params[i]).type==1) {
      channel = new ResonanceCentral(beam_params[i].parameters[0],
				     beam_params[i].parameters[1],
				     beam_params[i].parameters[3],
				     beam_params[i].parameters[4]);
      beamchannels->Add(channel);
      channel = new ResonanceForward(beam_params[i].parameters[0],
				     beam_params[i].parameters[1],
				     beam_params[i].parameters[2],
				     beam_params[i].parameters[3],
				     beam_params[i].parameters[4]);
      beamchannels->Add(channel);
      channel = new ResonanceBackward(beam_params[i].parameters[0],
				      beam_params[i].parameters[1],
				      beam_params[i].parameters[2],
				      beam_params[i].parameters[3],
				      beam_params[i].parameters[4]);
      beamchannels->Add(channel);
    }
    if ((beam_params[i]).type==3) {
      if ((psflavs[0].IsPhoton()) || (psflavs[1].IsPhoton())) {
	  channel = new LBSComptonPeakCentral(beam_params[i].parameters[0],
					      beam_params[i].parameters[1],
					      beam_params[i].parameters[2],
					      beam_params[i].parameters[3]);
	  beamchannels->Add(channel);
	  channel = new LBSComptonPeakForward(beam_params[i].parameters[0],
					      beam_params[i].parameters[1],
					      beam_params[i].parameters[2],
					      beam_params[i].parameters[3],
					      beam_params[i].parameters[4]);
	  beamchannels->Add(channel);
	  channel = new LBSComptonPeakBackward(beam_params[i].parameters[0],
					       beam_params[i].parameters[1],
					       beam_params[i].parameters[2],
					       beam_params[i].parameters[3],
					       beam_params[i].parameters[4]);
	  beamchannels->Add(channel);
      }
    }
  }
  return 1;
}


bool Phase_Space_Handler::CreateISRChannels()
{
  if (isr_params.size() < 1) return 0;

  Single_Channel * channel;   

  int length = isr_params.size();

  for (int i=0;i<length;i++) {
    if ((isr_params[i]).type==0) {
      // Maybe set also the exponent of the y - integral ???
      channel = new SimplePoleCentral(isr_params[i].parameters[0],
				      isr_params[i].parameters[2],
				      isr_params[i].parameters[3]);
      isrchannels->Add(channel);
      channel = new SimplePoleForward(isr_params[i].parameters[0],
				      isr_params[i].parameters[1],
				      isr_params[i].parameters[2],
				      isr_params[i].parameters[3]);
      isrchannels->Add(channel);
      channel = new SimplePoleBackward(isr_params[i].parameters[0],
				       isr_params[i].parameters[1],
				       isr_params[i].parameters[2],
				       isr_params[i].parameters[3]);
      isrchannels->Add(channel);
    } 
    if ((isr_params[i]).type==1) {
      channel = new ResonanceCentral(isr_params[i].parameters[0],
				     isr_params[i].parameters[1],
				     isr_params[i].parameters[3],
				     isr_params[i].parameters[4]);
      isrchannels->Add(channel);
      channel = new ResonanceForward(isr_params[i].parameters[0],
				     isr_params[i].parameters[1],
				     isr_params[i].parameters[2],
				     isr_params[i].parameters[3],
				     isr_params[i].parameters[4]);
      isrchannels->Add(channel);
      channel = new ResonanceBackward(isr_params[i].parameters[0],
				      isr_params[i].parameters[1],
				      isr_params[i].parameters[2],
				      isr_params[i].parameters[3],
				      isr_params[i].parameters[4]);
      isrchannels->Add(channel);
    }
    if (isr_params[i].type==3) {
      channel = new LLCentral(isr_params[i].parameters[0],
			      isr_params[i].parameters[2],
			      isr_params[i].parameters[3]);
      isrchannels->Add(channel);
      channel = new LLForward(isr_params[i].parameters[0],
			      isr_params[i].parameters[1],
			      isr_params[i].parameters[2],
			      isr_params[i].parameters[3]);
      isrchannels->Add(channel);
      channel = new LLBackward(isr_params[i].parameters[0],
			       isr_params[i].parameters[1],
			       isr_params[i].parameters[2],
			       isr_params[i].parameters[3]);
      isrchannels->Add(channel);
    }
  }
  return 1;
}

