#include "Process_Group.H"
#include "Single_Process.H"

#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

#include "Running_AlphaS.H"
#include "Combined_Selector.H"
#include "Primitive_Observable_Base.H"

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

using namespace AMEGIC;
using namespace PHASIC;
using namespace AORGTOOLS;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Process_Group::Process_Group() { 
  name  = "Empty_Group"; 
  fl    = flin  = flout = 0;
  pl    = plin  = plout = 0;
  isr   = 0; 
  beam  = 0;
  sel   = 0; 
  cuts  = 0;
  ps    = 0;
  moms  = 0;


  analyse  = 0;
  analysis = 0;

  n        = 0;
  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = 0.;
  lastlumi = 1.;

  atoms = 1;
}



Process_Group::Process_Group(int _nin,int _nout,Flavour *& _fl,
			     ISR::ISR_Handler * _isr,BEAM::Beam_Handler * _beam,
			     APHYTOOLS::Selector_Data * _seldata,
			     int _gen_str,int _kfactorscheme, int _scalescheme, 
			     Pol_Info * _pl) 
{
  msg.Tracking()<<"In Process_Group("<<_nin<<","<<_nout<<") : "<<gen_str<<endl;
  nin = _nin; nout = _nout; isr = _isr; beam = _beam;
  kfactorscheme = _kfactorscheme;
  scalescheme   = _scalescheme;


  flin    = new Flavour[nin];
  flout   = new Flavour[nout];  
  plin    = new Pol_Info[nin];
  plout   = new Pol_Info[nout]; 

  for (short int i=0;i<nin;i++) {
    flin[i]  = _fl[i];
    if (_pl!=0) plin[i] = _pl[i];
           else plin[i] = Pol_Info(flin[i]); 

  }
  for (short int i=0;i<nout;i++) { 
    flout[i] = _fl[i+nin]; 
    if (_pl!=0) plout[i] = _pl[i+nin];
           else plout[i] = Pol_Info(flout[i]); 
  }

  string stan,oli;
  GenerateNames(nin,flin,plin,nout,flout,plout,name,stan,oli);

  ConstructProcesses(_seldata);
  GroupProcesses();

  Initialize(_seldata);
}



Process_Group::~Process_Group()
{
  for(int i=procs.size();i>0;i--) {
    if (procs[i-1]) delete procs[i-1];
  }

  if (fl)       { delete [] fl;    fl       = 0; }
  if (flin)     { delete [] flin;  flin     = 0; }
  if (flout)    { delete [] flout; flout    = 0; }
  if (pl)       { delete [] pl;    pl       = 0; }
  if (plin)     { delete [] plin;  plin     = 0; }
  if (plout)    { delete [] plout; plout    = 0; }
  if (moms)     { delete [] moms;  moms     = 0; }
  if (sel)      { delete sel;      sel      = 0; }
  if (cuts)     { delete cuts;     cuts     = 0; }
  if (ps)       { delete ps;       ps       = 0; }
  if (analysis) { delete analysis; analysis = 0; }
}




/*----------------------------------------------------------------------------------
  
  Management of the processes included in the Process_Group

  ----------------------------------------------------------------------------------*/

void Process_Group::ConstructProcesses(APHYTOOLS::Selector_Data * _seldata) {
  int * flindex;
  flindex = new int[nin+nout];
  for (int i=0;i<nin+nout;i++) flindex[i] = 0;
  Flavour  * _flin, * _flout, * _fl;
  Pol_Info * _plin, *_plout, *_pl;
  _flin   = new Flavour[nin];
  _plin   = new Pol_Info[nin];
  _flout  = new Flavour[nout];
  _plout  = new Pol_Info[nout];
  _fl     = new Flavour[nin+nout];
  _pl     = new Pol_Info[nin+nout];

  string _name,_stan,_oli;
  bool flag = 1;
  bool take,overflow;
  for (;;) {
    if (!flag) break;
    for (int i=0;i<nin;i++) {
      if (flin[i].Size() != 1) _flin[i] = flin[i][flindex[i]]; 
                          else _flin[i] = flin[i];
      _plin[i]                          = plin[i];
    }
    for (int i=0;i<nout;i++) {
      if (flout[i].Size() != 1) _flout[i] = flout[i][flindex[nin+i]]; 
                           else _flout[i] = flout[i];
      _plout[i]                           = plout[i];
    }
    GenerateNames(nin,_flin,_plin,nout,_flout,_plout,_name,_stan,_oli);
    take = 1;
    for (int k=0;k<procs.size();k++) {
      if (_name == procs[k]->Name()) { take = 0; break; }
    }
    if (take) {
      if (CheckExternalFlavours(nin,_flin,nout,_flout)) {
	for (int i=0;i<nin;i++)  _fl[i]     = _flin[i];  
	for (int i=0;i<nout;i++) _fl[i+nin] = _flout[i]; 
	Add(new Single_Process(nin,nout,_fl,isr,beam,_seldata,2,
			       kfactorscheme,scalescheme,pl));
      }
    }
    overflow = 0;
    for (int i=1;i<nout+1;i++) {
      if (flout[nout-i].Size()-1 > flindex[nin+nout-i]) {
	flindex[nin+nout-i] = flindex[nin+nout-i]+1; 
	break;
      }
      else {
	if (i==nout) overflow = 1;
	flindex[nin+nout-i]   = 0;
      }
    }
    if (overflow) {
      overflow = 0;
      for (int i=1;i<nin+1;i++) {
	if (flin[nin-i].Size()-1 > flindex[nin-i]) {
	  flindex[nin-i] = flindex[nin-i]+1; 
	  break;
	}
	else {
	  if (i==nin) flag = 0;
	  flindex[nin-i]   = 0;
	}
      }
    }
  }
  delete [] _fl;
  delete [] _flin;
  delete [] _flout;
  delete [] _pl;
  delete [] _plin;
  delete [] _plout;
}


void Process_Group::GroupProcesses() {
  // First : Check for identical masses.
  double * massin  = new double[nin];
  double * massout = new double[nout];
  for (int i=0;i<nin;i++)  massin[i]  = flin[i].mass();
  for (int i=0;i<nout;i++) massout[i] = flout[i].mass();

  bool massok = 1;
  for (int i=0;i<procs.size();i++) {
    for (int j=0;j<procs[i]->Nin();j++) {
      if (!(AMATOOLS::IsEqual(massin[j],(procs[i]->Flavs()[j]).mass()))) {
	massok = 0; break;
      }
    }
    if (!massok) break;
    for (int j=0;j<procs[i]->Nout();j++) {
      if (!(AMATOOLS::IsEqual(massout[j],(procs[i]->Flavs()[j+procs[i]->Nin()]).mass()))) {
	massok = 0; break;
      }
    }
    if (!massok) break;
  }
  if (massok) msg.Debugging()<<"All processes in "<<name<<" have equal masses."<<endl;
  else {
    msg.Error()<<"Error in Process_Group : "<<name<<endl
	       <<"   Processes do not have equal masses. Abort."<<endl;
    abort();
  }
  delete [] massin;
  delete [] massout;

  std::vector<Process_Base *> singleprocs = procs;
  while (procs.size()>0) procs.pop_back();
  msg.Debugging()<<"Emptied process list             : "<<procs.size()<<endl
		 <<"Created list of Single_Processes : "<<singleprocs.size()<<endl;
  
  Process_Base  * sproc;
  Process_Group * group;
  string         help;
  Flavour        flav1,flav2;
  for (int i=0;i<singleprocs.size();i++) {
    sproc      = singleprocs[i];
    help       = string("SG_");
    if (sproc->Nin()==2) {
      flav1    = sproc->Flavs()[0]; 
      flav2    = sproc->Flavs()[1]; 
      if ( (flav1.isvector()) && (flav2.isvector()) ) {
	if ( (flav1.bar() == flav1) && (flav2.bar() == flav2) ) {
	  if (flav1 == flav2) help += string("V V -> ");
	                 else help += string("V V' -> ");
	}
	else if (flav1.bar() == flav1) { 
	  if (flav2.charge() > 0) help += string("V0 V+ -> "); 
	                     else help += string("V0 V- -> "); 
	}
	else if (flav2.bar() == flav2) { 
	  if (flav1.charge() > 0) help += string("V+ V0 -> "); 
	                     else help += string("V- V0 -> "); 
	}
	else {
	  if (flav1.charge() > 0) help += string("V+ "); 
	                     else help += string("V- "); 
	  if (flav2.charge() > 0) help += string("V+ -> "); 
	                     else help += string("V- -> "); 
	}
      }
      else if ( (flav1.isfermion()) && (flav2.isvector()) ) {
	if (flav1.isanti())          help += string("fb ");
	                        else help += string("f ");
	if (flav2==flav2.bar())      help += string("V -> ");
	else if (flav2.charge() > 0) help += string("V+ -> "); 
	                        else help += string("V- -> "); 
      }
      else if ( (flav1.isvector()) && (flav2.isfermion()) ) {
	if (flav1==flav1.bar())      help += string("V ");
	else if (flav1.charge() > 0) help += string("V+ "); 
	                        else help += string("V- "); 
	if (flav2.isanti())          help += string("fb -> ");
	                        else help += string("f -> ");
      }
      else if ( (flav1.isfermion()) && (flav2.isfermion()) ) {
	if ( (flav1.isanti()) && (flav2==flav1) )             help += string ("fb fb -> ");
	else if ( !(flav1.isanti()) && (flav2==flav1) )       help += string ("f f -> ");
	else if ( !(flav1.isanti()) && (flav2==flav1.bar()) ) help += string ("f fb -> ");
	else if ( (flav1.isanti()) && (flav2.isanti()) )      help += string ("fb fb' -> ");
	else if ( !(flav1.isanti()) && (flav2.isanti()) )     help += string ("fb f' -> ");
	else if ( (flav1.isanti()) && !(flav2.isanti()) )     help += string ("f fb' -> ");
	else if ( !(flav1.isanti()) && !(flav2.isanti()) )    help += string ("f f' -> ");
      }
    }
    else {
      msg.Error()<<"Error in Process_Group::GroupProcesses :"<<endl
		 <<"Right now process groups with one incoming particle are NOT enabled."<<endl;
      abort();
    }
    int scalars = 0,fermions = 0,vectors = 0;
    for (int j=0;j<sproc->Nout();j++) {
      if ((sproc->Flavs()[sproc->Nin()+j]).isscalar())  scalars++;
      if ((sproc->Flavs()[sproc->Nin()+j]).isfermion()) fermions++;
      if ((sproc->Flavs()[sproc->Nin()+j]).isvector())  vectors++;
    }
    char numb[20];
    sprintf(numb,"%i",scalars);
    help += string(numb) + string("S ");
    sprintf(numb,"%i",fermions);
    help += string(numb) + string("F ");
    sprintf(numb,"%i",vectors);
    help += string(numb) + string("V");


    bool found = 0;
    for (int j=0;j<procs.size();j++) {
      if (procs[j]->Name() == help) {
	procs[j]->Add(sproc);
	found = 1;
	break;
      }
    }
    if (!found) {
      group = new Process_Group();
      group->SetName(help);
      group->SetAtoms(0);
      group->SetBeam(beam);
      group->SetISR(isr);
      group->Add(sproc);
      procs.push_back(group);
    }
  }

  for (int i=0;i<procs.size();i++) {
    msg.Tracking()<<"Process_Group "<<procs[i]->Name()<<" : "<<procs[i]->Size()<<endl;
    for (int j=0;j<procs[i]->Size();j++) 
      msg.Tracking()<<"    "<<((*procs[i])[j])->Name()<<endl;
    msg.Tracking()<<"--------------------------------------------------"<<endl;
  }
}



void Process_Group::Add(Process_Base * _proc) 
{
  if (procs.size()==0) {
    nin     = _proc->Nin();
    nout    = _proc->Nout();
    nvec    = _proc->Nvec();
    nstrong = _proc->NStrong();
    neweak  = _proc->NEWeak();
  }
  else {
    if (_proc->Nvec() > nvec) nvec = _proc->Nvec();
  }
  if ( (nin != _proc->Nin()) || (nout != _proc->Nout())) {
    msg.Error()<<"Error : Cannot add process "<<_proc->Name()
	       <<" to group "<<name<<" ! "<<endl
	       <<"   Inconsistent number of external legs."<<endl
	       <<"  Before : ("<<nin<<" -> "<<nout<<" )"<<endl
	       <<"  Now    : ("<<_proc->Nin()<<" -> "<<_proc->Nout()<<" )"<<endl;
    return;
  }
  msg.Tracking()<<"Add process "<<_proc->Name()<<" to group "<<name<<" ! "<<endl; 
  procs.push_back(_proc);
}



void Process_Group::SelectOne()
{
  DeSelect();
  if (totalxs==0) selected = procs[int(Ran.get()*procs.size())];
  else {
    double disc = max * Ran.get();
    for (int i=0;i<procs.size();i++) {
      disc -= procs[i]->Max();
      if (disc<0.) {
	selected = procs[i];
	msg.Tracking()<<"Selected Process(_Group) : "<<selected->Name()<<endl;	
	selected->SelectOne();
	return;
      }
    }
    if (disc>0.) { 
      msg.Error()<<"Error in Process_Group::SelectOne() : ";
      msg.Error()<<"Total xsec, max = "<<totalxs<<", "<<max<<endl;
      return;
    }
  }
}



void Process_Group::DeSelect() {
  selected = 0;
  for (int i=0;i<procs.size();i++) procs[i]->DeSelect();
}



Process_Base * Process_Group::Selected() { 
  if (selected==this) return this;
  return selected->Selected(); 
}    



void Process_Group::Empty() {
  for (int i=0;i<procs.size();i++) procs[i]->Empty();
}



void Process_Group::SetResDir(std::string _resdir) {
  resdir = _resdir;
  for (int i=0;i<procs.size();i++) procs[i]->SetResDir(resdir);
}


void Process_Group::SetScale(double _scale)
{
  scale = _scale;
  for (int i=0;i<procs.size();i++) procs[i]->SetScale(scale); 
} 



void Process_Group::SetTables(bool _tables)
{
  tables = _tables;
  for (int i=0;i<procs.size();i++) procs[i]->SetTables(tables);
} 


void Process_Group::SetTotalXS(int tables)  { 
  if (analysis) analysis->FinishAnalysis(resdir+string("/Tab")+name,tables);
  totalxs  = totalsum/n; 
  totalerr = sqrt( (totalsumsqr/n - 
		    (AMATOOLS::sqr(totalsum)-totalsumsqr)/n/(n-1) )  / n); 
  if (sel) sel->Output();


  max = 0.;
  for (int i=0;i<procs.size();i++) {
    procs[i]->SetTotalXS(tables);
    max += procs[i]->Max();
  }
  msg.Events()<<"-----------------------------------------------------------------------"<<endl
	      <<"Total XS for "<<name<<"("<<procs.size()<<") : "<<totalxs*rpa.Picobarn()<<" pb"
	      <<" +/- "<<totalerr/totalxs*100.<<"%,"<<endl
	      <<"      max = "<<max<<endl;
}

void Process_Group::SetAtoms(bool _atoms) { atoms = _atoms; }

/*----------------------------------------------------------------------------------
  
  Initialization of the processes

  ----------------------------------------------------------------------------------*/

void Process_Group::Initialize(APHYTOOLS::Selector_Data * _seldata) {
  fl   = new Flavour[nvec];
  pl   = new Pol_Info[nvec];
  b    = new int[nvec];
  for (short int i=0;i<nin;i++) { 
    fl[i] = flin[i]; 
    pl[i] = plin[i]; 
    b[i]  = -1; 
  }
  for (short int i=nin;i<nin+nout;i++)  { 
    fl[i] = flout[i-nin]; 
    pl[i] = plout[i-nin]; 
    b[i]  = 1; 
  } 
  for (short int i=nin+nout;i<nvec;i++) { 
    fl[i] = Flavour(kf::pol); 
    b[i]  = 1; 
  }

  atoms     = 0;
  selected  = 0;
  analysis  = 0;
  analyse   = 0;

  InitCuts();

  if (_seldata) sel = new Combined_Selector(nin,nout,fl,_seldata);
  else {
    msg.Error()<<"Potential Error in Single_Process "<<name<<endl
	       <<"   No selection cuts specified. Init No_Selector !"<<endl;
    sel = new No_Selector();
  }

  moms = 0;

  totalxs  = totalerr = totalsum = totalsumsqr = 0.;
  last     = lastdxs  = max      = 0.;
  lastlumi = 1.;
}

int Process_Group::InitAmplitude(Topology * top,vec4d *& testmoms,
				  vector<double> & results,vector<Single_Process *> & links)
{
  int okay = 1;
  vector <string> deletethem;
  for (int i=0;i<procs.size();i++) {
    msg.Tracking()<<"========================================================="<<endl
		  <<"========================================================="<<endl
		  <<"Process_Group::InitAmplitude for "<<procs[i]->Name()<<endl;
    if (atoms) { delete [] testmoms; testmoms = 0; }
    switch (procs[i]->InitAmplitude(top,testmoms,results,links)) {
    case -2 : 
      msg.Error()<<"Error in creation of amplitude "<<procs[i]->Name()<<endl;
      deletethem.push_back(procs[i]->Name());
      break;
    case -1 : 
      msg.Events()<<"No diagrams or amplitudes for "<<procs[i]->Name()<<endl
		  <<"   delete it."<<endl;
      deletethem.push_back(procs[i]->Name());
      break;
    case 0 :
      okay = 0;      
      break;
    default :
      break;
    }
  }

  bool flag;
  for (int i=0;i<procs.size();i++) {
    flag = 0;
    if (procs[i]->Size() == 0) flag = 1;
    if (!flag) {
      for (int j=0;j<deletethem.size();j++) {
	if (procs[i]->Name() == deletethem[j]) flag = 1;
      }
    }
    if (flag) {
      delete procs[i];
      for (int j=i;j<procs.size()-1;j++) procs[j] = procs[j+1];
      procs.pop_back();
    }
  }

  msg.Tracking()<<"Process_Group::Initialize Amplitude for "<<name;
  if (okay) msg.Tracking()<<" successful."<<endl;
       else msg.Tracking()<<" failed."<<endl;
  return okay;
}



bool Process_Group::SetUpIntegrator()
{
  bool okay = 1;
  if (atoms) {
    for (int i=0;i<procs.size();i++) {
      if (procs[i]->Partner() == NULL) {
	if (!(procs[i]->SetUpIntegrator())) okay = 0;
      }
    }
    return okay;
  }
  
  sel->BuildCuts(cuts);
  if (nin==2) {
    if ( (fl[0].mass() != rpa.gen.Beam1().mass()) ||
	 (fl[1].mass() != rpa.gen.Beam2().mass()) ) isr->SetPartonMasses(fl);
  }
  ps  = new Phase_Space_Handler(this,isr,beam);
  ps->CollectChannels();
  if (!ps->CreateIntegrators()) return 0;
  for (int i=0;i<procs.size();i++) procs[i]->Empty();
  return 1;
}



void Process_Group::InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *> _obs) {
  analysis = new APHYTOOLS::Primitive_Analysis(this->Name());//check thsi
  for (int i=0;i<_obs.size();i++) {
    analysis->AddObservable(_obs[i]->GetCopy());
  }
  for (int i=0;i<procs.size();i++) procs[i]->InitAnalysis(_obs);
  analyse  = 1;
}



void Process_Group::InitCuts() {
  cuts = new Cut_Data();
  cuts->Init(nin+nout,fl);
}

/*----------------------------------------------------------------------------------
  
  Evaluation of the processes included in All_Processes

  ----------------------------------------------------------------------------------*/

bool Process_Group::CalculateTotalXSec()
{
  msg.Tracking()<<"Process_Group::CalculateTotalXSec()"<<endl;
  if (atoms) {
    bool okay = 1;
    for (int i=0;i<procs.size();i++) {
      msg.Tracking()<<"Process_Group::CalculateTotalXSec for "<<procs[i]->Name()<<endl;
      if (!(procs[i]->CalculateTotalXSec())) okay = 0;
    }
    return okay;
  }
  else {
    if (nin==2) {
      if ( (fl[0].mass() != rpa.gen.Beam1().mass()) ||
	   (fl[1].mass() != rpa.gen.Beam2().mass()) ) isr->SetPartonMasses(fl);
    }
    sel->BuildCuts(cuts);
    tables  = 0;
    totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn(); 
    if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl
		 <<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
    }
    SetTotalXS(0);
    if (totalxs>0.) return 1;
  }
}


void  Process_Group::RescaleXSec(double fac) {
  Process_Base::RescaleXSec(fac);
  for (int i=0;i<procs.size();i++) {
    procs[i]->RescaleXSec(fac);
  }
}

bool Process_Group::LookUpXSec(double ycut,bool calc,string obs) {
  msg.Events()<<"Process_Group::LookUpXSec() for "<<name<<" in "<<resdir<<endl;
  bool okay = 1;
  if (atoms) {
    for (int i=0;i<procs.size();i++) {
      if (!(procs[i]->LookUpXSec(ycut,1,obs))) okay = 0;
    }
    return okay;
  }
  else {
    for (int i=0;i<procs.size();i++) {
      if (!(procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
    }
    if (okay) {                    
      totalxs = 0; max = 0;
      for (int i=0;i<procs.size();i++) {
	totalxs += procs[i]->Total();
	max     += procs[i]->Max();
      }
      msg.Events()<<name<<" : Set total xsec and max at ycut = "<<ycut
		  <<" : "<<endl<<"   "<<totalxs<<" / "<<max<<endl;

      ps->ReadIn(resdir+string("/MC_")+name);
      if (ps->BeamIntegrator() != 0) ps->BeamIntegrator()->Print();
      if (ps->ISRIntegrator() != 0)  ps->ISRIntegrator()->Print();
      if (ps->FSRIntegrator() != 0)  ps->FSRIntegrator()->Print();
      return 1;
    }
    if (calc) {
      if (!(PrepareXSecTables())) return 0;
      okay = 1;
      for (int i=0;i<procs.size();i++) {
	if (!(procs[i]->LookUpXSec(ycut,0,obs))) okay = 0;
      }
      if (okay) {                    
	totalxs = 0; max = 0;
	for (int i=0;i<procs.size();i++) {
	  totalxs += procs[i]->Total();
	  max     += procs[i]->Max();
	}
	msg.Out()<<name<<" : Set total xsec and max at ycut = "<<ycut
		    <<" : "<<endl<<"   "<<totalxs<<" / "<<max<<endl;
	return 1;
      }
    }
    return 0;
  }
}



bool Process_Group::PrepareXSecTables()
{
  msg.Tracking()<<"Process_Group::PrepareXSecTables()"<<endl;

  if (atoms) {
    bool okay = 1;
    for (int i=0;i<procs.size();i++) {
      msg.Tracking()<<"Process_Group::PrepareXSecTables() for "<<procs[i]->Name()<<endl;
      if (!(procs[i]->PrepareXSecTables())) okay = 0;
    }
    return okay;
  }
  else {
    if (nin==2) {
      if ( (fl[0].mass() != rpa.gen.Beam1().mass()) ||
	   (fl[1].mass() != rpa.gen.Beam2().mass()) ) isr->SetPartonMasses(fl);
    }
    sel->BuildCuts(cuts);
    totalxs = ps->Integrate()/AORGTOOLS::rpa.Picobarn(); 
    if (!(AMATOOLS::IsZero((n*totalxs-totalsum)/(n*totalxs+totalsum)))) {
      msg.Error()<<"Result of PS-Integrator and internal summation do not coincide!"<<endl;
      msg.Error()<<"  "<<name<<" : "<<totalxs<<" vs. "<<totalsum/n<<endl;
    }
    SetTotalXS(1);
    ps->WriteOut(resdir+string("/MC_")+name);
    if (totalxs>0.) return 1;
  }
}



void Process_Group::AddPoint(const double value) 
{
  n++;
  totalsum    += value;
  totalsumsqr += value*value;
  if (value>max) max = value;

  if (analyse) analysis->DoAnalysis(value*rpa.Picobarn());
  for (int i=0;i<procs.size();i++) {
    if (dabs(last)>0.) {
      procs[i]->AddPoint(value*procs[i]->Last()/last);
    }
    else {
      procs[i]->AddPoint(0.);
    }  
  }
}






double Process_Group::Differential(vec4d * p)
{
  last = 0;
  for (int i=0;i<procs.size();i++) {
    last += procs[i]->DSigma(p,1);
  }
  //ControlOutput(p);

  if ((!(last<=0)) && (!(last>0))) {
    msg.Error()<<"---- Process_Group::Differential -------------------"<<endl;
    for (int i=0;i<nin+nout;i++) 
      msg.Error()<<i<<" th Momentum "<<p[i]<<endl;
    PrintDifferential();
  }
  return last;
}



double Process_Group::Differential2()
{
  if (isr->On()==0) return 0.;
  double tmp = 0.;
  for (int i=0;i<procs.size();i++) tmp += procs[i]->DSigma2();

  if ((!(tmp<=0)) && (!(tmp>0))) {
    msg.Error()<<"---- Process_Group::Differential -------------------"<<endl;
    PrintDifferential();
  }
  last += tmp;
  return tmp;
}



double Process_Group::DSigma(vec4d * p,bool lookup)
{
  last = 0;
  for (int i=0;i<procs.size();i++) last += procs[i]->DSigma(p,lookup);
  // ControlOutput(p);
  return last;
}



double Process_Group::DSigma2()
{
  double tmp = 0.;
  for (int i=0;i<procs.size();i++) tmp += procs[i]->DSigma2();
  last += tmp;
  return tmp;
}


bool Process_Group::OneEvent() {
  if (atoms) {
    SelectOne();
    return selected->OneEvent();
  }
  else return ps->OneEvent();
}

bool Process_Group::SameEvent() {
  if (atoms) {
    if (selected)
      return selected->SameEvent();
    msg.Error()<<" ERROR in bool Process_Group::SameEvent() "<<endl;
    return 0;
  }
  else return ps->SameEvent();
}

double Process_Group::WeightedEvent() {
  if (atoms) {
    SelectOne();
    return selected->WeightedEvent();
  }
  else return ps->WeightedEvent();
}




/*----------------------------------------------------------------------------------
  
  Helpers

  ----------------------------------------------------------------------------------*/

void Process_Group::PrintDifferential()
{
  if (!(rpa.gen.Tracking())) return;
  last = 0;
  msg.Out()<<"--------------------------------------------------------"<<endl;
  msg.Out()<<"--------------------------------------------------------"<<endl;
  msg.Out()<<"--------------------------------------------------------"<<endl;
  for (int i=0;i<procs.size();i++) {
    msg.Out()<<"====================================================="<<endl;
    procs[i]->PrintDifferential();
    last += procs[i]->Last();
    msg.Out()<<"====================================================="<<endl;
  }
  msg.Out()<<"--------------------------------------------------------"<<endl;
  msg.Out()<<" Total : "<<last<<endl;
  msg.Out()<<"--------------------------------------------------------"<<endl;
  msg.Out()<<"--------------------------------------------------------"<<endl;
  msg.Out()<<"--------------------------------------------------------"<<endl;
}

void Process_Group::ControlOutput(vec4d * p)
{ 
  msg.Out()<<"***************************************************************************"<<endl;
  msg.Out()<<"***************************************************************************"<<endl;
  double s   = (p[0]+p[1]).abs2();
  double t   = (p[0]-p[2]).abs2();
  double u   = (p[0]-p[3]).abs2();
  double a_s = as->AsFixed();
  msg.Out()<<"-------- Process_Group : "<<name<<" : DSigma -------------------"<<endl;
  msg.Out()<<"         scale = "<<scale<<" = "<<2.*s*t*u/(s*s+u*u+t*t)<<endl;
  msg.Out()<<"         s,t,u = "<<s<<", "<<t<<", "<<u<<" : "<<sqrt(4.*M_PI*a_s)<<endl;
  msg.Out()<<"-----------------------------------------------------------------------"<<endl;
  double g4  = sqr(4.*M_PI*a_s);
  if (name == string("gg -> gg"))
    msg.Out()<<"gg   -> gg   : "
	     <<g4 * ((3.-(t*u)/(s*s)-(s*u)/(t*t)-(t*s)/(u*u)))<<endl;
  if (name == string("qg -> qg"))
    msg.Out()<<"qg   -> qg   : "
	     <<g4 * (-4./9.*(s*s+u*u)/(s*u)+(s*s+u*u)/(t*t))<<endl;
  if (name == string("qbg -> qbg"))
    msg.Out()<<"qg   -> qg   : "
	     <<g4 * (-4./9.*(s*s+u*u)/(s*u)+(s*s+u*u)/(t*t))<<endl;
  if (name == string("qbq'b -> qbq'b"))
    msg.Out()<<"qq'  -> qq'   : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;

  // if (name == string("SG_f f' -> 0S 2F 0V"))
    msg.Out()<<"qq'  -> qq'   : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;
  if (name == string("qq'b -> qq'b"))
    msg.Out()<<"qq'b -> qq'b  : "
	     <<g4 * (4.*(s*s+u*u)/(9.*t*t))<<endl;
  // if (name == string("qq -> qq"))
    msg.Out()<<"qq   -> qq    : "
	     <<g4 * ((4./9.*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-
		       8./27.*(s*s)/(u*t))/2.)<<endl;
  if (name == string("SG_f f -> 0S 2F 0V"))
    msg.Out()<<"qq   -> qq    : "
	     <<g4 * ((4./9.*((s*s+u*u)/(t*t)+(s*s+t*t)/(u*u))-
		       8./27.*(s*s)/(u*t))/2.)<<endl;
  if (name == string("qqb -> q'q'b"))
    msg.Out()<<"qqb  -> q'q'b : "
	     <<g4 * (4.*(t*t+u*u)/(9.*s*s))<<endl;
  if (name == string("qqb -> qqb"))
    msg.Out()<<"qqb  -> qqb   : "
	     <<g4 * (4./9.*((s*s+u*u)/(t*t)+(u*u+t*t)/(s*s))-
		      8./27.*(u*u)/(s*t))<<endl;
  if (name == string("qqb -> gg"))
    msg.Out()<<"qqb  -> gg    : "
	     <<g4 * ((32./27.*(t*t+u*u)/(t*u)-
		       8./3.*(t*t+u*u)/(s*s))/2.)<<endl;
  if (name == string("gg -> qqb"))
    msg.Out()<<"gg   -> qqb   : "
	     <<g4 * (1./6.*(t*t+u*u)/(t*u)-3./8.*(t*t+u*u)/(s*s))<<endl;
  msg.Out()<<"-----------------------------------------------------------------------"<<endl;
}

double Process_Group::DSigma(double s, double t, double u) {
  // this is a dummy method for the use with XS'
  AORGTOOLS::msg.Error()<<"Error : Process_Base::Dsigma(s,t,u) called in Single_Process"<<std::endl;
}
