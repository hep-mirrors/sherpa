#include "Amegic.H"
#include "All_Processes.H"
#include "Single_Process.H"

#include "Model_Handler.H"
#include "Vertex.H"
#include "Topology.H"
#include "Pol_Info.H"

#include "Run_Parameter.H"
#include "Message.H"

#include <iomanip>


using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace BEAM;
using namespace ISR;
using namespace std;



/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Model * mo =0;

Amegic::Amegic(string _path,ISR_Handler * _isr,Beam_Handler * _beam) :
  path(_path), isr(_isr), beam(_beam)
{
  gen_str = 2; //strings + libraries
  resdir  = string(".");

  polbunches = 0;

  Model_Handler mh;
  mo = mh.GetModel();
  mo->Init();
  mo->Init_Vertex();

  p_fifo = new ofstream("fifotest.out");
}
 

Amegic::~Amegic() {
  if (mo)        { delete mo;        mo        = 0; }
  if (top)       { delete top;       top       = 0; }
  if (procs)     { delete procs;     procs     = 0; }
  if (seldata)   { delete seldata;   seldata   = 0; }
  if (beamtypes) { delete beamtypes; beamtypes = 0; }
  if (isrtypes)  { delete isrtypes;  isrtypes  = 0; }
  if (splimits)  { delete splimits;  splimits  = 0; }
  if (bunches)   { delete bunches;   bunches   = 0; }
  if (beams)     { delete beams;     beams     = 0; }
  if (partons)   { delete partons;   partons   = 0; }
  if (beam)      { delete beam;      beam      = 0; }
  if (isr)       { delete isr;       isr       = 0; }
  if (p_fifo)    { delete p_fifo;    p_fifo       = 0; }
  msg.Tracking()<<"Amegic regularly finished."<<endl;
}



/*----------------------------------------------------------------------------------
  
  Process initialization

  ----------------------------------------------------------------------------------*/

bool Amegic::InitializeProcesses(int _runmode) 
{
  runmode = _runmode;
  msg.Debugging()<<"In Amegic::InitializeProcesses() "<<endl;
  procs = new All_Processes();
  procs->SetName("All_Processes");
  procs->SetAtoms(1);

  seldata      = new Selector_Data(path);

  beamon = isron = 0;
  arrows = 1;
  msg.Debugging()<<"Open file "<<path+string("/ISR.dat")<<endl;
  Data_Read dr(path+string("/ISR.dat"));
  beamtypes    = new int[2];
  beamtypes[0] = dr.GetValue<Beam_Type::code>("BEAM1");
  beamtypes[1] = dr.GetValue<Beam_Type::code>("BEAM2");
  if ((beamtypes[0] != Beam_Type::No) || (beamtypes[1] != Beam_Type::No)) { beamon = 1; ++arrows; }
  isrtypes     = new int[2];
  isrtypes[0]  = dr.GetValue<ISR_Type::code>("ISR1");
  isrtypes[1]  = dr.GetValue<ISR_Type::code>("ISR2");
  if ((isrtypes[0] != ISR_Type::No) || (isrtypes[1] != ISR_Type::No))     { isron = 1;  ++arrows; }
  double s     = sqr(AORGTOOLS::rpa.gen.Ecms());
  splimits     = new double[2];
  splimits[0]  = s*sqr(dr.GetValue<double>("SMIN"));
  splimits[1]  = s*sqr(dr.GetValue<double>("SMAX"));

  bunches      = new Flavour[2];
  beams        = new Flavour[2];
  partons      = new Flavour[2];
  for (int i=0;i<2;i++) {
    bunches[i] = Flavour(kf::none);
    beams[i]   = Flavour(kf::none);
    partons[i] = Flavour(kf::none);
  }
  int Nmax     = ReadProcesses(path);

  top          = new Topology(Max(Nmax,4));

  vector<double>           results;
  vector<Single_Process *> links;
  Vec4D * moms  = 0;

  switch (procs->InitAllProcesses(top,moms,results,links)) { 
  case 1  : 
    if (runmode==XS_MODE) {
      bool okay = 1;
      for (int j=0; j<links.size(); j++) if (links[j]->NumberOfDiagrams()!=IS_XS_FLAG) okay = 0;
      return okay;
    }
    return 1;
  case 0  : 
    if (runmode==XS_MODE) {
      msg.Error()<<"Some differentials are not available as fast function ! Try \"-P\"."<<endl;
      return 0;
    }
    msg.Error()<<"Some libraries were missing ! Type \"make install\" and rerun."<<endl;
    return 0;
  default :
    return 0;
  }
  return 1;
}

int Amegic::ReadProcesses(string path)
{ 
  msg.Debugging()<<"Open file "<<path+string("/Processes.dat")<<endl;
  ifstream from((path+string("/Processes.dat")).c_str());
  if (!from) {
    msg.Error()<<(path+string("/Processes.dat")).c_str()<<" not found !"<<endl;
    abort();
  }

  int nmax  = 0;
  int count = 0;
  Flavour  * FS, * flavs;
  Pol_Info * plFS, * plavs;
  char buffer[100];
  for(;from;) {
    bool error = 0;
    int  nIS = 0, nFS = 0;
    from.getline(buffer,100);
    if (buffer[0] != '\%' && strlen(buffer)>0) {
      string buf(buffer);
      msg.Debugging()<<"Analyse : "<<buf<<endl;
      // First : Check number of arrows
      int  position;
      int  flag  = 0;
      if ( (!beamon) && (isron) ) flag=1;
      if ( (!beamon) && (!isron)) flag=2;
      while (flag > -1) {
	position = buf.find(string("->"));
	if (position > 0) {
	  string ini = buf.substr(0,position);
	  buf        = buf.substr(position+2);
	  switch (flag) {
	  case 0: 
	    if (beamon) {
	      msg.Out()<<"Try to extract flavours for beams."<<endl;
	      if (ExtractFlavours(bunches,polbunches,ini) != 2) {
		msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
			   <<"Wrong number of bunches in "<<string(buffer)<<endl;
	      } 
	      if (!isron) flag++;
	      msg.Out()<<"found : "<<bunches[0]<<";"<<bunches[1]<<endl;
	    }
	    break;
	  case 1: 
	    if (isron) {
	      if (ExtractFlavours(beams,polbeams,ini) != 2) {
		msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
			   <<"Wrong number of beams in "<<string(buffer)<<endl;
	      } 
	    }
	    break;
	  case 2:
	    nIS = ExtractFlavours(partons,plpartons,ini);
	    if ((nIS< 1) || (nIS > 2)) {
	      msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
			 <<"Wrong number of partons in "<<string(buffer)<<endl;
	    } 
	    break;
	  default:
	    msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
		       <<"Funny number of arrows in "<<string(buffer)<<endl
		       <<"Have "<<flag<<" arrows with ISR : "<<isron
		       <<", Beam : "<<beamon<<endl;
	    error = 1;
	  }
	  flag++;
	  if (flag>3) {
	    msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
		       <<"Funny number of arrows in "<<string(buffer)<<endl
		       <<"Have "<<flag<<" arrows with ISR : "<<isron
		       <<", Beam : "<<beamon<<endl;
	    error = 1;
	  }
	}
	else {
	  if (flag<arrows) {
	    msg.Error()<<"Error in Amegic::ReadProcesses("<<path<<")."<<endl
		       <<"Funny number of arrows in "<<string(buffer)<<endl
		       <<"Have "<<flag<<" arrows with ISR : "<<isron
		       <<", Beam : "<<beamon<<endl;
	    error = 1;
	  }
	  flag = -1;
	}
      }
      if (!error) {
	if (!beam) {
	  beam = new Beam_Handler(beamtypes,bunches,polbunches,beams,polbeams,splimits);
	  if (beam->On()>0) {
	    if (isr) {
	      if (!(beam->CheckConsistency(bunches,beams))) {
		msg.Error()<<"Error in initialising Beam_Handler with ISR."<<endl
			   <<" "<<bunches[0]<<" -> "<<beams[0]<<", "
			   <<" "<<bunches[1]<<" -> "<<beams[1]<<endl
			   <<"  Delete it and ignore the process."<<endl;
		delete beam;
		error = 1;
	      }
	    }
	    else {
	      if (!(beam->CheckConsistency(bunches,partons))) {
		msg.Error()<<"Error in initialising Beam_Handler without ISR."<<endl
			   <<" "<<bunches[0]<<" -> "<<partons[0]<<", "
			   <<" "<<bunches[1]<<" -> "<<partons[1]<<endl
			   <<"  Delete it and ignore the process."<<endl;
		delete beam;
		error = 1;
	      }
	    }
	  }
	}
	if (!isr)  {
	  isr  = new ISR_Handler(isrtypes,beams,partons,splimits);
	  if (isr->On()>0) {
	    if (!(isr->CheckConsistency(beams,partons))) {
	      msg.Error()<<"Error in initialising ISR_Handler."<<endl
			 <<" "<<beams[0]<<" -> "<<partons[0]<<", "
			 <<" "<<beams[1]<<" -> "<<partons[1]<<endl
			 <<"  Delete it and ignore the process."<<endl;
	      delete isr;
	      error = 1;
	    }
	  }
	}
	if (beam->On()>0) {
	  if (isr->On()>0) {
	    if (!(beam->CheckConsistency(bunches,beams))) {
	      error = 1;
	    }
	  }
	  else {
	    if (!(beam->CheckConsistency(bunches,partons))) {
	      error = 1;
	    }
	  }
	}
	if (isr->On()>0) {
	  if (!(isr->CheckConsistency(beams,partons))) {
	      msg.Error()<<"Error in initialising ISR_Handler."<<endl
			 <<" "<<beams[0]<<" -> "<<partons[0]<<", "
			 <<" "<<beams[1]<<" -> "<<partons[1]<<endl
			 <<"  Ignore the process."<<endl;
	    error = 1;
	  }
	}
      }
      if (!error) {
	nFS = ExtractFlavours(FS,plFS,buf);
	++count;
	msg.Tracking()<<"Init process : "<<partons[0]<<" "<<partons[1]<<" -> ";
	for (int j=0;j<nFS;j++) msg.Tracking()<<FS[j]<<" ";
	msg.Tracking()<<"  (Check : "<<buf<<" )"<<endl;
	
	if (nIS+nFS > nmax) nmax = nIS+nFS;
	flavs = new Flavour[nIS+nFS];
	plavs = new Pol_Info[nIS+nFS];
	for (int i=0;i<nIS;i++) { flavs[i]     = partons[i]; plavs[i]     = plpartons[i]; }
	for (int i=0;i<nFS;i++) { flavs[i+nIS] = FS[i];      plavs[i+nIS] = plFS[i]; }
	bool single = 1;
	for (int i=0;i<nIS+nFS;i++) {
	  msg.Debugging()<<"Flavour "<<i<<" : "<<flavs[i]<<" : "<<flavs[i].Size()<<endl;
	  if (flavs[i].Size()>1) { single = 0; break; }
	} 
	if (single) {
	  if (!(procs->CheckExternalFlavours(2,partons,nFS,FS))) {
	    msg.Tracking()<<"Mismatch of flavours. Cannot initialize this process."<<endl
			  <<"flavours are "<<buf<<endl;
	  }
	  else {
	    procs->Add(new Single_Process(nIS,nFS,flavs,isr,beam,seldata,2,
					  rpa.me.KFactorScheme(),
					  rpa.me.ScaleScheme(),plavs,runmode));
	  }
	}
	else procs->Add(new Process_Group(nIS,nFS,flavs,isr,beam,seldata,2,
					  rpa.me.KFactorScheme(),
					  rpa.me.ScaleScheme(),plavs, runmode));
	delete [] flavs;
      }
    }
  }
  from.close();
  msg.Tracking()<<count<<" process(es) !"<<endl;
  return nmax;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int Amegic::ExtractFlavours(Flavour*& fl,Pol_Info*& pl,string buf)
{
  // kill initial spaces
  int ii[20];
  char pc[20],pp[20];
  double pd[20],angle[20];
  int count = 0;
  buf += string(" "); 
  int int_tag;
  short int i;
  for (i=0;i<20;i++) { pp[i]=' '; pd[i]=0.; pc[i]=' '; angle[i]=0.;}

  while(buf.length()>1) {

    if (buf[0]==' ') buf = buf.substr(1);
    else {
      int next = buf.find(string(" "));
      string pn = buf.substr(0,next);
      int nxt = pn.find(string("("));
      string number;
      if (nxt==string::npos) number = pn.substr(0,pn.length());
      else {
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]==')'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	  int lh=pn.find(string("l"));

	  if(lh!=string::npos) {
	    pc[count]='l';
	    pp[count]='+';
	    string ha = pn.substr(lh+1,pn.length()-lh);
	    std::strstream astream;
	    astream<<ha;
	    astream>>angle[count];
	    msg.Debugging()<<"*****Extract_Flavours:angle:"<<pn<<";"<<ha<<";"<<angle[count]<<endl;
	  }
	  else {
	    pp[count]=pn[pn.length()-1];
	    pc[count]=pn[pn.length()-2];
	    pn.erase(pn.length()-2,2);
	  }
	  std::strstream pstream;
	  pstream<<pn;
	  pstream>>pd[count];	  
	}
      }
      std::strstream sstream;
      sstream<<number;
      sstream>>ii[count];

      count++;
      buf = buf.substr(next);
    }
  }
  fl = new Flavour[count];
  pl = new Pol_Info[count];
  for (i=0;i<count;i++) {
    fl[i] = Flavour(kf::code(int(abs(double(ii[i])))));
    if (ii[i]<0) fl[i] = fl[i].Bar();

#ifdef Explicit_Pols
    int t1=mt::p_m, t2=mt::p_p;
    if(pc[i]=='l') { t1=mt::p_l0; t2=mt::p_l1; }
    pl[i].p_type = pc[i];
    pl[i].angle  = angle[i];
 
    int dof      = 1;
    int type;
    switch(pp[i]){
      case '-' : type = t1;     break;
      case '+' : type = t2;     break;
      case '0' : type = mt::p_l;break;
      default  : type = t1;
    }
    if(fl[i].IsFermion()) dof=2;
    if(fl[i].IsVector() &&  AMATOOLS::IsZero(fl[i].Mass())) dof=2;
    if(fl[i].IsVector() && !AMATOOLS::IsZero(fl[i].Mass())) dof=3;
    if(fl[i].IsTensor()) dof=5;
 
    if (AMATOOLS::IsZero(pd[i]-1.)) pl[i].init(1);
                               else pl[i].init(dof);

    if(!fl[i].IsTensor()){
      int tf[3]={t1,t2,mt::p_l};
      if(pl[i].num==1) {
	pl[i].type[0]=type;
	pl[i].factor[0]=dof;}
      else{
	for(int j=0;j<pl[i].num;j++){
	  pl[i].type[j]=tf[j];
	  if(pl[i].type[j]==type)
	    pl[i].factor[j]=1.+pd[i]*(dof-1.);
	  else pl[i].factor[j]=1.-pd[i];
	}
      }
    }
    else {
      pl[i].type[0]=mt::p_t1;pl[i].factor[0]=1.;
      pl[i].type[1]=mt::p_t2;pl[i].factor[1]=1.;
      pl[i].type[2]=mt::p_t3;pl[i].factor[2]=1.;
      pl[i].type[3]=mt::p_t4;pl[i].factor[3]=1.;
      pl[i].type[4]=mt::p_t5;pl[i].factor[4]=1.;
    }
  
    msg.Debugging()<<"*****Extract_Flavours:Pol:  "
		   <<pl[i].num<<" "<<pl[i].type[0]<<" "<<pl[i].type[1]<<" "
		   <<pl[i].factor[0]<<" "<<pl[i].factor[1]<<endl
		   <<"*****Extract_Flavours:Poltype:"<<pc[i]<<pl[i].p_type<<endl;
    
#else
    pl[i] = Pol_Info(fl[i]); 
#endif
  }

  return count;
}


int Amegic::ExtractFlavours(Flavour*& fl,double*& plfac,string buf)
{
  // kill initial spaces
  int ii[20];
  char pc[20],pp[20];
  double pd[20],angle[20];
  int count = 0;
  buf += string(" "); 
  int int_tag;
  short int i;
  for (i=0;i<20;i++) { pp[i]=' '; pd[i]=0.; pc[i]=' '; angle[i]=0.;}

  while(buf.length()>1) {
    if (buf[0]==' ') buf = buf.substr(1);
    else {
      int next = buf.find(string(" "));
      string pn = buf.substr(0,next);
      int nxt = pn.find(string("("));
      string number;
      if (nxt==string::npos) number = pn.substr(0,pn.length());
      else {
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]==')'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	  int lh=pn.find(string("l"));

	  if(lh!=string::npos) {
	    pc[count]='l';
	    pp[count]='+';
	    string ha = pn.substr(lh+1,pn.length()-lh);
	    std::strstream astream;
	    astream<<ha;
	    astream>>angle[count];
	    msg.Debugging()<<"*****Extract_Flavours:angle:"<<pn<<";"<<ha<<";"<<angle[count]<<endl;
	  }
	  else {
	    pp[count]=pn[pn.length()-1];
	    pc[count]=pn[pn.length()-2];
	    pn.erase(pn.length()-2,2);
	  }
	  std::strstream pstream;
	  pstream<<pn;
	  pstream>>pd[count];
	  
	  cout<<"The pp's : "<<pp[count]<<endl;
	  cout<<"The pd's : "<<pd[count]<<endl; 
	 }
      }
      std::strstream sstream;
      sstream<<number;
      sstream>>ii[count];

      count++;
      buf = buf.substr(next);
    }
  }
  fl  = new Flavour[count];
  plfac = new double[count];
  for (i=0;i<count;i++) {
    fl[i] = Flavour(kf::code(int(abs(double(ii[i])))));
    if (ii[i]<0) fl[i] = fl[i].Bar();
    plfac[i] = pd[i]; 
    if (pp[i]=='-') plfac[i]=-plfac[i];
  }
  return count;
}


/*----------------------------------------------------------------------------------
  
  Process evaluation

  ----------------------------------------------------------------------------------*/


bool Amegic::CalculateTotalXSec() {
  return procs->CalculateTotalXSec();
}

bool Amegic::LookUpXSec(double ycut,bool calc,string obs) {
  procs->SetResDir(resdir);
  return procs->LookUpXSec(ycut,calc,obs);
}


bool Amegic::PrepareXSecTables() {
  procs->SetResDir(resdir);
  return procs->PrepareXSecTables();
}

void Amegic::FifoOutput(double wt) 
{
  int nin  = procs->Selected()->Nin();
  int nout = procs->Selected()->Nout();


  ostream &  fifo = *p_fifo;
  fifo<<" event"<<endl;
  fifo<<"   "<<nout<<" number of particles"<<endl;
  fifo<<setiosflags(std::ios::scientific);
  fifo<<setiosflags(std::ios::uppercase);
  fifo<<std::setprecision(9);
  fifo<<" "<<std::setw(16)<<wt<<"  event weight"<<endl;  // unweighted

  //  fifo<<std::setiosflags(std::ios::scientific);
    //<<std::setprecision(9)<<std::setw(16);

  for (int j = 0;j<nin; j++) {
    fifo<<" "<<std::setw(3)<<int(procs->Selected()->Flavs()[j])<<" ";
    for (int k=1;k<4;++k) {
      fifo<<std::setw(16)<<procs->Selected()->Momenta()[j][k]<<" ";
    }
    fifo<<std::setw(16)<<procs->Selected()->Momenta()[j][0]<<endl;
  }
  msg.Debugging()<<"                      -> "<<endl;
  for (int j = nin;j<nin+nout; j++) {
    fifo<<" "<<std::setw(3)<<int(procs->Selected()->Flavs()[j])<<" ";
    for (int k=1;k<4;++k)
      fifo<<std::setw(16)<<procs->Selected()->Momenta()[j][k]<<" ";
    fifo<<std::setw(16)<<procs->Selected()->Momenta()[j][0]<<endl;
  }

  // "conti" or "endss"
}


void Amegic::SingleEvents() {
  for (int i=1;i<=rpa.gen.NumberOfEvents();++i) {
    msg.Debugging()<<"------------------------------------------------------------"<<endl
		   <<"----------------"<<i<<" th Event --------------------------"<<endl
		   <<"------------------------------------------------------------"<<endl;
    if (procs->OneEvent()) {
      FifoOutput();
      if (i==rpa.gen.NumberOfEvents()) {
	(*p_fifo)<<"endss"<<endl;
      }
      else {
	(*p_fifo)<<"conti"<<endl;
      }

      msg.Debugging()<<"OneEvent for "<<procs->Name()<<" successful !"<<endl
		     <<"    Selected "<<procs->Selected()->Name()<<" as subprocess."<<endl
		     <<"    Found "<<procs->Selected()->NumberOfDiagrams()
		     <<" Feynman diagrams."<<endl;
      for (int j = 0;j<procs->Selected()->Nin(); j++) {
	msg.Debugging()<<procs->Selected()->Flavs()[j]<<" : "
		       <<procs->Selected()->Momenta()[j]<<endl;
      }
      msg.Debugging()<<"                      -> "<<endl;
      for (int j = 0;j<procs->Selected()->Nout(); j++) {
	msg.Debugging()<<procs->Selected()->Flavs()[j+procs->Selected()->Nin()]<<" : "
		       <<procs->Selected()->Momenta()[j+procs->Selected()->Nin()]<<endl;
      }
      msg.Debugging()<<endl;
    }
  }
}
