#include "Amegic.H"

#include "Interaction_Model_Handler.H"
#include "Message.H"
#include "Topology.H"
#include "MyStrStream.H"

#include <iomanip>

//#define _DECAYS_

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Amegic::Amegic(std::string _path,std::string _file,
	       MODEL::Model_Base * _model) :
  m_path(_path), m_file(_file), m_nmax(0),m_maxjet(0), 
  p_procs(NULL), p_decs(NULL), p_model(NULL), p_top(NULL), p_fifo(NULL),
  p_dataread(NULL), p_seldata(NULL), p_beam(NULL), p_isr(NULL)
{
  p_dataread          = new Data_Read(m_path+m_file);
  InitializeInteractionModel(_model);

  rpa.SetPath(m_path);
  
  bool use_fifo=0;
  if (rpa.gen.NumberOfEvents()>0 && use_fifo)
    p_fifo = new ofstream("flap.dat");
  else 
    p_fifo =0;
}
 


double Amegic::OverflowStatistics(Process_Base * proc,int level)
{
  if (level==0) msg_Tracking()<<"Amegic::OverflowStatistics for : "<<std::endl;
  if (proc==NULL) {
    if (p_procs==NULL) return 0.;
    proc=p_procs;
  }
  double sum=0.;
  if ((*proc)[0]==proc) {
    sum=proc->Overflow()/proc->Max();
    for (int i=0;i<level;++i) msg_Tracking()<<"  ";
    msg_Tracking()<<" "<<proc->Name()<<" "<<sum<<std::endl;
  }
  else {
    for (size_t i=0; i<proc->Size();++i) {
      sum+=OverflowStatistics((*proc)[i],level+1);
    }
    for (int i=0;i<level;++i) msg_Tracking()<<"  ";
    msg_Tracking()<<" "<<proc->Name()<<" "<<sum<<std::endl;
  }
  if (level==0 && sum!=0. && !msg.LevelIsTracking()) {
    int save_msg = msg.Level();
    msg.SetLevel(4);
    msg_Tracking()<<"Amegic::OverflowStatistics for : "<<std::endl;
    sum=0.;
    if ((*proc)[0]==proc) {
      sum=proc->Overflow()/proc->Max();
      for (int i=0;i<level;++i) msg_Tracking()<<"  ";
      msg_Tracking()<<" "<<proc->Name()<<" "<<sum<<std::endl;
    }
    else {
      for (size_t i=0; i<proc->Size();++i) {
	sum+=OverflowStatistics((*proc)[i],level+1);
      }
      for (int i=0;i<level;++i) msg_Tracking()<<"  ";
      msg_Tracking()<<" "<<proc->Name()<<" "<<sum<<std::endl;
    }
    msg.SetLevel(save_msg);
  }
  return sum;
}

Amegic::~Amegic() {
  OverflowStatistics();

  if (p_dataread) { delete p_dataread; p_dataread = NULL; }
  if (p_fifo)     { delete p_fifo;     p_fifo     = 0;    }
  if (p_model)    { delete p_model;    p_model    = 0;    }
  if (p_procs)    { delete p_procs;    p_procs    = 0;    }
  if (p_decs)     { delete p_decs;     p_decs     = 0;    }
  if (p_top)      { delete p_top;      p_top      = 0;    }
  if (p_seldata)  { delete p_seldata;  p_seldata  = 0;    }
}


/*---------------------------------------------------------------------------------

  Process initialization 

  ----------------------------------------------------------------------------------*/


bool Amegic::InitializeProcesses(BEAM::Beam_Spectra_Handler * _beam,PDF::ISR_Handler * _isr) {
  p_beam              = _beam; 
  p_isr               = _isr;
  string processfile  = p_dataread->GetValue<string>("PROCESSFILE",string("Processes.dat"));
  p_procs             = new All_Processes();
  p_procs->SetName("All_Processes");
  p_procs->SetAtoms(1);

  string selfile      = p_dataread->GetValue<string>("SELECTORFILE",string("Selector.dat"));
  p_seldata           = new Selector_Data(m_path+selfile);


  ReadInProcessfile(processfile);

  m_count            = p_procs->Size();

  m_nmax             = ATOOLS::Max(m_nmax,4);
  p_top              = new Topology(m_nmax);

  ATOOLS::Vec4D * moms  = 0;

  msg.Out()<<"Amegic::InitializeProcesses : \n"
	   <<"   Process initialization started; new libraries may be created."<<std::endl;
  switch (p_procs->InitAllProcesses(p_model,p_top,moms)) { 
  case 1  : 
    msg.Out()<<"   No new libraries have been created."<<std::endl;
    return 1;
  case 0  : 
    msg.Error()<<"Amegic::InitializeProcesses : "<<std::endl
	       <<"   Some new libraries were created and have to be compiled and linked."<<std::endl
	       <<om::bold<<"   Type \"./makelibs\" in '"
	       <<rpa.gen.Variable("SHERPA_CPP_PATH")<<"' and rerun."<<om::reset<<endl;
    return 0;
  default :
    return 0;
  }
  return 1;
}


bool Amegic::InitializeDecays(bool constructall) {
  int maxnumber = p_dataread->GetValue<int>("DECAY_PRODUCTS",3);
  if (maxnumber>3) {
    msg.Error()<<"ERROR in Amegic::InitializeDecays()."<<endl
	       <<"   Number of potential decay products larger than three."<<endl
	       <<"   This has not been implemented yet. Reduce number of decay products"<<endl
	       <<"   to three and continue run."<<endl;
    maxnumber = 3;
  }
  if (p_top) {
    if (m_nmax<maxnumber) {
      msg.Error()<<"ERROR in Amegic::InitializeDecays()."<<endl
		 <<"   Potential inconsistency in number of legs of processes and decays."<<endl
		 <<"   Processes : "<<m_nmax<<", Decays : "<<1+maxnumber<<endl
		 <<"   Reduce number of legs for decay accordingly."<<endl; 
    } 
  }
  else p_top    = new Topology(1+maxnumber);
  p_decs        = new All_Decays(p_model,p_top);
  if (constructall) return p_decs->InitializeDecayTables();
  return 1;
}


void Amegic::InitializeInteractionModel(MODEL::Model_Base * _model)
{
  string modeltype   = p_dataread->GetValue<string>("SIGNAL_MODEL",string("SM"));
  string cplscheme   = p_dataread->GetValue<string>("COUPLING_SCHEME",string("Running"));
  string massscheme  = p_dataread->GetValue<string>("YUKAWA_MASSES",string("Running"));
  string widthscheme = p_dataread->GetValue<string>("WIDTH_SCHEME",string("Fixed"));

  Interaction_Model_Handler mh(_model);
  p_model = mh.GetModel(modeltype,cplscheme,massscheme);
  p_model->Init_Vertex();
}

void Amegic::ReadInProcessfile(string file) 
{
  int    _scale_scheme   = p_dataread->GetValue<int>("SCALE_SCHEME",0);
  int    _kfactor_scheme = p_dataread->GetValue<int>("KFACTOR_SCHEME",0);
  double _scale_factor   = p_dataread->GetValue<double>("SCALE_FACTOR",1.);
  double _scale          = p_dataread->GetValue<double>("FIXED_SCALE",sqr(rpa.gen.Ecms()));


  ifstream from((m_path+file).c_str());
  if (!from) {
    msg.Error()<<"ERROR in Amegic::InitializeProcesses : "<<endl
	       <<"   Process data file : "<<(m_path+file).c_str()<<" not found."<<std::endl
	       <<"   Abort program execution."<<endl;
    abort();
  }

  char buffer[100];

  int         flag,position,njets;
  string      buf,ini,fin;
  int         nIS,   nFS,    nex;
  Flavour   * IS,  * FS,   * excluded, * flavs;
  Pol_Info  * plIS,* plFS, * pldummy,  * plavs;
  int         order_ew,order_strong,scale_scheme,kfactor_scheme; 
  double      scale_factor,fixed_scale;
  double      enhance_factor=1.,maxreduction_factor=1.;
  bool        print_graphs=false;
  string      selectorfile;
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '%' && strlen(buffer)>0) {
      msg.LogFile()<<buffer<<std::endl;
      buf        = string(buffer);
      position   = buf.find(string("Process :")); 
      flag       = 0;
      if (position>-1 && position<(int)buf.length()) {
	flag     = 1;
	buf      = buf.substr(position+9);
	position = buf.find(string("->"));
	if (position > 0) {
	  ini    = buf.substr(0,position);
	  fin    = buf.substr(position+2);
	  nIS    = ExtractFlavours(IS,plIS,ini);
	  nFS    = ExtractFlavours(FS,plFS,fin);
	  njets  = 0;
	  if (nFS>m_maxjet) m_maxjet = nFS;
	  if ((nIS< 1) || (nIS > 2)) {
	    msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
		       <<"   Wrong number of partons in "<<string(buffer)<<endl;
	    flag = 0;
	  }
	  if (nIS==2) {
	    if (!(p_procs->CheckExternalFlavours(2,IS,nFS,FS))) {
	      msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
			 <<"   Mismatch of flavours. Cannot initialize this process."<<endl
			 <<"   flavours are "<<buf<<endl;
	      flag = 0;
	    }
	    if (!(p_isr->CheckConsistency(IS))) {
	      msg.Error()<<"Error in initialising ISR_Handler."<<endl
			 <<" "<<p_isr->Flav(0)<<" -> "<<IS[0]<<", "
			 <<" "<<p_isr->Flav(1)<<" -> "<<IS[1]<<endl
			 <<"  Delete it and ignore the process."<<endl;
	      flag = 0;
	    }
	  }
	  if (flag==0) {
	    delete [] IS;
	    delete [] plIS;
	    delete [] FS;
	    delete [] plFS;
	  }
	  else {
	    order_ew = order_strong = -1;
	    selectorfile   = string("");
	    scale_scheme   = _scale_scheme;
	    scale_factor   = _scale_factor;
	    kfactor_scheme = _kfactor_scheme;
	    fixed_scale    = _scale;
	    order_ew       = 99;
	    order_strong   = 99;
	    nex            = 0;
	    do {
	      from.getline(buffer,100);
	      if (buffer[0] != '%' && strlen(buffer)>0) {
		msg.LogFile()<<buffer<<std::endl;
		buf      = string(buffer);
		position = buf.find(string("Excluded particles :"));
		if (position > -1) {
		  buf    = buf.substr(position+20);
		  nex    = ExtractFlavours(excluded,pldummy,buf);
		}

		position = buf.find(string("Order electroweak :"));
		if (position > -1) {
		  MyStrStream str;
		  buf    = buf.substr(position+19);
		  Shorten(buf);
		  str<<buf;
		  str>>order_ew;
		}

		position     = buf.find(string("Order strong :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf    = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>order_strong;
		}

		position       = buf.find(string("Selector file :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+15);
		  Shorten(buf);
		  selectorfile = buf;
		}

		position       = buf.find(string("Scale scheme :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>scale_scheme;
		}

		position       = buf.find(string("Scale factor :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>scale_factor;
		}

		position       = buf.find(string("KFactor scheme :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+16);
		  Shorten(buf);
		  str<<buf;
		  str>>kfactor_scheme;
		}

		position       = buf.find(string("Fixed scale :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+13);
		  Shorten(buf);
		  str<<buf;
		  str>>fixed_scale;
		}

		position       = buf.find(string("Enhance_Factor :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  str<<buf;
		  str>>enhance_factor;
		}
		position       = buf.find(string("Max_Reduction :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  str<<buf;
		  str>>maxreduction_factor;
		}

		position       = buf.find(string("Print_Graphs"));
		if (position > -1) {
		  print_graphs=true;
		}

		position     = buf.find(string("N_Max :"));
		if (position > -1) {
		  int nmax;
		  MyStrStream str;      
		  buf    = buf.substr(position+7);
		  Shorten(buf);
		  str<<buf;
		  str>>nmax;
		  if (nmax>m_maxjet) {
		    msg.Out()<<" WARNING: setting max n to "<<nmax<<std::endl;
		    m_maxjet = nmax;
		  }		  
		}

		position       = buf.find(string("End process"));
		if (!from) {
		  msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
			     <<"   End of file reached without 'End process'-tag."<<endl
			     <<"   Continue and hope for the best."<<endl;
		  position     = 0;
		}
	      }
	    }
	    while (position==-1);

	    if (nIS+nFS>m_nmax) m_nmax = nIS+nFS;
	    flavs              = new Flavour[nIS+nFS];
	    plavs              = new Pol_Info[nIS+nFS];
	    for (int i=0;i<nIS;i++) { flavs[i]     = IS[i]; plavs[i]     = plIS[i]; } 
	    for (int i=0;i<nFS;i++) { flavs[i+nIS] = FS[i]; plavs[i+nIS] = plFS[i]; }
	    bool single        = 1;
	    for (int i=0;i<nIS+nFS;i++) {
	      if (flavs[i].Size()>1) { single = 0; break; }
	    } 

	    // for beam
	    int bsh_pol=p_beam->Polarisation();
	    int beam_is_poled[2]={bsh_pol&1,bsh_pol&2};
	    for (int i=0;i<nIS;++i) {
	      if (plavs[i].DoFNumber()>1 && beam_is_poled[i]) { single = 0; break; }
	    } 
	    

	    double summass = 0.;
	    for (int i=0;i<nFS;i++)
	      summass += flavs[i+nIS].Mass();

	    if (summass<rpa.gen.Ecms()) {
	      Process_Base * proc=NULL;
	      if (single) proc = new Single_Process(nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
						     order_strong,order_ew,
						     -kfactor_scheme,-scale_scheme,scale_factor,fixed_scale,
						     plavs,nex,excluded);
	      else proc = new Process_Group(nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
					    order_strong,order_ew,
					    -kfactor_scheme,-scale_scheme,scale_factor,fixed_scale,
					    plavs,nex,excluded);
	      proc->SetEnhance(enhance_factor,maxreduction_factor);
	      if (print_graphs) proc->SetPrintGraphs();
	      p_procs->Add(proc);
	      print_graphs=false;
	    }
	    else {
	      msg.Out()<<"Ignored process: ";
	      for (short int i=0;i<nIS;i++) msg.Out()<<" "<<IS[i].Name();
	      msg.Out()<<" -> ";
	      for (short int i=0;i<nFS;i++) msg.Out()<<FS[i].Name()<<" ";
	      msg.Out()<<", kinematically not allowed."<<endl;
	    }	    

	    delete [] flavs;
	    delete [] plavs;
	  }
	}
      }
    }
  }
  p_procs->SetMaxJetNumber(m_maxjet);
}


int Amegic::ExtractFlavours(Flavour*& fl,Pol_Info*& pl,string buf)
{
  // kill initial spaces
  int ii[20];
  char pc[20],pp[20];
  double pd[20],angle[20];
  int count = 0;
  buf += string(" "); 
  short int i;
  for (i=0;i<20;i++) { pp[i]=' '; pd[i]=0.; pc[i]=' '; angle[i]=0.;}

  while(buf.length()>1) {
    if (buf[0]==' ') buf = buf.substr(1);
    else {
      int next = buf.find(string(" "));
      string pn = buf.substr(0,next);
      int nxt = pn.find(string("("));
      string number;
      if (nxt==(int)string::npos) number = pn.substr(0,pn.length());
      else {
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]==')'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	  int lh=pn.find(string("l"));

	  if(lh!=(int)string::npos) {
	    pc[count]='l';
	    pp[count]='+';
	    string ha = pn.substr(lh+1,pn.length()-lh);
	    MyStrStream astream;
	    astream<<ha;
	    astream>>angle[count];
	  }
	  else {
	    pp[count]=pn[pn.length()-1];
	    pc[count]=pn[pn.length()-2];
	    pn.erase(pn.length()-2,2);
	  }
	  MyStrStream pstream;
	  pstream<<pn;
	  pstream>>pd[count];	  
	}
      }
      MyStrStream sstream;
      sstream<<number;
      sstream>>ii[count];

      count++;
      buf = buf.substr(next);
    }
  }

  fl = new Flavour[count];
  pl = new Pol_Info[count];
  
  for (i=0;i<count;i++) {
    fl[i] = Flavour(kf::code(iabs(ii[i])));
    if (ii[i]<0) fl[i] = fl[i].Bar();
    pl[i]=Pol_Info(fl[i]);

#ifdef Explicit_Pols
    int t1=mt::p_m, t2=mt::p_p;
    if(pc[i]=='l') { t1=mt::p_l0; t2=mt::p_l1; }
    if (pc[i]!=' ') pl[i].pol_type = pc[i];
    pl[i].angle    = angle[i];
 
    int type;
    switch(pp[i]){
      case '-' : type = t1;     break;
      case '+' : type = t2;     break;
      case '0' : type = mt::p_l;break;
      default  : type = t1;
    }
    if(!fl[i].IsTensor()){
      int tf[3] = {t1,t2,mt::p_l};
      if (ATOOLS::IsZero(pd[i]-1.)) {
	pl[i].type[0]=type;
	pl[i].factor[0]=pl[i].num;
	pl[i].num=1;
      }
      else{
	for (int j=0;j<pl[i].num;j++){
	  pl[i].type[j]=tf[j];
	  if(pl[i].type[j]==type)  pl[i].factor[j] = 1.+pd[i]*(pl[i].num-1.);
	                     else  pl[i].factor[j] = 1.-pd[i];
	}
      }
    }
  
#else
    pl[i] = Pol_Info(fl[i]); 
#endif
  }

  return count;
}

void Amegic::Shorten(std::string& str) {
  //kill initial spaces
  for (;;) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  //kill final spaces
  for (;;) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
}


/*----------------------------------------------------------------------------------
  
  Process evaluation

  ----------------------------------------------------------------------------------*/


bool Amegic::CalculateTotalXSec(string _resdir) {
  bool success = p_procs->CalculateTotalXSec(_resdir);
  if (success)   p_procs->SetupEnhance();
  return success;
}

bool Amegic::CalculateBranchingWidths(string _resdir) {
  return p_decs->CalculateWidths(_resdir);
}

void Amegic::SetResDir(string _respath) {
  m_respath = _respath;
  p_procs->SetResDir(_respath);
}

bool Amegic::LookUpXSec(double ycut,bool calc,string obs) {
  return p_procs->LookUpXSec(ycut,calc,obs);
}

bool Amegic::PrepareXSecTables() {
  return p_procs->PrepareXSecTables();
}

void Amegic::FifoOutput(double wt) 
{
  if (p_fifo) {
    int nin  = p_procs->Selected()->NIn();
    int nout = p_procs->Selected()->NOut();


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
      fifo<<" "<<std::setw(3)<<int(p_procs->Selected()->Flavours()[j])<<" ";
      for (int k=1;k<4;++k) {
	fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][k]<<" ";
      }
      fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][0]<<endl;
    }
    for (int j = nin;j<nin+nout; j++) {
      fifo<<" "<<std::setw(3)<<int(p_procs->Selected()->Flavours()[j])<<" ";
      for (int k=1;k<4;++k)
	fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][k]<<" ";
      fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][0]<<endl;
    }
  }
  // "conti" or "endss"
}

bool Amegic::SameEvent()
{
  return p_procs->SameEvent();
}


bool Amegic::UnweightedEvent()
{
  if (p_procs->OneEvent()) return 1;
  return 0;
}


void Amegic::SingleEvents() {
  for (int i=1;i<=rpa.gen.NumberOfEvents();++i) {
    msg_Events()<<"------------------------------------------------------------"<<endl
		<<"----------------"<<i<<" th Event --------------------------"<<endl
		<<"------------------------------------------------------------"<<endl;
    if (p_procs->OneEvent()) {
      if (p_fifo) {
	FifoOutput();
	if (i==rpa.gen.NumberOfEvents()) {
	  (*p_fifo)<<" endss"<<endl;
	}
	else {
	  (*p_fifo)<<" conti"<<endl;
	}
      }
    }
  }
}


ATOOLS::Blob_Data_Base * Amegic::SameWeightedEvent()
{
  return p_procs->SameWeightedEvent();
}

ATOOLS::Blob_Data_Base *  Amegic::WeightedEvent()
{
  return p_procs->WeightedEvent();
  //  if (weight>0.) return weight;
  //  return 0.;
}
