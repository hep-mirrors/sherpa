#include "Amegic.H"
#include "Decay.H"

#include "Interaction_Model_Handler.H"
#include "Message.H"
#include "Topology.H"
#include "MyStrStream.H"

#include <iomanip>

//#define _DECAYS_

using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Amegic::Amegic(std::string _path,std::string _file,MODEL::Model_Base * _model) :
  m_path(_path), m_file(_file), m_maxjet(0), m_nmax(0),
  p_procs(NULL), p_decs(NULL), p_model(NULL), p_top(NULL), p_fifo(NULL),
  p_dataread(NULL), p_seldata(NULL), p_beam(NULL), p_isr(NULL)
{
  p_dataread         = new Data_Read(m_path+m_file);
  InitializeInteractionModel(_model);

  rpa.SetPath(m_path);
  
  bool use_fifo=0;
  if (rpa.gen.NumberOfEvents()>0 && use_fifo)
    p_fifo = new ofstream("flap.dat");
  else 
    p_fifo =0;
}
 

Amegic::~Amegic() {
  if (p_dataread) { delete p_dataread; p_dataread = NULL; }
  if (p_fifo)     { delete p_fifo;     p_fifo     = 0;    }
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
  msg.Tracking()<<"All together : "<<m_count<<" processes with up to "<<m_nmax<<" legs."<<endl;

  m_nmax             = AMATOOLS::Max(m_nmax,4);
  p_top              = new Topology(m_nmax);

  vector<double>           results;
  vector<Single_Process *> links;
  AMATOOLS::Vec4D * moms  = 0;

  switch (p_procs->InitAllProcesses(p_model,p_top,moms,results,links)) { 
  case 1  : 
    return 1;
  case 0  : 
    msg.Error()<<"Some libraries were missing ! Type make install and rerun."<<endl;
    return 0;
  default :
    return 0;
  }
  return 1;
}


bool Amegic::InitializeDecays() {
  int maxnumber = p_dataread->GetValue<int>("DECAY_PRODUCTS",3);
  if (p_top) {
    if (m_nmax<maxnumber) {
      msg.Error()<<"Error in Amegic::InitializeDecays()."<<endl
		 <<"   Potential inconsistency in number of legs of processes and decays."<<endl
		 <<"   Processes : "<<m_nmax<<", Decays : "<<1+maxnumber<<endl
		 <<"   Reduce number of legs for decay accordingly."<<endl; 
    } 
  }
  else p_top    = new Topology(1+maxnumber);
  p_decs        = new All_Decays(p_model,p_top);
  return p_decs->InitializeDecays();
}


void Amegic::InitializeInteractionModel(MODEL::Model_Base * _model)
{
  string modeltype   = p_dataread->GetValue<string>("SIGNAL_MODEL",string("SM"));
  string cplscheme   = p_dataread->GetValue<string>("COUPLING SCHEME",string("Running"));
  string massscheme  = p_dataread->GetValue<string>("YUKAWA MASSES",string("Running"));
  string widthscheme = p_dataread->GetValue<string>("WIDTH SCHEME",string("Complex"));

  Interaction_Model_Handler mh(_model);
  p_model = mh.GetModel(modeltype,cplscheme,massscheme);
  p_model->Init_Vertex();
}

void Amegic::ReadInProcessfile(string file) 
{
  int    _scale_scheme   = p_dataread->GetValue<int>("SCALE SCHEME",0);
  int    _kfactor_scheme = p_dataread->GetValue<int>("KFACTOR SCHEME",0);
  double _scale_factor   = p_dataread->GetValue<double>("SCALE FACTOR",1.);


  ifstream from((m_path+file).c_str());
  if (!from) {
    msg.Error()<<"Error in Amegic::InitializeProcesses : "<<endl
	       <<"   File : "<<(m_path+file).c_str()<<" not found ! Abort program execution."<<endl;
    abort();
  }

  char buffer[100];

  int         flag,position,njets;
  string      buf,ini,fin;
  int         nIS,   nFS,    nex;
  Flavour   * IS,  * FS,   * excluded, * flavs;
  Pol_Info  * plIS,* plFS, * pldummy,  * plavs;
  int         order_ew,order_strong,scale_scheme,kfactor_scheme; 
  double      scale_factor;
  string      selectorfile;
  MyStrStream str;      
  for(;from;) {
    from.getline(buffer,100);
    if (buffer[0] != '\%' && strlen(buffer)>0) {
      buf        = string(buffer);
      position   = buf.find(string("Process :")); 
      flag       = 0;
      if (position>-1 && position<buf.length()) {
	flag     = 1;
	buf      = buf.substr(position+9);
	position = buf.find(string("->"));
	if (position > 0) {
	  ini    = buf.substr(0,position);
	  fin    = buf.substr(position+2);
	  nIS    = ExtractFlavours(IS,plIS,ini);
	  nFS    = ExtractFlavours(FS,plFS,fin);
	  njets  = 0;
	  /*
	  for (int i=0;i<nFS;i++) { if (FS[i].Strong()) njets++; } 
	  if (njets>m_maxjet) {
	    for (int i=0;i<nFS;i++) msg.Out()<<FS[i]<<" ";
	    msg.Out()<<" -> "<<njets<<endl;
	    m_maxjet = njets;
	  }
	  */
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
	    order_ew       = 99;
	    order_strong   = 99;
	    nex            = 0;
	    do {
	      from.getline(buffer,100);
	      if (buffer[0] != '\%' && strlen(buffer)>0) {
		buf      = string(buffer);
		position = buf.find(string("Excluded particles :"));
		if (position > -1) {
		  buf    = buf.substr(position+20);
		  nex    = ExtractFlavours(excluded,pldummy,buf);
		}

		position = buf.find(string("Order electroweak :"));
		if (position > -1) {
		  buf    = buf.substr(position+19);
		  Shorten(buf);
		  str<<buf;
		  str>>order_ew;
		}

		position     = buf.find(string("Order strong :"));
		if (position > -1) {
		  buf    = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>order_strong;
		}

		position       = buf.find(string("Selector file :"));
		if (position > -1) {
		  buf          = buf.substr(position+15);
		  Shorten(buf);
		  selectorfile = buf;
		}

		position       = buf.find(string("Scale scheme :"));
		if (position > -1) {
		  buf          = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>scale_scheme;
		}

		position       = buf.find(string("Scale factor :"));
		if (position > -1) {
		  buf          = buf.substr(position+14);
		  Shorten(buf);
		  str<<buf;
		  str>>scale_factor;
		}

		position       = buf.find(string("KFactor scheme :"));
		if (position > -1) {
		  buf          = buf.substr(position+17);
		  Shorten(buf);
		  str<<buf;
		  str>>kfactor_scheme;
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
	    msg.Tracking()<<"Read in process :";
	    for (short int i=0;i<nIS;i++) msg.Tracking()<<" "<<IS[i].Name();
	    msg.Tracking()<<" -> ";
	    for (short int i=0;i<nFS;i++) msg.Tracking()<<FS[i].Name()<<" ";
	    msg.Tracking()<<" EW("<<order_ew<<"), QCD("<<order_strong<<")"<<endl;
	    if (nex>0) {
	      msg.Tracking()<<" Excluded particles : ";
	      for (short int i=0;i<nex;i++) msg.Tracking()<<excluded[i].Name()<<" ";
	      msg.Tracking()<<endl;
	    }

	    if (nIS+nFS>m_nmax) m_nmax = nIS+nFS;
	    flavs              = new Flavour[nIS+nFS];
	    plavs              = new Pol_Info[nIS+nFS];
	    for (int i=0;i<nIS;i++) { flavs[i]     = IS[i]; plavs[i]     = plIS[i]; } 
	    for (int i=0;i<nFS;i++) { flavs[i+nIS] = FS[i]; plavs[i+nIS] = plFS[i]; }
	    bool single        = 1;
	    for (int i=0;i<nIS+nFS;i++) {
	      if (flavs[i].Size()>1) { single = 0; break; }
	    } 

	    if (single) p_procs->Add(new Single_Process(nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
							order_strong,order_ew,
							kfactor_scheme,scale_scheme,scale_factor,
							plavs,nex,excluded));
 	    else p_procs->Add(new Process_Group(nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
						order_strong,order_ew,
						kfactor_scheme,scale_scheme,scale_factor,
						plavs,nex,excluded));
	    delete [] flavs;
	    delete [] plavs;
	  }
	}
      }
    }
  }
}


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

  msg.Debugging()<<"Count in ExtractFlavours : "<<count<<" : ";
  for (i=0;i<count;i++) msg.Debugging()<<ii[i]<<" ";
  msg.Debugging()<<endl;
  
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
 
    if (AMATOOLS::IsZero(pd[i]-1.)) pl[i].Init(1);
                               else pl[i].Init(dof);

    if(!fl[i].IsTensor()){
      int tf[3] = {t1,t2,mt::p_l};
      if(pl[i].num==1) {
	pl[i].type[0]=type;
	pl[i].factor[0]=dof;
      }
      else{
	for (int j=0;j<pl[i].num;j++){
	  pl[i].type[j]=tf[j];
	  if(pl[i].type[j]==type)  pl[i].factor[j] = 1.+pd[i]*(dof-1.);
	                     else  pl[i].factor[j] = 1.-pd[i];
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
  return p_procs->CalculateTotalXSec(_resdir);
}

bool Amegic::CalculateBranchingWidths(string _resdir) {
  return p_decs->CalculateBranchingWidths(_resdir);
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
    int nin  = p_procs->Selected()->Nin();
    int nout = p_procs->Selected()->Nout();


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
      fifo<<" "<<std::setw(3)<<int(p_procs->Selected()->Flavs()[j])<<" ";
      for (int k=1;k<4;++k) {
	fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][k]<<" ";
      }
      fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][0]<<endl;
    }
    msg.Debugging()<<"                      -> "<<endl;
    for (int j = nin;j<nin+nout; j++) {
      fifo<<" "<<std::setw(3)<<int(p_procs->Selected()->Flavs()[j])<<" ";
      for (int k=1;k<4;++k)
	fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][k]<<" ";
      fifo<<std::setw(16)<<p_procs->Selected()->Momenta()[j][0]<<endl;
    }
  }
  // "conti" or "endss"
}


bool Amegic::UnweightedEvent()
{
  if (p_procs->OneEvent()) {
    msg.Debugging()<<"OneEvent for "<<p_procs->Name()<<" successful !"<<endl
		   <<"    Selected "<<p_procs->Selected()->Name()<<" as subprocess."<<endl
		   <<"    Found "<<p_procs->Selected()->NumberOfDiagrams()
		   <<" Feynman diagrams."<<endl;
    for (int j = 0;j<p_procs->Selected()->Nin(); j++) {
      msg.Debugging()<<p_procs->Selected()->Flavs()[j]<<" : "
		     <<p_procs->Selected()->Momenta()[j]<<endl;
    }
    msg.Debugging()<<"                      -> "<<endl;
    for (int j = 0;j<p_procs->Selected()->Nout(); j++) {
      msg.Debugging()<<p_procs->Selected()->Flavs()[j+p_procs->Selected()->Nin()]<<" : "
		     <<p_procs->Selected()->Momenta()[j+p_procs->Selected()->Nin()]<<endl;
    }
    msg.Debugging()<<endl;
    return 1;
  }
  return 0;
}


void Amegic::SingleEvents() {
  for (int i=1;i<=rpa.gen.NumberOfEvents();++i) {
    msg.Debugging()<<"------------------------------------------------------------"<<endl
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
      msg.Debugging()<<"OneEvent for "<<p_procs->Name()<<" successful !"<<endl
		     <<"    Selected "<<p_procs->Selected()->Name()<<" as subprocess."<<endl
		     <<"    Found "<<p_procs->Selected()->NumberOfDiagrams()
		     <<" Feynman diagrams."<<endl;
      for (int j = 0;j<p_procs->Selected()->Nin(); j++) {
	msg.Debugging()<<p_procs->Selected()->Flavs()[j]<<" : "
		       <<p_procs->Selected()->Momenta()[j]<<endl;
      }
      msg.Debugging()<<"                      -> "<<endl;
      for (int j = 0;j<p_procs->Selected()->Nout(); j++) {
	msg.Debugging()<<p_procs->Selected()->Flavs()[j+p_procs->Selected()->Nin()]<<" : "
		       <<p_procs->Selected()->Momenta()[j+p_procs->Selected()->Nin()]<<endl;
      }
      msg.Debugging()<<endl;
    }
  }
}
