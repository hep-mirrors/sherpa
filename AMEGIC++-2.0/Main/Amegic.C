#include "Amegic.H"

#include "Interaction_Model_Handler.H"
#include "Message.H"
#include "Topology.H"
#include "Algebra_Interpreter.H"
#include "MyStrStream.H"
#include "Process_Info.H"
#include <iomanip>
#include "Single_Process_MHV2.H"

//#define _DECAYS_

using namespace AMEGIC;
using namespace ATOOLS;
using namespace std;

/*----------------------------------------------------------------------------------
  
  Constructors

  ----------------------------------------------------------------------------------*/

Amegic::Amegic(std::string _path,std::string _file,
	       MODEL::Model_Base * _model) :
  m_path(_path), m_file(_file), m_nmax(0),m_minqcdjet(99), m_maxqcdjet(0), 
  m_maxjet(0), 
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
  rpa.gen.SetVariable("PROCESSFILE",processfile);
  p_procs             = new All_Processes();
  p_procs->SetName("All_Processes");
  p_procs->SetAtoms(1);

  string selfile      = p_dataread->GetValue<string>("SELECTORFILE",string("Selector.dat"));
  rpa.gen.SetVariable("SELECTORFILE",selfile);
  p_seldata           = new Selector_Data(m_path+selfile);

  ReadInProcessfile(processfile);

  m_count            = p_procs->Size();

  m_nmax             = ATOOLS::Max(m_nmax,4);
  p_top              = new Topology(m_nmax);

  ATOOLS::Vec4D * moms  = 0;

  msg.Events()<<"Amegic::InitializeProcesses : \n"
	   <<"   Process initialization started; new libraries may be created."<<std::endl;
  switch (p_procs->InitAllProcesses(p_model,p_top,moms)) { 
  case 1  : 
    msg.Events()<<"   No new libraries have been created."<<std::endl;
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
  Data_Read::SetTags(PHASIC::Integrable_Base::ScaleTags());
  PHASIC::scl::scheme _sc = (PHASIC::scl::scheme)(p_dataread->GetValue<int>("SCALE_SCHEME",0));
  ATOOLS::Data_Read::ResetTags();
  std::string _facscale   = p_dataread->GetValue<std::string>("FACTORIZATION_SCALE","");
  int    _kfactor_scheme  = p_dataread->GetValue<int>("KFACTOR_SCHEME",0);
  double _scale           = p_dataread->GetValue<double>("FIXED_SCALE",sqr(rpa.gen.Ecms()));
  int usepi               = p_dataread->GetValue<int>("PI",0);
  double scale_factor     = p_dataread->GetValue<double>("SCALE_FACTOR",1.);
  double factorization_scale_factor   = scale_factor*p_dataread->GetValue<double>("FACTORIZATION_SCALE_FACTOR",1.);
  double renormalization_scale_factor = scale_factor*p_dataread->GetValue<double>("RENORMALIZATION_SCALE_FACTOR",1.);
  rpa.gen.SetScaleFactors(factorization_scale_factor,renormalization_scale_factor);

  p_procs->SetUsePI(usepi);

  ifstream from((m_path+file).c_str());
  if (!from) {
    msg.Error()<<"ERROR in Amegic::InitializeProcesses : "<<endl
	       <<"   Process data file : "<<(m_path+file).c_str()<<" not found."<<std::endl
	       <<"   Abort program execution."<<endl;
    abort();
  }

  int         flag,position,njets;
  string      buf,ini,fin;
  int         nIS,   nFS,    nex;
  Flavour   * IS,  * FS,   * excluded, * flavs, *iflb, *fflb;
  Pol_Info  * plIS,* plFS, * pldummy,  * plavs, *iplb, *fplb;
  Process_Info* pinfo=0;
  int         order_ew,order_strong,kfactor_scheme;
  PHASIC::scl::scheme scale_scheme; 
  double      fixed_scale;
  vector<double>  venhance_factor,vmaxreduction_factor,vmaxredepsilon,vycut,vmaxerror;
  std::string enhance_function="1", factorization_scale;
  bool        print_graphs=false;
  int         enable_mhv=0; 
  string      selectorfile;
  while (from) {
    getline(from,buf);
    position = -1;
    if (buf.length()>0 && buf[0] != '%') {
      msg.LogFile()<<buf<<std::endl;
      vector<Process_Info*> AppPI;
      vector<int> AppKf;
      vector<int> AppNum;
      AppPI.clear();AppKf.clear();AppNum.clear();
      position   = buf.find(string("Process :")); 
      flag       = 0;
      if (position>-1 && position<(int)buf.length()) {
	pinfo    = new Process_Info(0,0);
	flag     = 1;
	buf      = buf.substr(position+9);
	position = buf.find(string("->"));
	if (position > 0) {
	  ini    = buf.substr(0,position);
	  fin    = buf.substr(position+2);
	  nIS    = ExtractFlavours(IS,plIS,ini);
	  nFS    = ExtractFlavours(FS,plFS,fin,&AppKf,&AppNum);
	  pinfo->AddSubList(nFS,FS,plFS);
	  if (AppKf.size()>AppPI.size()) AppPI.push_back(pinfo);
	  njets  = 0;
	  if ((nIS< 1) || (nIS > 2)) {
	    msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
		       <<"   Wrong number of partons in "<<buf<<endl;
	    flag = 0;
	  }
	  if (nIS==2) {
	    if (!(CF.ValidProcess(2,IS,nFS,FS))) {
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
	    if (pinfo) delete pinfo;
	  }
	  else {
	    order_ew = order_strong = -1;
	    selectorfile        = string("");
	    scale_scheme        = _sc;
	    factorization_scale = _facscale;
	    kfactor_scheme      = _kfactor_scheme;
	    fixed_scale         = _scale;
	    order_ew            = 99;
	    order_strong        = 99;
	    nex                 = 0;

	    vycut.clear();                vycut.push_back(-1.);
	    venhance_factor.clear();      venhance_factor.push_back(1.);
	    vmaxreduction_factor.clear(); vmaxreduction_factor.push_back(1.);
	    vmaxredepsilon.clear();       vmaxredepsilon.push_back(0.);
	    vmaxerror.clear();            vmaxerror.push_back(-1.);

	    enhance_function    = "1";
	    do {
	      getline(from,buf);
	      position = -1;
	      if (buf.length()>0 && buf[0] != '%') {
		msg.LogFile()<<buf<<std::endl;
		position = buf.find(string("Decay :"));
		if (position > -1) {
		  buf    = buf.substr(position+7);
		  position = buf.find(string("->"));
		  if (position > 0) {
		    ini    = buf.substr(0,position);
		    fin    = buf.substr(position+2);
		    int c=ExtractFlavours(iflb,iplb,ini);
		    if (iflb[0].Size()>1 || c!=1) {
		      msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
				 <<"   Invalid subsequent decay channel: "<<buf<<endl;
		      abort();
		    }

		    Process_Info *phelp = pinfo->FindDM(iplb[0].type[0]);
		    if (phelp) {
		      if (iflb[0]!=*(phelp->Flav())) {
			msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
				   <<"   Mismatch of decay flavour and identifier: "<<iflb[0]<<" "<<*(phelp->Flav())<<endl;
			abort();
		      }

		      int c=ExtractFlavours(fflb,fplb,fin,&AppKf,&AppNum);
		      if (c<2) {
			msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
				   <<"   Wrong number of partons in decay process: "<<buf<<endl;
			abort();
		      }
		      if (!(CF.ValidProcess(1,iflb,c,fflb))) {
			msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
				   <<"   Mismatch of flavours in decay process. Cannot initialize this process."<<endl
				   <<"   flavours are "<<buf<<endl;
			abort();
		      }
		      phelp->AddSubList(c,fflb,fplb);
		      if (AppKf.size()>AppPI.size()) AppPI.push_back(phelp);
		      delete [] fflb;
		      delete [] fplb;		      
		    }
		    delete [] iflb;
		    delete [] iplb;
		  }
		}

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
		  int helpsc;
		  str<<buf;
		  str>>helpsc;
		  scale_scheme = (PHASIC::scl::scheme)(helpsc);
		}

		position       = buf.find(string("Factorization scale :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(position+21);
		  Shorten(buf);
		  factorization_scale = buf;
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
		  ExtractMPvalues(buf,venhance_factor);
		}
		position       = buf.find(string("Enhance_Function :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  enhance_function=buf;
		}
		position       = buf.find(string("Max_Reduction :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  ExtractMPvalues(buf,vmaxreduction_factor);
		}
		position       = buf.find(string("Max_Epsilon :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  ExtractMPvalues(buf,vmaxredepsilon);
		}

		position       = buf.find(string("Integration_Error :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  ExtractMPvalues(buf,vmaxerror);
		}

		position       = buf.find(string("YCUT :"));
		if (position > -1) {
		  MyStrStream str;      
		  buf          = buf.substr(buf.find(":",position)+1);
		  Shorten(buf);
		  ExtractMPvalues(buf,vycut);
		}

		position       = buf.find(string("Print_Graphs"));
		if (position > -1) {
		  print_graphs=true;
		}
		position       = buf.find(string("Enable_MHV"));
		if (position > -1) {
		  enable_mhv=1;
		}
		position       = buf.find(string("Enable_MHV_ONLY"));
		if (position > -1) {
		  enable_mhv=4;
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
		    m_maxqcdjet = nmax;
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
	    } while (from && position==-1);
	    if (!pinfo || !pinfo->CheckCompleteness()) {
	      msg.Error()<<"Error in Amegic::InitializeProcesses("<<m_path+file<<")."<<endl
			 <<"   Missing decay processes! "<<endl;
	      if (pinfo) pinfo->Print();
	      abort();
	    }

	    pinfo->SetQCDjetNums();
	    vector<int> ii;
	    ii.clear();
	    for (size_t i=0;i<AppNum.size();i++) ii.push_back(0);
	    for (size_t i=0;i<AppNum.size();i++) {
	      AppPI[i]->m_maxqcdjets+=AppNum[i]-1;
	    }

	    do {
	      for (size_t i=0;i<AppNum.size();i++) {
		Flavour afl=Flavour(kf::code(iabs(AppKf[i])));
		Pol_Info apl=Pol_Info(afl);
		for (int j=0;j<ii[i];j++) AppPI[i]->m_sublist[0].push_back(new Process_Info(&afl,&apl));
	      }
	      Process_Info* pcinfo=new Process_Info(pinfo);

	      for (size_t i=0;i<AppNum.size();i++) {
		for (int j=0;j<ii[i];j++) {
		  delete AppPI[i]->m_sublist[0][AppPI[i]->m_sublist[0].size()-1];
		  AppPI[i]->m_sublist[0].pop_back();
		}
	      }
	      
	      nFS = pcinfo->TotalNout();
	      double  enhance_factor      = venhance_factor[0];
	      double  maxreduction_factor = vmaxreduction_factor[0];
	      double  maxredepsilon       = vmaxredepsilon[0];
	      double  ycut                = vycut[0];
	      double  maxerror            = vmaxerror[0];
	      if ((int)venhance_factor.size()>nFS&&venhance_factor[nFS]!=-1.) enhance_factor=venhance_factor[nFS];
	      if ((int)vmaxreduction_factor.size()>nFS&&vmaxreduction_factor[nFS]!=-1.) maxreduction_factor=vmaxreduction_factor[nFS];
	      if ((int)vmaxredepsilon.size()>nFS&&vmaxredepsilon[nFS]!=-1.) maxredepsilon=vmaxredepsilon[nFS];
	      if ((int)vycut.size()>nFS&&vycut[nFS]!=-1.) ycut=vycut[nFS];
	      if ((int)vmaxerror.size()>nFS&&vmaxerror[nFS]!=-1.) maxerror=vmaxerror[nFS];

	      if (nFS>m_maxjet) m_maxjet = nFS;
	      m_nmax = Max(m_nmax,pcinfo->Nmax(nIS));
	      flavs              = new Flavour[nIS+nFS];
	      plavs              = new Pol_Info[nIS+nFS];
	      for (int i=0;i<nIS;i++) { flavs[i]     = IS[i]; plavs[i]     = plIS[i]; }
	      pcinfo->GetFlavList(&flavs[nIS]);
	      pcinfo->GetPolList(&plavs[nIS]);
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
	    
	      int qcdjets(0);
	      double summass = 0.;
	      for (int i=0;i<nFS;i++) {
		summass += flavs[i+nIS].Mass();
	      }

	      Flavour *stcf = new Flavour[nFS];
	      int nst=pcinfo->GetStableFlavList(stcf);
	      for (int i=0;i<nst;i++) {
		if (stcf[i].Strong()) ++qcdjets;
	      }
	      delete[] stcf;
	      m_minqcdjet=Min(qcdjets,m_minqcdjet);
	      m_maxqcdjet=ATOOLS::Max(m_maxqcdjet,qcdjets);
	      if (summass<rpa.gen.Ecms()) {
		Process_Base * proc=NULL;
		if (single) {
		  if ((enable_mhv==1||enable_mhv==4) && CF.MHVCalculable(nIS,flavs,nFS,flavs+nIS))
		    proc = new Single_Process_MHV2(pcinfo,nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
						   order_strong,order_ew,
						   kfactor_scheme,scale_scheme,
						   fixed_scale,
						   plavs,nex,excluded,usepi,ycut,maxerror,enhance_function);
		  else if (enable_mhv!=4)
		    proc = new Single_Process(pcinfo,nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
					      order_strong,order_ew,
					      kfactor_scheme,scale_scheme,fixed_scale,
					      plavs,nex,excluded,usepi,ycut,maxerror,enhance_function);
		}
		else proc = new Process_Group(pcinfo,nIS,nFS,flavs,p_isr,p_beam,p_seldata,2,
					      order_strong,order_ew,
					      kfactor_scheme,scale_scheme,fixed_scale,
					      plavs,nex,excluded,usepi,ycut,maxerror,enhance_function,enable_mhv);
		if (proc) {
		  proc->SetEnhance(enhance_factor,maxreduction_factor,maxredepsilon);
		  proc->SetFactorizationScale(factorization_scale);
		  if (print_graphs) proc->SetPrintGraphs();
		  p_procs->Add(proc);
		  print_graphs=false;
		}
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
	      size_t k=0;
	      while (k<ii.size()) {
		ii[k]++;
		if (ii[k]==AppNum[k]) {
		  k++;
		  if (k<AppNum.size()) ii[k-1]=0;
		}
		else break;
	      }
	    } while (ii.size()!=0 && ii[ii.size()-1]<AppNum[AppNum.size()-1]);
	    if (pinfo) delete pinfo;
	  }
	}
      }
    }
  }
  p_procs->SetMaxJetNumber(m_maxjet);
}

void Amegic::ExtractMPvalues(string& str,vector<double>& dv)
{
  Shorten(str);
  int position;
  Algebra_Interpreter inter;
  inter.AddTag("E_CMS",ToString(rpa.gen.Ecms()));
  position = str.find(string("{"));
  if (position==-1) {
    dv[0]=ToType<double>(inter.Interprete(str));
    return;
  }
  string hstr = str.substr(0,position);
  Shorten(hstr);  
  double value = ToType<double>(inter.Interprete(hstr));
  str = str.substr(position+1,str.length()-position-2);
  do {
    position = str.find(string(","));
    int pos = -1;
    if (position>-1) {
      hstr = str.substr(0,position);
      Shorten(hstr);
      str = str.substr(position+1);
    }
    else hstr=str;
    MyStrStream strm;      
    strm<<hstr;
    strm>>pos;
    if (pos>0) {
      if (pos>=(int)dv.size()) dv.resize(pos+1,-1.);
      dv[pos]=value;
    }
  } while (position>-1);
}

int Amegic::ExtractFlavours(Flavour*& fl,Pol_Info*& pl,string buf,std::vector<int>* kf,vector<int>* num)
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
      nxt = number.find(string("["));
      if (nxt!=(int)string::npos) {
	pn = number;
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]==']'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	  pc[count]='d';
	  pp[count]=pn[0];
	}
      }
 
     nxt = number.find(string("{"));
      if (nxt!=(int)string::npos) {
	pn = number;
	number = pn.substr(0,nxt);
	pn.erase(0,nxt);
	if(pn[pn.length()-1]=='}'){
	  pn.erase(0,1);
	  pn.erase(pn.length()-1,1);
	}
	MyStrStream sstream;
	int nnum,nkf;
 	sstream<<pn;
 	sstream>>nnum;
 	MyStrStream kstream;
 	kstream<<number;
 	kstream>>nkf;	
	num->push_back(nnum+1);
	kf->push_back(nkf);
      }
      else {
	MyStrStream sstream;
	sstream<<number;
	sstream>>ii[count];

	count++;
      }
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
    if (pc[i]=='l') { t1=mt::p_l0; t2=mt::p_l1; }
    if (pc[i]!=' ') pl[i].pol_type = pc[i];
    pl[i].angle    = angle[i];
 
    int type;
    switch(pp[i]){
      case '-' : type = t1;     break;
      case '+' : type = t2;     break;
      case '0' : type = mt::p_l;break;
      default  : type = t1;
    }

    if (pc[i]=='d') pl[i].type[0] = pp[i];
    else {
      if(!fl[i].IsTensor()){
	int tf[3] = {t1,t2,mt::p_l};
	if (ATOOLS::IsZero(pd[i]-1.)) {
	  pl[i].type[0]=type;
	  pl[i].factor[0]=pl[i].num;
	  pl[i].num=1;
	}
	else{
	  for (int j=0;j<pl[i].num;j++) {
	    pl[i].type[j]=tf[j];
	    if(pl[i].type[j]==type)  pl[i].factor[j] = 1.+pd[i]*(pl[i].num-1.);
	    else  pl[i].factor[j] = 1.-pd[i];
	  }
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
  while (str.length()>0) {    
    if (int(str[0])==32 || int(str[0])==9) str = str.substr(1);
    else break;
  }
  //kill final spaces
  while (str.length()>0) {    
    if (int(str[str.length()-1])==32 ||
	//Tabulator
	int(str[str.length()-1])==9) str = str.substr(0,str.length()-1);
    else break;
  }
}


/*----------------------------------------------------------------------------------
  
  Process evaluation

  ----------------------------------------------------------------------------------*/


bool Amegic::CalculateTotalXSec(string _resdir,int mode) 
{
  bool success = p_procs->CalculateTotalXSec(_resdir);
  if (success) {
    p_procs->SetupEnhance();
    p_procs->SetWEventMode(mode);
    if (mode<0) {
      double nevents=rpa.gen.NumberOfEvents()/p_procs->TotalXS();
      p_procs->SetEvents(nevents);
    }
  }
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

ATOOLS::Blob_Data_Base *Amegic::SameEvent()
{
  return p_procs->SameEvent();
}


ATOOLS::Blob_Data_Base *Amegic::UnweightedEvent()
{
  return p_procs->OneEvent();
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


void Amegic::AddEvent(const double xs,const double validxs,const int ncounts)
{
  p_procs->AddEvent(xs,validxs,ncounts);
}

