#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "PDF/Main/Shower_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#ifdef USING__GZIP
#include "ATOOLS/Org/Gzip_Stream.H"
#endif

#include <unistd.h>

using namespace SHERPA;
using namespace PHASIC;
using namespace PDF;
using namespace ATOOLS;

Matrix_Element_Handler::Matrix_Element_Handler
(const std::string &dir,const std::string &file,
 const std::string &processfile,const std::string &selectorfile):
  m_gens(dir, file),
  p_proc(NULL), p_beam(NULL), p_isr(NULL), p_model(NULL),
  m_path(dir), m_file(file), m_processfile(processfile),
  m_selectorfile(selectorfile), m_eventmode(0),
  p_shower(NULL), m_totalxs(0.0), 
  m_ranidx(0), p_ranin(NULL), p_ranout(NULL)
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  if (!read.ReadFromFile(m_respath,"RESULT_DIRECTORY")) m_respath="./Results";
  if (rpa.gen.Variable("PATH_PIECE")!="")
    m_respath=rpa.gen.Variable("PATH_PIECE")+"/"+m_respath;
  std::string evtm;
  if (!read.ReadFromFile(evtm,"EVENT_GENERATION_MODE")) evtm="Unweighted";
  if (evtm=="Unweighted" || evtm=="U") m_eventmode=1;
  else m_eventmode=0;
  if (!read.ReadFromFile(m_seedmode,"EVENT_SEED_MODE")) m_seedmode=0;
  else msg_Info()<<METHOD<<"(): Set seed mode "<<m_seedmode<<"."<<std::endl;
  std::string seedfile;
  if (!read.ReadFromFile(seedfile,"EVENT_SEED_FILE")) 
    seedfile="ran.stat."+ToString(rpa.gen.RandomSeed());
  else msg_Info()<<METHOD<<"(): Set seed file "<<seedfile<<"."<<std::endl;
#ifdef USING__GZIP
  seedfile+=".gz";
#endif
  if (m_seedmode==1) {
#ifdef USING__GZIP
    p_ranin = new ATOOLS::igzstream(seedfile.c_str());
#else
    p_ranin = new std::ifstream(seedfile.c_str());
#endif
    if (p_ranin!=NULL && !p_ranin->good()) THROW
      (fatal_error,"Cannot initialize random generator status file");
  }
  else if (m_seedmode==2) {
#ifdef USING__GZIP
    p_ranout = new ATOOLS::ogzstream(seedfile.c_str());
#else
    p_ranout = new std::ofstream(seedfile.c_str());
#endif
    if (p_ranout!=NULL && !p_ranout->good()) THROW
      (fatal_error,"Cannot initialize random generator status file");
  }
  if (!read.ReadFromFile(m_nlomode,"NLO_Mode")) m_nlomode=1;
}

Matrix_Element_Handler::~Matrix_Element_Handler()
{
  if (p_ranin) delete p_ranin;
  if (p_ranout) delete p_ranout;
}

bool Matrix_Element_Handler::CalculateTotalXSecs() 
{
  for (ME_Generators::const_iterator it=m_gens.begin();it!=m_gens.end();++it) {
    MakeDir(m_respath+"/"+(*it)->Name(),false);
  }

  m_totalxs=0.0;
  bool okay(true);
  for (size_t i=0;i<m_procs.size();++i) {
    m_procs[i]->SetUpThreading();
    m_procs[i]->SetLookUp(true);
    std::string name=m_procs[i]->Generator()?m_procs[i]->Generator()->Name():"";
    if (!m_procs[i]->CalculateTotalXSec(m_respath+"/"+name,false)) okay=false;
    m_procs[i]->SetLookUp(false);
    m_totalxs+=m_procs[i]->Integrator()->TotalXS();
    m_procs[i]->Integrator()->SetUpEnhance();
  }
  return okay;
}

void Matrix_Element_Handler::SetRandomSeed()
{
  if (m_seedmode==1) {
    m_ranidx=ran.ReadInStatus(*p_ranin,m_ranidx);
    if (m_ranidx==std::string::npos) {
      msg_Error()<<METHOD<<"(): Status file read error. Abort."<<std::endl;
      abort();
    }
  }
  else if (m_seedmode==2) {
    m_ranidx=ran.WriteOutStatus(*p_ranout,m_ranidx);
    if (m_ranidx==std::string::npos) {
      msg_Error()<<METHOD<<"(): Status file write error. Abort."<<std::endl;
      abort();
    }
  }
}

bool Matrix_Element_Handler::GenerateOneEvent() 
{
  SetRandomSeed();
  p_isr->SetPDFMember();
  if (m_eventmode==1) return GenerateUnweightedEvent(); 
  return GenerateWeightedEvent();
}

bool Matrix_Element_Handler::GenerateUnweightedEvent() 
{
  p_proc=NULL;
  double sum(0.0);
  for (size_t i(0);i<m_procs.size();++i)
    sum+=m_procs[i]->Integrator()->SelectionWeight();
  double disc(sum*ran.Get());
  sum=0.0;
  Process_Base *proc(NULL);
  for (size_t i(0);i<m_procs.size();++i) {
    if ((sum+=m_procs[i]->Integrator()->SelectionWeight())>=disc) {
      proc=m_procs[i];
      break;
    }
  }
  if (proc==NULL) THROW(fatal_error,"No process selected");
  PHASIC::Weight_Info *info=proc->OneEvent();
  if (info==NULL) THROW(fatal_error,"No weight information.");
  m_evtinfo=*info;
  delete info;
  p_proc=proc->Selected();
  m_evtinfo.m_ntrial=1.0/p_proc->Integrator()->EnhanceFactor();
  double wf(m_totalxs*rpa.Picobarn()/
	    p_proc->Integrator()->EnhanceFactor());
  if (p_proc->GetSubevtList()) 
    (*p_proc->GetSubevtList())*=wf/m_evtinfo.m_weight;
  m_evtinfo.m_weight=wf;
  return m_evtinfo.m_weight!=0.0;
}

bool Matrix_Element_Handler::GenerateWeightedEvent() 
{
  double addn(0.0), sum(0.0);
  for (size_t i(0);i<m_procs.size();++i)
    sum+=m_procs[i]->Integrator()->SelectionWeight();
  while (true) {
    double disc(sum*ran.Get()), csum(0.0);
    Process_Base *proc(NULL);
    for (size_t i(0);i<m_procs.size();++i) {
      if ((csum+=m_procs[i]->Integrator()->SelectionWeight())>=disc) {
	proc=m_procs[i];
	break;
      }
    }
    if (proc==NULL) THROW(fatal_error,"No process selected");
    PHASIC::Weight_Info *info=proc->WeightedEvent();
    p_proc=proc->Selected();
    double sw(p_proc->Integrator()->SelectionWeight()/sum);
    addn+=1.0;
    if (info) {
      m_evtinfo=*info;
      delete info;
      double wf(rpa.Picobarn()/sw);
      m_evtinfo.m_weight*=wf;
      if (p_proc->GetSubevtList()) {
	(*p_proc->GetSubevtList())*=wf;
	p_proc->GetSubevtList()->MultMEwgt(wf);
      }
      if (p_proc->GetMEwgtinfo()) (*p_proc->GetMEwgtinfo())*=wf;
      m_evtinfo.m_ntrial=addn;
      return true;
    }
  }
  return false;
}

std::vector<Process_Base*> Matrix_Element_Handler::InitializeProcess(const Process_Info &pi)
{
  DEBUG_FUNC(pi);
  std::vector<Process_Base*> procs;
  if (pi.m_fi.NLOType()==nlo_type::lo) {
    Process_Base *proc(m_gens.InitializeProcess(pi, true));
    if (proc) {
      m_procs.push_back(proc);
      procs.push_back(proc);
    }
    return procs;
  }
  else {
    if (m_nlomode==1) {
      if (pi.m_fi.m_nloqcdtype==nlo_type::lo) {
	if (pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop|nlo_type::born)) {
	  Process_Info rpi(pi);
	  rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop|
						 nlo_type::polecheck|nlo_type::born));
	  procs.push_back(m_gens.InitializeProcess(rpi,true));
	}
      }
      else if (pi.m_fi.m_nloewtype==nlo_type::lo) {
	if (pi.m_fi.NLOType()&nlo_type::born) {
	  Process_Info rpi(pi);
	  rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&nlo_type::born);
	  procs.push_back(m_gens.InitializeProcess(rpi,true));
	}
	if (pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop)) {
	  Process_Info rpi(pi);
	  rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&(nlo_type::vsub|nlo_type::loop|nlo_type::polecheck));
	  procs.push_back(m_gens.InitializeProcess(rpi,true));
	}
      }
      if (pi.m_fi.NLOType()&nlo_type::real || pi.m_fi.NLOType()&nlo_type::rsub){
        // if real or rsub is requested, the extra jet is not yet contained
        // in the process info, but has to be added here
        Process_Info rpi(pi);
	rpi.m_fi.SetNLOType(pi.m_fi.NLOType()&(nlo_type::real|nlo_type::rsub));
        if (pi.m_fi.m_nloqcdtype==nlo_type::lo) {
          rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
        }
        else if (pi.m_fi.m_nloewtype==nlo_type::lo) {
          rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
        }
        procs.push_back(m_gens.InitializeProcess(rpi,true));
      }
      for (size_t i(0);i<procs.size();i++) m_procs.push_back(procs[i]);
    }
  }
  return procs;
}

bool Matrix_Element_Handler::InitializeProcesses
(MODEL::Model_Base *model,
 BEAM::Beam_Spectra_Handler *beam,PDF::ISR_Handler *isr)
{
  /*
    This is the basis for all CKKW and process interplay.
    Don't even try to think about modifying 
    either this routine or any of its dependencies !!!
  */
  p_beam=beam; p_isr=isr; p_model=model;
  if (!m_gens.InitializeGenerators(model,beam,isr)) return false;
#ifdef USING__Threading
  double rbtime(ATOOLS::rpa.gen.Timer().RealTime());
#endif
  double btime(ATOOLS::rpa.gen.Timer().UserTime());
  BuildProcesses();
  if (msg_LevelIsTracking()) msg_Info()<<"Process initialization";
#ifdef USING__Threading
  double retime(ATOOLS::rpa.gen.Timer().RealTime());
#endif
  double etime(ATOOLS::rpa.gen.Timer().UserTime());
  Data_Reader psread;
  psread.SetAddCommandLine(false);
  psread.AddWordSeparator("\t");
  psread.SetInputPath("/proc/"+ToString(getpid())+"/");
  psread.SetInputFile("status");
  std::vector<std::string> musage;
  if (psread.VectorFromFile(musage,"VmSize:"))
    for (size_t i(1);i<musage.size();++i) musage.front()+=" "+musage[i];
#ifdef USING__Threading
  msg_Info()<<" done ( "<<(musage.size()>0?musage.front():"?")<<", "
	    <<FormatTime(size_t(retime-rbtime))<<" / "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
#else
  msg_Info()<<" done ( "<<(musage.size()>0?musage.front():"?")<<", "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
#endif
  if (m_procs.empty() && m_gens.size()>0)
    THROW(normal_exit,"No hard process found");
  msg_Info()<<METHOD<<"(): Performing tests "<<std::flush;
#ifdef USING__Threading
  rbtime=retime;
#endif
  btime=etime;
  bool res(m_gens.PerformTests());
#ifdef USING__Threading
  retime=ATOOLS::rpa.gen.Timer().RealTime();
#endif
  etime=ATOOLS::rpa.gen.Timer().UserTime();
  psread.RereadInFile();
  if (psread.VectorFromFile(musage,"VmSize:"))
    for (size_t i(1);i<musage.size();++i) musage.front()+=" "+musage[i];
#ifdef USING__Threading
  msg_Info()<<" done ( "<<(musage.size()>0?musage.front():"?")<<", "
	    <<FormatTime(size_t(retime-rbtime))<<" / "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
#else
  msg_Info()<<" done ( "<<(musage.size()>0?musage.front():"?")<<", "
	    <<FormatTime(size_t(etime-btime))<<" )."<<std::endl;
#endif
  msg_Debugging()<<METHOD<<"(): Processes {\n";
  msg_Debugging()<<"  m_procs:\n";
  for (size_t i(0);i<m_procs.size();++i) 
    msg_Debugging()<<"    "<<m_procs[i]->Name()<<" -> "<<m_procs[i]<<"\n";
  msg_Debugging()<<"}\n";
  return res;
}

void Matrix_Element_Handler::BuildProcesses()
{
  Data_Reader read(" ",";","!","=");
  read.AddComment("#");
  read.AddWordSeparator("\t");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  // set color scheme
  read.SetTags(cls::ColorSchemeTags());
  cls::scheme cls((cls::scheme)read.GetValue<int>("COLOUR_SCHEME",1));
  read.SetTags(std::map<std::string,std::string>());
  // set helicity scheme
  read.SetTags(hls::HelicitySchemeTags());
  hls::scheme hls((hls::scheme)read.GetValue<int>("HELICITY_SCHEME",1));
  read.SetTags(std::map<std::string,std::string>());
  // set kfactor scheme
  std::string kfactor=read.GetValue<std::string>("KFACTOR","NO");
  // set scale scheme
  std::string scale=read.GetValue<std::string>("SCALES","METS");
  std::vector<std::string> helpsv;
  if (!read.VectorFromFile(helpsv,"COUPLINGS"))
    helpsv.push_back("Alpha_QCD 1");
  std::string coupling(MakeString(helpsv,0));
  // init processes
  msg_Info()<<METHOD<<"(): Looking for processes "<<std::flush;
  if (msg_LevelIsTracking()) msg_Info()<<"\n";
  std::vector<std::vector<std::string> > procdata;
  Data_Reader pread(" ",";","%",":");
  pread.AddComment("#");
  pread.AddWordSeparator("\t");
  pread.SetAddCommandLine(false);
  pread.SetInputPath(m_path);
  pread.SetInputFile(m_processfile);
  if (m_gens.size()>0 && !pread.MatrixFromFile(procdata,""))
    THROW(missing_input,"No data in "+m_path+m_processfile+"'.");
  for (size_t nf(0);nf<procdata.size();++nf) {
    std::vector<std::string> &cur(procdata[nf]);
    if (cur.size()<2) continue;
    if (cur[0]=="Process") {
      Process_Info pi;
      pi.m_scale=scale;
      pi.m_coupling=coupling;
      pi.m_kfactor=kfactor;
      pi.m_cls=cls;
      pi.m_hls=hls;
      std::string proc(MakeString(cur,1));
      size_t pos(proc.find("->"));
      if (pos==std::string::npos) continue;
      Processblock_Info pbi;
      std::string ini(proc.substr(0,pos));
      std::string fin(proc.substr(pos+2));
      std::vector<std::string> dectags;
      for (size_t ng(nf);ng<procdata.size();++ng) {
	std::vector<std::string> &cur(procdata[ng]);
	if (cur.size()<2) continue;
	if (cur[0]=="Decay") dectags.push_back(MakeString(cur,1));
	if (cur[0]=="DecayOS") dectags.push_back("Z"+MakeString(cur,1));
	if (cur[0]=="Scales") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vscale,nf);
	}
	if (cur[0]=="Couplings") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vcoupl,nf);
	}
	if (cur[0]=="KFactor") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vkfac,nf);
	}
	if (cur[0]=="CKKW") {
	  if (p_shower==NULL || p_shower->GetShower()==NULL)
	    THROW(fatal_error,"Invalid shower generator");
	  pi.m_ckkw=1;
	  pbi.m_gycut=cur[1];
	}
	if (cur[0]=="Y_Cut") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vycut,nf);
	}
	if (cur[0]=="Selector_File") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vsfile,nf);
	}
	if (cur[0]=="Order_EW") pi.m_oew=ToType<int>(cur[1]);
	if (cur[0]=="Order_QCD") pi.m_oqcd=ToType<int>(cur[1]);
	if (cur[0]=="Max_Order_EW") pi.m_maxoew=ToType<int>(cur[1]);
	if (cur[0]=="Max_Order_QCD") pi.m_maxoqcd=ToType<int>(cur[1]);
	if (cur[0]=="Max_N_Quarks") pi.m_nmaxq=ToType<int>(cur[1]);
	if (cur[0]=="Print_Graphs") pi.m_gpath=cur[1];
	if (cur[0]=="Presample_MC") pi.m_psmc=ToType<int>(cur[1]);
	if (cur[0]=="Enable_MHV") pi.m_amegicmhv=ToType<int>(cur[1]);	
	if (cur[0]=="Integration_Error") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vmaxerr,nf);
	}
	if (cur[0]=="Max_Epsilon") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vmaxeps,nf);
	}
	if (cur[0]=="Enhance_Factor") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vefac,nf);
	}
	if (cur[0]=="Enhance_Function") {
	  std::string cb(MakeString(cur,1));
	  ExtractMPvalues(cb,pbi.m_vefunc,nf);
	}
	if (cur[0]=="NLO_QCD_Part") {
	  pi.m_fi.m_nloqcdtype=ToType<nlo_type::code>(cur[1]);
	  pi.m_nlomode=m_nlomode;
	}
	if (cur[0]=="NLO_EW_Part") {
	  pi.m_fi.m_nloewtype=ToType<nlo_type::code>(cur[1]);
	  pi.m_nlomode=m_nlomode;
	}
	if (cur[0]=="Loop_Generator") pi.m_loopgenerator=cur[1];
        pi.p_gens=&m_gens;
	if (cur[0]=="End" && cur[1]=="process") break;
      }
      BuildSingleProcessList(pi,pbi,ini,fin,dectags);
      if (msg_LevelIsDebugging()) {
        msg_Indentation(4);
        msg_Out()<<m_procs.size()<<" process(es) found ..."<<std::endl;
        for (unsigned int i=0; i<m_procs.size(); ++i) {
          msg_Out()<<m_procs[i]->Name();
          if (m_procs[i]->IsGroup())
            msg_Out()<<" has subprocesses ...";
          msg_Out()<<std::endl;
        }
      }
    }
  }
}

void Matrix_Element_Handler::BuildDecays
(Subprocess_Info &ACFS,const std::vector<std::string> &dectags)
{
  for (size_t i(0);i<dectags.size();++i) {
    std::string dec(dectags[i]);
    int osf=0;
    if (dec[0]=='Z') {
      dec=dec.substr(1);
      osf=1;
    }
    size_t pos(dec.find("->"));
    if (pos==std::string::npos) continue;
    Subprocess_Info ACDIS, ACDFS;
    std::string ini(dec.substr(0,pos));
    std::string fin(dec.substr(pos+2));
    ExtractFlavours(ACDIS,ini);
    ExtractFlavours(ACDFS,fin);
    if (ACDIS.m_ps.empty() || ACDFS.m_ps.empty())
      THROW(fatal_error,"Wrong decay specification");
    Subprocess_Info &CIS(ACDIS.m_ps.front());
    size_t oldsize(ACFS.m_ps.size()), cdsize(ACDFS.m_ps.size());
    ACFS.m_ps.resize(oldsize*cdsize);
    for (size_t cfss(1);cfss<cdsize;++cfss)
      for (size_t acfsi(0);acfsi<oldsize;++acfsi)
	ACFS.m_ps[cfss*oldsize+acfsi]=
	  ACFS.m_ps[(cfss-1)*oldsize+acfsi];
    for (size_t acfsi(0);acfsi<oldsize;++acfsi) {
      for (size_t cfss(0);cfss<cdsize;++cfss) {
	Subprocess_Info &CFS(ACDFS.m_ps[cfss]);
	msg_Debugging()<<METHOD<<"(): Init decay {\n"<<CIS<<CFS<<"}\n";
	if (CIS.NExternal()!=1 || CFS.NExternal()<2)
	  THROW(fatal_error,"Wrong number of particles in decay");
	if (!ACFS.m_ps[cfss*oldsize+acfsi].AddDecay(CIS,CFS,osf))
	  THROW(fatal_error,"No match for decay "+dec);
      }
    }
  }
}

void Matrix_Element_Handler::BuildSingleProcessList
(Process_Info &pi,Processblock_Info &pbi,
 const std::string &ini,const std::string &fin,
 const std::vector<std::string> &dectags)
{
  Subprocess_Info AIS, AFS;
  ExtractFlavours(AIS,ini);
  ExtractFlavours(AFS,fin);
  std::vector<Process_Base*> procs;
  for (size_t fss(0);fss<AFS.m_ps.size();++fss) {
    Subprocess_Info ACFS;
    ACFS.m_ps.push_back(AFS.m_ps[fss]);
    BuildDecays(ACFS,dectags);
    for (size_t afsi(0);afsi<ACFS.m_ps.size();++afsi) {
      msg_Debugging()<<METHOD<<"(): Check N_max ("
		     <<fss<<"): {\n"<<ACFS.m_ps[afsi]<<"}\n";
      pi.m_fi.GetNMax(ACFS.m_ps[afsi]);
    }
  }
  for (size_t iss(0);iss<AIS.m_ps.size();++iss) {
    Subprocess_Info &IS(AIS.m_ps[iss]);
    for (size_t fss(0);fss<AFS.m_ps.size();++fss) {
      Subprocess_Info &FS(AFS.m_ps[fss]);
      msg_Debugging()<<METHOD<<"(): Init core ("<<iss<<","
		     <<fss<<"): {\n"<<IS<<FS<<"}\n";
      std::vector<Flavour> flavs;
      IS.GetExternal(flavs);
      if (flavs.size()>1) {
        if (!p_isr->CheckConsistency(&flavs.front())) {
          msg_Error()<<METHOD<<"(): Error in initialising isr ("
                     <<p_isr->Flav(0)<<" -> "<<flavs[0]<<") x ("
                     <<p_isr->Flav(1)<<" -> "<<flavs[1]
                     <<"). Ignoring process."<<std::endl;
          continue;
        }
      }
      Subprocess_Info ACFS;
      ACFS.m_ps.push_back(FS);
      BuildDecays(ACFS,dectags);
      for (size_t afsi(0);afsi<ACFS.m_ps.size();++afsi) {
	Subprocess_Info &CFS(ACFS.m_ps[afsi]);
	msg_Debugging()<<METHOD<<"(): Init process ("<<iss<<","
		       <<fss<<"): {\n"<<IS<<CFS<<"}\n";
	std::vector<Flavour> flavs;
	IS.GetExternal(flavs);
	CFS.GetExternal(flavs);
	size_t nis(IS.NExternal()), nfs(CFS.NExternal());
	double inisum=0.0, finsum=0.0, dd(0.0);
	for (size_t i(0);i<nis;++i) inisum+=flavs[i].Mass();
	for (size_t i(0);i<nfs;++i) finsum+=flavs[i+nis].Mass();
	if (inisum>=rpa.gen.Ecms() || finsum>=rpa.gen.Ecms()) continue;
	std::string pnid(CFS.MultiplicityTag()), ds;
	Process_Info cpi(pi);
	cpi.m_ii=IS;
	cpi.m_fi=CFS;
        cpi.m_fi.m_nloqcdtype=pi.m_fi.m_nloqcdtype;
        cpi.m_fi.m_nloewtype=pi.m_fi.m_nloewtype;
	cpi.m_fi.SetNMax(pi.m_fi);
	if (GetMPvalue(pbi.m_vscale,nfs,pnid,ds)) cpi.m_scale=ds;
	if (GetMPvalue(pbi.m_vcoupl,nfs,pnid,ds)) cpi.m_coupling=ds;
	if (GetMPvalue(pbi.m_vkfac,nfs,pnid,ds)) cpi.m_kfactor=ds;
	if (GetMPvalue(pbi.m_vsfile,nfs,pnid,ds)) cpi.m_selectorfile=ds;
	std::vector<Process_Base*> proc=InitializeProcess(cpi);
	for (size_t i(0);i<proc.size();i++) {
	  procs.push_back(proc[i]);
	  proc[i]->Integrator()->
	    SetISRThreshold(ATOOLS::Max(inisum,finsum));
	  if (GetMPvalue(pbi.m_vefunc,nfs,pnid,ds))
	    proc[i]->Integrator()->SetEnhanceFunction(ds);
	  if (GetMPvalue(pbi.m_vefac,nfs,pnid,dd))
	    proc[i]->Integrator()->SetEnhanceFactor(dd);
	  if (GetMPvalue(pbi.m_vmaxeps,nfs,pnid,dd))
	    proc[i]->Integrator()->SetMaxEpsilon(dd);
	  double maxerr(-1.0);
	  if (GetMPvalue(pbi.m_vmaxerr,nfs,pnid,dd)) maxerr=dd;
	  proc[i]->Integrator()->SetPSHandler
	    (new Phase_Space_Handler
	     (proc[i]->Integrator(),p_isr,p_beam,maxerr));
	}
      }
    }
  }
  if (pi.m_ckkw && rpa.gen.NumberOfEvents()==0)
    THROW(fatal_error,"Number of events cannot be zero in CKKW mode");
  for (size_t i(0);i<procs.size();++i) {
    Process_Info &cpi(procs[i]->Info());
    Selector_Key skey(NULL,new Data_Reader(),true);
    std::string sfile(cpi.m_selectorfile!=""?
		      cpi.m_selectorfile:m_selectorfile);
    while (sfile[sfile.length()-1]==' ') sfile.erase(sfile.length()-1,1);
    skey.ReadData(m_path,sfile);
    if (pi.m_ckkw&1) {
      std::vector<std::string> jfargs(1,pbi.m_gycut);
      GetMPvalue(pbi.m_vycut,cpi.m_fi.NExternal(),
		 cpi.m_fi.MultiplicityTag(),jfargs[0]);
      if (i==0) jfargs.push_back("LO");
      if (procs.size()>1) skey.SetData("METS",jfargs);
    }
    procs[i]->SetSelector(skey);
    if (pi.m_ckkw&1) {
      cpi.m_coupling="Alpha_QCD 1";
      if (procs.size()>1) {
	if (cpi.m_scale.rfind('}')==std::string::npos)
	  cpi.m_scale+="{MU_F2}";
	cpi.m_scale+="{"+p_shower->GetShower()->GetKT2("Q2_CUT")+"}";
      }
    }
    procs[i]->SetScale
      (Scale_Setter_Arguments(p_model,cpi.m_scale,cpi.m_coupling));
    procs[i]->SetKFactor
      (KFactor_Setter_Arguments(cpi.m_kfactor));
    procs[i]->SetShower(p_shower->GetShower());
  }
}

size_t Matrix_Element_Handler::ExtractFlavours(Subprocess_Info &info,std::string buffer)
{
  info.m_ps.resize(1);
  info.m_ps.front().m_ps.clear();
  while(true) {
    while (buffer.length()>0 && 
	   (buffer[0]==' ' || buffer[0]=='\t')) buffer.erase(0,1);
    if (buffer.length()==0) break;
    size_t pos(Min(buffer.find(' '),buffer.length()));
    std::string cur(buffer.substr(0,pos));
    buffer=buffer.substr(pos);
    pos=cur.find('(');
    std::string polid;
    if (pos!=std::string::npos) {
      polid=cur.substr(pos);
      cur=cur.substr(0,pos);
      polid=polid.substr(1,polid.find(')')-1);
    }
    pos=cur.find('[');
    std::string decid;
    if (pos!=std::string::npos) {
      decid=cur.substr(pos);
      cur=cur.substr(0,pos);
      decid=decid.substr(1,decid.find(']')-1);
    }
    int n(-1);
    pos=cur.find('{');
    if (pos!=std::string::npos) {
      std::string nid(cur.substr(pos));
      cur=cur.substr(0,pos);
      n=ToType<size_t>(nid.substr(1,nid.find('}')-1));
    }
    int kfc(ToType<int>(cur));
    Flavour cfl((kf_code)abs(kfc));
    if (kfc<0) cfl=cfl.Bar();
    if (n==-1) {
      for (size_t i(0);i<info.m_ps.size();++i)
	info.m_ps[i].m_ps.push_back(Subprocess_Info(cfl,decid,polid));
    }
    else {
      size_t oldsize(info.m_ps.size());
      info.m_ps.resize(oldsize*(n+1));
      for (int j(1);j<=n;++j) {
	for (size_t i(0);i<oldsize;++i) {
	  info.m_ps[j*oldsize+i]=info.m_ps[(j-1)*oldsize+i];
	  info.m_ps[j*oldsize+i].m_ps.push_back(Subprocess_Info(cfl,"",polid));
	}
      }
    }
  }
  return info.m_ps.back().m_ps.size();
}

namespace SHERPA {

  template <> double Matrix_Element_Handler::ExtractMPvalue(const std::string& str)
  {
    Algebra_Interpreter inter;
    inter.AddTag("E_CMS",ToString(rpa.gen.Ecms()));
    return ToType<double>(inter.Interprete(str));
  }

  template <> std::string Matrix_Element_Handler::ExtractMPvalue(const std::string& str)
  {
    return str;
  }

  template <typename Type>
  void Matrix_Element_Handler::AddMPvalue(std::string lstr,std::string rstr,const Type &val,
			  std::map<std::string,std::pair<int,Type> >& dv,
			  const int nfs,const int &priority)
  {
    if (rstr.length()==0) {
      if (nfs==0 && 
	  (dv.find(lstr)==dv.end() || dv[lstr].first>priority)) {
	msg_Debugging()<<METHOD<<"(): adding '"<<val
		       <<"' {"<<lstr<<"}("<<priority<<")\n";
	dv[lstr]=std::pair<int,Type>(priority,val);
      }
      return;
    }
    size_t pos(rstr.find('-')), ltp(rstr.find('['));
    if (pos==std::string::npos || ltp<pos-1) {
      if (ltp!=std::string::npos) {
	size_t rtp(rstr.find(']',ltp));
	AddMPvalue(lstr+rstr.substr(0,rtp+1),rstr.substr(rtp+1),val,dv,
		   nfs-ToType<int>(rstr.substr(ltp+1,rtp-ltp-1)),priority);
	return;
      }
      AddMPvalue(lstr+rstr,"",val,dv,nfs-ToType<int>(rstr),priority);
      return;
    }
    std::string rlstr(rstr.substr(0,pos)), rrstr(rstr.substr(pos+1)), rmstr;
    if (pos>0 && ltp==pos-1) {
      rmstr="]";
      rrstr=rrstr.substr(1);
    }
    for (int i(0);i<=nfs;++i)
      AddMPvalue(lstr+rlstr+ToString(i)+rmstr,rrstr,val,dv,nfs-i,priority);
  }

  template void Matrix_Element_Handler::AddMPvalue
  (std::string lstr,std::string rstr,const double &val,
   std::map<std::string,std::pair<int,double> >& dv,
   const int nfs,const int &priority);
  template void Matrix_Element_Handler::AddMPvalue
  (std::string lstr,std::string rstr,const std::string &val,
   std::map<std::string,std::pair<int,std::string> >& dv,
   const int nfs,const int &priority);

  template <typename Type>
  bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,Type> >& dv,
   const int nfs,const std::string &pnid,Type &rv)
  {
    std::map<std::string,std::pair<int,Type> > cdv(dv);
    for (typename std::map<std::string,std::pair<int,Type> >::const_iterator 
	   dit(dv.begin());dit!=dv.end();++dit) { 
      AddMPvalue<Type>("",dit->first,dit->second.second,
		       dv,nfs,dit->second.first);
    }
    if (dv.find(pnid)!=dv.end()) {
      rv=dv[pnid].second;
      return true;
    }
    std::string nfstag(ToString(nfs));
    if (dv.find(nfstag)!=dv.end()) {
      rv=dv[nfstag].second;
      return true;
    }
    return false;
  }

  template bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,double> >& dv,
   const int nfs,const std::string &pnid,double &rv);
  template bool Matrix_Element_Handler::GetMPvalue
  (std::map<std::string,std::pair<int,std::string> >& dv,
   const int nfs,const std::string &pnid,std::string &rv);

  template <typename Type>
  void Matrix_Element_Handler::ExtractMPvalues(std::string& str,std::map
			       <std::string,std::pair<int,Type> >& dv,
			       const int &priority)
  {
    int position;
    position = str.find("{");
    while (position>0 && str[position-1]!=' ' && str[position-1]!='\t')
      position=str.find('{',position+1);
    if (position==-1) {
      msg_Debugging()<<METHOD<<"(): adding '"<<str
		     <<"' {-}("<<priority<<")\n";
      dv["-"]=std::pair<int,Type>(priority,ExtractMPvalue<Type>(str));
      return;
    }
    std::string hstr = str.substr(0,position);
    Type value = ExtractMPvalue<Type>(hstr);
    str = str.substr(position+1,str.length()-position-2);
    do {
      position = str.find(",");
      if (position>-1) {
	hstr = str.substr(0,position);
	str = str.substr(position+1);
      }
      else hstr=str;
      if (hstr.length()>0) {
	msg_Debugging()<<METHOD<<"(): adding '"<<value
		       <<"' {"<<hstr<<"}("<<priority<<")\n";
	dv[hstr]=std::pair<int,Type>(priority,value);
      }
    } while (position>-1);
  }

  template void Matrix_Element_Handler::ExtractMPvalues
  (std::string& str,std::map<std::string,std::pair<int,double> >& dv,
   const int &priority);
  template void Matrix_Element_Handler::ExtractMPvalues
  (std::string& str,std::map<std::string,std::pair<int,std::string> >& dv,
   const int &priority);

}

std::string Matrix_Element_Handler::MakeString
(const std::vector<std::string> &in,const size_t &first)
{
  std::string out(in.size()>first?in[first]:"");
  for (size_t i(first+1);i<in.size();++i) out+=" "+in[i];
  return out;
}


