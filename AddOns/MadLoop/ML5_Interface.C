#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#define ML5_BOOST_TO_CMS

inline int mp(const int id,const int i)
{ return i+id*4; }

inline void GetMom(double *m,const int n,const ATOOLS::Vec4D &p)
{ m[mp(n,0)]=p[0]; for (int i(1);i<4;++i) m[mp(n,i)]=p[i]; }

extern "C" {

  inline void MakeFortranString
  (char *output,std::string input,unsigned int length)
  {
    for (unsigned int i=0;i<length;++i) 
      output[i]=(char)32;
    for (size_t j=0;j<input.length();++j) 
      output[j]=(char)input[j];
  }

  void setpara2_(char *name);
  void update_as_param2_(double *mur,double *as);
  void setmadlooppath_(char *name);
  void setparamlog_(int *log);
  void printout_();

}

using namespace PHASIC;
using namespace ATOOLS;

namespace ML5 {

  class ML5_Interface: public PHASIC::ME_Generator_Base,
		       public ATOOLS::Terminator_Object {
  private:

    static std::string s_path, s_model;

    static int s_init, s_mode;
    static size_t s_pid;

  public :

    // constructor
    ML5_Interface(): ME_Generator_Base("ML5")
    {
      rpa->gen.AddCitation(1,"NLO ME from \\cite{Hirschi:2011pa}.");
      exh->AddTerminatorObject(this);
    }

    ~ML5_Interface()
    {
      PrepareTerminate();
    }

    void PrepareTerminate();

    // member functions
    bool Initialize(const std::string &path,const std::string &file,
		    MODEL::Model_Base *const model,
		    BEAM::Beam_Spectra_Handler *const beam,
		    PDF::ISR_Handler *const isr)
    {
      Data_Reader read(" ",";","!","=");
      read.SetInputPath(path);
      read.SetInputFile(file);
      std::string filename(file);
      if (filename.find('.')!=std::string::npos)
	filename=filename.substr(0,filename.find('.'));
      s_path=read.GetValue<std::string>("ML5_LIBPATH","Process/OLP/"+filename);
      s_model=read.GetValue<std::string>
	("ML5_MODEL",Flavour(kf_b).IsMassive()?"loop_sm":"loop_sm-no_b_mass");
      s_mode=read.GetValue<int>("ML5_MODE",0);
      s_init=!FileExists(s_path+".db");
      My_In_File::OpenDB(s_path+"/");
#ifdef USING__MPI
      if (MPI::COMM_WORLD.Get_rank()==0) {
#endif
      if (ML5_Interface::Init()) {
	Remove(s_path+".mg5");
	std::ofstream card((s_path+".mg5").c_str(),std::ios::app);
	if (s_model.substr(0,s_model.find('-'))=="loop_sm")
	  card<<"set complex_mass_scheme True\n";
	card<<"import model "<<s_model<<"\n";
      }
      std::string fname(rpa->gen.Variable("SHERPA_RUN_PATH"));
      fname+="/"+s_path+"/SubProcesses/MadLoop5_resources";
      if (model->Name()=="SM" && DirectoryExists(fname)) {
	std::ofstream param((fname+"/param_card.dat").c_str());
	param.precision(12);
	param<<"Block LOOP\n";
	param<<"    1 "<<Flavour(kf_Z).Mass()<<"\n";
	std::string mymodel(s_model.substr(0,s_model.find('-')));
	if (mymodel=="loop_sm") {
	  param<<"Block SMINPUTS\n";
	  param<<"    1 "<<1.0/model->ScalarConstant("alpha_QED")<<"\n";
	  param<<"    3 "<<model->ScalarConstant("alpha_S")<<"\n";
	}
	else if (mymodel=="loop_qcd_qed_sm") {
	  param<<"Block SMINPUTS\n";
	  param<<"    1 "<<1.0/model->ScalarConstant("alpha_QED")<<"\n";
	  param<<"    2 "<<model->ScalarConstant("alpha_S")<<"\n";
	}
	else {
	  msg_Error()<<METHOD<<"(): "<<om::red
		     <<"Unknown model. Please update param card.\n"<<om::reset;
	}
	param<<"Block MASS\n";
	for (KF_Table::const_iterator kit(s_kftable.begin());
	     kit!=s_kftable.end();++kit) {
	  if (kit->first==0) continue;
	  if (kit->first>25) break;
	  param<<"    "<<kit->first<<" "<<Flavour(kit->first).Mass()<<"\n";
	}
	param<<"Block YUKAWA\n";
	for (KF_Table::const_iterator kit(s_kftable.begin());
	     kit!=s_kftable.end();++kit) {
	  if (Flavour(kit->first).Mass()==0.0) continue;
	  if (kit->first>22) break;
	  param<<"    "<<kit->first<<" "<<Flavour(kit->first).Yuk()<<"\n";
	}
	for (KF_Table::const_iterator kit(s_kftable.begin());
	     kit!=s_kftable.end();++kit) {
	  if (kit->first==0) continue;
	  if (kit->first>25) break;
	  if (kit->second->m_width || kit->first==23 || kit->first==24)
	    param<<"DECAY "<<kit->first<<" "<<Flavour(kit->first).Width()<<"\n";
	}
	param<<"Block QNUMBERS 82\n";
	param<<"    1 0\n";
	param<<"    2 1\n";
	param<<"    3 8\n";
	param<<"    4 1\n";
      }
#ifdef USING__MPI
      }
      int dummy;
      MPI::COMM_WORLD.Bcast(&dummy,1,MPI::INT,0);
#endif  
      s_loader->AddPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+s_path+"/lib");
      s_loader->AddPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+s_path+"/lib/collier_lib");
      s_loader->AddPath(rpa->gen.Variable("SHERPA_RUN_PATH")+"/"+s_path+"/lib/ninja_lib");
      s_loader->LoadLibrary("collier");
      s_loader->LoadLibrary("ninja");
      s_loader->LoadLibrary("cts");
      s_loader->LoadLibrary("iregi");
      s_loader->LoadLibrary("ML5_DHELAS");
      s_loader->LoadLibrary("ML5_MODEL");
      return true;
    }

    PHASIC::Process_Base *InitializeProcess
    (const PHASIC::Process_Info &pi, bool add) { return NULL; }
    int  PerformTests() { return 1; }
    bool NewLibraries() { return s_init; }

    void SetClusterDefinitions(PDF::Cluster_Definitions_Base *const defs) {}

    ATOOLS::Cluster_Amplitude *ClusterConfiguration
    (PHASIC::Process_Base *const proc,const size_t &mode)
    { return NULL; }

    inline static std::string Path() { return s_path; }

    inline static int Init() { return s_init; }
    inline static int Mode() { return s_mode; }

    inline static size_t &PId() { return s_pid; }

  }; // end of class ML5_Interface
 
  std::string ML5_Interface::s_path="OLP";
  std::string ML5_Interface::s_model="loop";

  int ML5_Interface::s_init=0;
  int ML5_Interface::s_mode=0;

  size_t ML5_Interface::s_pid=0;

  typedef void (*ME_Function)(double *,double *,double *,double *,int *);
  typedef void (*MESO_Function)(int *);
  typedef void (*Path_Function)(char *);
  typedef void (*Void_Function)(void);

  class ML5_Process: public PHASIC::Virtual_ME2_Base {
  protected:
    ME_Function p_me;
    double *p_p, *p_res, *p_prec, m_prec;
    size_t m_id;
    std::string m_libname;
  public:
    ML5_Process(const PHASIC::Process_Info& pi,
		const ATOOLS::Flavour_Vector& flavs,int mode):
      Virtual_ME2_Base(pi,flavs), p_me(NULL),
      p_res(NULL), p_prec(NULL), m_prec(-1.0)
    {
      m_mode=1;
      std::string cn(ML5_Interface::Path()), cpn;
      std::string pn(Process_Base::GenerateName(pi.m_ii,pi.m_fi));
      if (ML5_Interface::Init()) {
#ifdef USING__MPI
	if (MPI::COMM_WORLD.Get_rank())
	  THROW(fatal_error,"Initialization not possible in MPI mode");
#endif
	size_t nin(pi.m_ii.m_ps.size());
	std::ofstream card((cn+".mg5").c_str(),std::ios::app);
	card<<"add process ";
	for (size_t i(0);i<nin;++i)
	  card<<(long int)(flavs[i])<<" ";
	card<<"> ";
	for (size_t i(nin);i<flavs.size();++i)
	  card<<(long int)(flavs[i])<<" ";
	int loqcd=pi.m_fi.m_nloqcdtype==nlo_type::loop;
	int loew=pi.m_fi.m_nloewtype==nlo_type::loop;
	if (pi.m_maxcpl[0]!=99) card<<"QCD="<<(pi.m_maxcpl[0]-loqcd)<<" ";
	if (pi.m_maxcpl[1]!=99) card<<"QED="<<(pi.m_maxcpl[1]-loew)<<" ";
	if (loqcd) card<<"[virt=QCD]";
	if (loew) card<<"[virt=QED]";
	if (pi.m_maxcpl[0]!=99) card<<" QCD="<<2*pi.m_maxcpl[0];
	if (pi.m_maxcpl[1]!=99) card<<" QED="<<2*pi.m_maxcpl[1];
	card<<" @"<<ML5_Interface::PId()<<"\n";
	My_In_File::ExecDB(cn+"/","begin");
	My_Out_File mapfile(cn+"/"+pn+".map");
	if (!mapfile.Open()) THROW(fatal_error,"Cannot write map file");
	*mapfile<<pn<<" "<<ML5_Interface::PId()<<"\n";
	mapfile.Close();
	My_In_File::ExecDB(cn+"/","commit");
	ML5_Interface::PId()++;
      }
      My_In_File mapfile(cn+"/"+pn+".map");
      if (!mapfile.Open()) THROW(fatal_error,"Cannot read map file");
      *mapfile>>cpn>>m_id;
      if (cpn!=pn) msg_Info()<<METHOD<<"():"<<om::red<<" Process '"<<pn
			     <<"' is mapped onto '"<<cpn<<"'!\n"<<om::reset;
      msg_Info()<<"!";
      p_p = new double[4*flavs.size()];
      std::string bp=rpa->gen.Variable("SHERPA_RUN_PATH");
      bp+="/"+cn+"/SubProcesses/MadLoop5_resources";
      m_libname="ML5_P"+ToString(m_id);
      void *module(s_loader->LoadLibrary(m_libname));
      if (module) {
	static bool init(false);
	if (!init) {
	  init=true;
	  int log=0;
	  char name[512];
	  MakeFortranString(name,bp+"/param_card.dat",512);
	  setparamlog_(&log);
	  setpara2_(name);
	  MakeFortranString(name,bp,512);
	  setmadlooppath_(name);
	}
	int onoff(1);
	((MESO_Function)s_loader->GetLibraryFunction
	 ("","ml5_"+ToString(m_id)+"_force_stability_check_",module))(&onoff);
	int info;
	((MESO_Function)s_loader->GetLibraryFunction
	 ("","ml5_"+ToString(m_id)+"_get_nsqso_loop_",module))(&info);
	p_res = new double[4*(1+info)];
	p_prec = new double[1+info];
	if (msg_LevelIsDebugging()) printout_();
	if (mode&1) s_loader->UnloadLibrary(m_libname,module);
	else p_me = (ME_Function)s_loader->GetLibraryFunction
	       ("","ml5_"+ToString(m_id)+"_sloopmatrix_thres_",module);
      }
    }

    ~ML5_Process()
    {
      delete [] p_p;
      if (p_res) delete [] p_res;
      if (p_prec) delete [] p_prec;
    }

    void Calc(const ATOOLS::Vec4D_Vector &ip)
    {
      void *module(NULL);
      if (p_me==NULL) {
	module=s_loader->LoadLibrary(m_libname);
	if (module==NULL) THROW(normal_exit,"Missing loop library");
	p_me = (ME_Function)s_loader->GetLibraryFunction
	  ("","ml5_"+ToString(m_id)+"_sloopmatrix_thres_",module);
	if (p_me==NULL) THROW(normal_exit,"Missing loop ME");
      }
#ifndef ML5_BOOST_TO_CMS
      const ATOOLS::Vec4D_Vector &p(ip);
#else
      Poincare cms(ip[0]+ip[1]);
      ATOOLS::Vec4D_Vector p(ip);
      for (size_t n(0);n<p.size();++n) cms.Boost(p[n]);
#endif
      for (size_t n(0);n<p.size();++n) GetMom(p_p,n,p[n]);
      double mur(sqrt(m_mur2)), as((*MODEL::as)(m_mur2));
      update_as_param2_(&mur,&as);
      int retcode;
      p_me(p_p,p_res,&m_prec,p_prec,&retcode);
      if (module) {
	s_loader->UnloadLibrary(m_libname,module);
	p_me=NULL;
      }
      if (retcode/100==4) {
	msg_Error()<<METHOD<<"(): Unstable point {\n";
	msg_Error()<<"  Process "<<m_flavs<<"\n";
	std::cout.precision(16);
	for (size_t i(0);i<p.size();++i)
	  msg_Error()<<"  p_lab["<<i<<"]=Vec4D"<<p[i]<<";\n";
	msg_Error()<<"}"<<std::endl;
	if (ML5_Interface::Mode()&1) abort();
      }
      double norm((2.0*M_PI)/(*MODEL::as)(m_mur2));
      m_res.Finite()=p_res[1]*norm/m_norm;
      m_res.IR()=p_res[2]*norm/m_norm;
      m_res.IR2()=p_res[3]*norm/m_norm;
      m_born=p_res[0];
    }

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.0*M_PI;
    }

  };// end of class ML5_Process

  class ML5_LoopSquared: public PHASIC::Tree_ME2_Base {
  protected:
    ME_Function p_me;
    double *p_p, *p_res, *p_prec, m_prec;
    size_t m_id;
    std::string m_libname;
  public:

    ML5_LoopSquared(const Process_Info& pi,
	     const Flavour_Vector& flavs):
      Tree_ME2_Base(pi,flavs), p_me(NULL),
      p_res(NULL), p_prec(NULL), m_prec(-1.0)
    {
      std::string cn(ML5_Interface::Path()), cpn;
      std::string pn(Process_Base::GenerateName(pi.m_ii,pi.m_fi));
      if (ML5_Interface::Init()) {
#ifdef USING__MPI
	if (MPI::COMM_WORLD.Get_rank())
	  THROW(fatal_error,"Initialization not possible in MPI mode");
#endif
	size_t nin(pi.m_ii.m_ps.size());
	std::ofstream card((cn+".mg5").c_str(),std::ios::app);
	card<<"add process ";
	for (size_t i(0);i<nin;++i)
	  card<<(long int)(flavs[i])<<" ";
	card<<"> ";
	for (size_t i(nin);i<flavs.size();++i)
	  card<<(long int)(flavs[i])<<" ";
	int loqcd=pi.m_fi.m_nloqcdtype==nlo_type::loop;
	int loew=pi.m_fi.m_nloewtype==nlo_type::loop;
	if (pi.m_maxcpl[0]!=99) card<<"QCD="<<(pi.m_maxcpl[0]-loqcd)<<" ";
	if (pi.m_maxcpl[1]!=99) card<<"QED="<<(pi.m_maxcpl[1]-loew)<<" ";
	card<<"[virt=QCD]";
	if (pi.m_maxcpl[0]!=99) card<<" QCD="<<2*pi.m_maxcpl[0];
	if (pi.m_maxcpl[1]!=99) card<<" QED="<<2*pi.m_maxcpl[1];
	card<<" @"<<ML5_Interface::PId()<<"\n";
	My_In_File::ExecDB(cn+"/","begin");
	My_Out_File mapfile(cn+"/"+pn+".map");
	if (!mapfile.Open()) THROW(fatal_error,"Cannot write map file");
	*mapfile<<pn<<" "<<ML5_Interface::PId()<<"\n";
	mapfile.Close();
	My_In_File::ExecDB(cn+"/","commit");
	ML5_Interface::PId()++;
      }
      My_In_File mapfile(cn+"/"+pn+".map");
      if (!mapfile.Open()) THROW(fatal_error,"Cannot read map file");
      *mapfile>>cpn>>m_id;
      if (cpn!=pn) msg_Info()<<METHOD<<"():"<<om::red<<" Process '"<<pn
			     <<"' is mapped onto '"<<cpn<<"'!\n"<<om::reset;
      msg_Info()<<"!";
      p_p = new double[4*flavs.size()];
      std::string bp=rpa->gen.Variable("SHERPA_RUN_PATH");
      bp+="/"+cn+"/SubProcesses/MadLoop5_resources";
      m_libname="ML5_P"+ToString(m_id);
      void *module(s_loader->LoadLibrary(m_libname));
      if (module) {
	static bool init(false);
	if (!init) {
	  init=true;
	  int log=0;
	  char name[512];
	  MakeFortranString(name,bp+"/param_card.dat",512);
	  setparamlog_(&log);
 	  setpara2_(name);
	  MakeFortranString(name,bp,512);
	  setmadlooppath_(name);
	}
	int onoff(1);
	((MESO_Function)s_loader->GetLibraryFunction
	 ("","ml5_"+ToString(m_id)+"_force_stability_check_",module))(&onoff);
	int info;
	((MESO_Function)s_loader->GetLibraryFunction
	 ("","ml5_"+ToString(m_id)+"_get_nsqso_loop_",module))(&info);
	p_res = new double[4*(1+info)];
	p_prec = new double[1+info];
	if (msg_LevelIsDebugging()) printout_();
	p_me = (ME_Function)s_loader->GetLibraryFunction
	  ("","ml5_"+ToString(m_id)+"_sloopmatrix_thres_",module);
      }
    }

    ~ML5_LoopSquared()
    {
      delete [] p_p;
      if (p_res) delete [] p_res;
      if (p_prec) delete [] p_prec;
    }

    double Calc(const Vec4D_Vector &ip)
    {
      void *module(NULL);
      if (p_me==NULL) {
	if (ML5_Interface::Init()) return 1.0;
	module=s_loader->LoadLibrary(m_libname);
	if (module==NULL) THROW(normal_exit,"Missing loop library");
	p_me = (ME_Function)s_loader->GetLibraryFunction
	  ("","ml5_"+ToString(m_id)+"_sloopmatrix_thres_",module);
	if (p_me==NULL) THROW(normal_exit,"Missing loop ME");
      }
      if (ML5_Interface::Init()) return 1.0;
#ifndef ML5_BOOST_TO_CMS
      const ATOOLS::Vec4D_Vector &p(ip);
#else
      Poincare cms(ip[0]+ip[1]);
      ATOOLS::Vec4D_Vector p(ip);
      for (size_t n(0);n<p.size();++n) cms.Boost(p[n]);
#endif
      for (size_t n(0);n<p.size();++n) GetMom(p_p,n,p[n]);
      double mur(p_aqcd->Scale());
      if (mur<0.0) mur=rpa->gen.Ecms();
      else mur=sqrt(mur);
      double as((*MODEL::as)(mur*mur));
      update_as_param2_(&mur,&as);
      int retcode;
      p_me(p_p,p_res,&m_prec,p_prec,&retcode);
      if (module) {
	s_loader->UnloadLibrary(m_libname,module);
	p_me=NULL;
      }
      if (retcode%10==0) {
	msg_Error()<<METHOD<<"(): Unstable point {\n";
	msg_Error()<<"  Process "<<m_flavs<<"\n";
	std::cout.precision(16);
	for (size_t i(0);i<p.size();++i)
	  msg_Error()<<"  p_lab["<<i<<"]=Vec4D"<<p[i]<<";\n";
	msg_Error()<<"}"<<std::endl;
	if (ML5_Interface::Mode()&1) abort();
      }
      return p_res[1]/m_norm;
    }

    int OrderQCD(const int &id) { return 99; }

    int OrderEW(const int &id) { return 99; }

  };// end of class ML5_Tree

  void ML5_Interface::PrepareTerminate()
  {
    My_In_File::CloseDB(s_path+"/");
    if (ML5_Interface::Init()) {
      std::ofstream card((s_path+".mg5").c_str(),std::ios::app);
      card<<"output "<<s_path<<"\n";
      std::string loc(MADLOOP_PREFIX);
      loc+="/Template/loop_material/StandAlone/Cards";
      std::ifstream mlpin((loc+"/MadLoopParams.dat").c_str());
      std::ofstream mlpout((s_path+"_ML5Params.dat").c_str());
      bool tnow=false, fnow=false;
      std::string line;
      for (getline(mlpin,line);mlpin.good();getline(mlpin,line)) {
	if (tnow) {
	  tnow=false;
	  line=".TRUE.";
	}
	if (fnow) {
	  fnow=false;
	  line=".FALSE.";
	}
	if (line.find("#DoubleCheckHelicityFilter")!=std::string::npos) fnow=true;
	if (line.find("#WriteOutFilters")!=std::string::npos) tnow=true;
	if (line.find("#UseLoopFilter")!=std::string::npos) tnow=true;
	mlpout<<line<<"\n";
      }
      mlpout.close();
      mlpin.close();
#ifdef ML5_CALL_LAUNCH
      card<<"launch\n0\n";
      for (size_t i(0);i<s_pid;++i) card<<"n\n";
#endif
      std::ofstream script("makeloops");
      script<<"#!/bin/bash\ncpwd=$PWD\n";
      script<<"test -d "<<s_path<<" && exit 1\n";
      script<<MADLOOP_PREFIX<<"/bin/mg5_aMC < "<<s_path<<".mg5\n";
      script<<rpa->gen.Variable("SHERPA_SHARE_PATH")
	    <<"/sconsloops "<<s_path<<"\n";
      script<<"scons install\n";
      ChMod("makeloops",0755);
      msg_Out()<<om::red<<"Run './makeloops' to build loop library"
	       <<om::reset<<std::endl;
    }
  }

  class ML5D_Process {};

}// end of namespace ML5

using namespace ML5;

DECLARE_GETTER(ML5_Interface,"ML5",ME_Generator_Base,ME_Generator_Key);
ME_Generator_Base *ATOOLS::Getter
<ME_Generator_Base,ME_Generator_Key,ML5_Interface>::
operator()(const ME_Generator_Key &key) const
{
  return new ML5_Interface();
}
void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,ML5_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Interface to the MadLoop loop ME generator"; 
}

DECLARE_VIRTUALME2_GETTER(ML5_Process,"ML5_Process")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,ML5_Process>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="ML5") return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  return new ML5_Process(pi,fl,0);
}

DECLARE_VIRTUALME2_GETTER(ML5D_Process,"ML5D_Process")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,ML5D_Process>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="ML5D") return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  return new ML5_Process(pi,fl,1);
}

DECLARE_TREEME2_GETTER(ML5_Process,"ML5_Process")
Tree_ME2_Base *ATOOLS::Getter
<Tree_ME2_Base,Process_Info,ML5_Process>::
operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator.find("ML5")!=0) return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo) return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  return new ML5_LoopSquared(pi,fl);
}
