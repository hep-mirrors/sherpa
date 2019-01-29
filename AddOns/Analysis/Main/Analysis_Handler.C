#include "AddOns/Analysis/Main/Analysis_Handler.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Variable.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/KF_Table.H"

#ifdef PROFILE__all
#define PROFILE__Analysis_Handler
#endif
#ifdef PROFILE__Analysis_Handler
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

#define Observable_Getter_Function \
  ANALYSIS::Primitive_Observable_Base::Observable_Getter_Function

using namespace ANALYSIS;
using namespace SHERPA;
using namespace ATOOLS;

Analysis_Handler::Analysis_Handler():
  Analysis_Interface("Internal"),
  m_weighted(0), m_write(false)
{
  if(s_kftable.find(kf_bjet)==s_kftable.end()) // if not initialized yet
  s_kftable[kf_bjet] = new
    Particle_Info(kf_bjet,0.,0.,0,1, 2,1,1,1,0,"bj","bj","bj","bj",1,1);
}

Analysis_Handler::~Analysis_Handler()
{
  Clean();
}

void Analysis_Handler::Clean()
{
  while (m_analyses.size()>0) {
    delete m_analyses.back();
    m_analyses.pop_back();
  }
}

void Analysis_Handler::ShowSyntax(const int i)
{
  ATOOLS::Variable_Base<double>::ShowVariables(i);
  ATOOLS::Particle_Qualifier_Base::ShowQualifiers(i);
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<METHOD;
  msg_Out()<<"(): {\n\n"
	   <<"   You can give analyses as elements of the INTERNAL_ANALYSES\n"
	   <<"   yaml sequence in your yaml configuration. Each analysis is a\n"
	   <<"   yaml mapping with the following keys (among others):\n"
	   <<"\n"
	   <<"   LEVEL: [...]\n"
	   <<"   PATH_PIECE: ...\n"
	   <<"   OBSERVABLES: [...]\n"
	   <<"   ANALYSES_OBJECTS: [...]\n"
	   <<"\n"
	   <<"   LEVEL can be a list of the following keys: MENLO, ME, MI\n"
	   <<"   Shower, Hadron."
	   <<"\n"
	   <<"   Each observable/analysis object is itself a mapping with a\n"
	   <<"   single key-value pair. The key gives the name of the\n"
	   <<"   observable/object and the value is either a scalar, mapping\n"
	   <<"   or a sequence, giving more details on its evaluation.\n"
	   <<"   In the following we list the possible keys and the\n"
	   <<"   corresponding value syntax for each observable and\n"
	   <<"   analysis object.\n"
	   <<"\n"
	   <<"   Observables:\n"
	   <<"\n";
  Observable_Getter_Function::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n\n   Analysis objects:\n\n";
  Object_Getter_Function::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"}"<<std::endl;
}

bool Analysis_Handler::Init()
{
  msg_Info()<<"Analysis_Handler::ReadIn(): {\n";
  bool success=false;

  Scoped_Settings analyses_settings{
    Settings::GetMainSettings()["INTERNAL_ANALYSES"] };
  analyses_settings.DeclareVectorSettingsWithEmptyDefault({ "LEVEL" });

  for (auto& s : analyses_settings.GetItems()) {
    const auto levels = s["LEVEL"].GetVector<std::string>();
    if (levels.empty())
      break;
    auto split = false;
    const auto splitsh = s["SPLITSH"].SetDefault(true).Get<bool>();
    int mode = ANALYSIS::fill_all
               | ANALYSIS::split_vars
               | ANALYSIS::splitt_jetseeds
               | (splitsh ? ANALYSIS::split_sh : 0);
    for (const auto& level : levels) {
      if (split)
        mode = mode | ANALYSIS::splitt_phase;
      if (level == "MENLO")
        mode = mode | ANALYSIS::do_menlo;
      else if (level == "ME")
        mode = mode | ANALYSIS::do_me;
      else if (level == "MI")
        mode = mode | ANALYSIS::do_mi;
      else if (level == "Shower")
        mode = mode | ANALYSIS::do_shower;
      else if (level == "Hadron")
        mode = mode | ANALYSIS::do_hadron;
      else {
        msg_Error()
          << "Analysis_Handler::ReadIn(): "
          << "Invalid analysis mode '" << level << "'" << std::endl;
        continue;
      }
    }
    success  = true;
    const auto outpath = s["PATH_PIECE"].SetDefault("").Get<std::string>();
    msg_Info() << "   new Primitive_Analysis(\"" << outpath << "\") -> " << levels[0];
    for (size_t j{ 1 }; j < levels.size(); ++j) msg_Info() << ", " << levels[j];
    msg_Info() << "\n";
    msg_Tracking() << "   new Primitive_Analysis(..) {\n";
    mode = mode | m_weighted;
    m_analyses.push_back(new Primitive_Analysis(this, ToString(s.GetIndex()), mode));
    m_analyses.back()->SetOutputPath(outpath);
    const auto usedb = s["USE_DB"].SetDefault(false).Get<bool>();
    m_analyses.back()->SetUseDB(usedb);
    const auto& maxjettag = s["NMAX_JETS"].SetDefault("").Get<std::string>();
    m_analyses.back()->SetMaxJetTag(maxjettag);
    const auto& splitjetconts = s["JETCONTS"].SetDefault(1).Get<int>();
    m_analyses.back()->SetSplitJetConts(splitjetconts);
    auto obssettings = s["OBSERVABLES"];
    for (auto& singleobssettings : obssettings.GetItems()) {
      const auto& obsnames = singleobssettings.GetKeys();
      if (obsnames.size() != 1)
        THROW(fatal_error,
              "Each observable setting must be a single key-value pair.");
      const auto& obsname = obsnames.front();
      Analysis_Key anakey{
        singleobssettings[obsname], m_analyses.back() };
      ANALYSIS::Primitive_Observable_Base* observable{
        Observable_Getter_Function::GetObject(obsname, anakey) };
      if (observable != nullptr) {
        m_analyses.back()->AddObject(observable);
      }
    }
    auto anasettings = s["ANALYSES_OBJECTS"];
    for (auto& singleanasettings : anasettings.GetItems()) {
      const auto& ananames = singleanasettings.GetKeys();
      if (ananames.size() != 1)
        THROW(fatal_error,
              "Each observable setting must be a single key-value pair.");
      const auto& ananame = ananames.front();
      Analysis_Key anakey{
        singleanasettings[ananame], m_analyses.back() };
      ANALYSIS::Analysis_Object* object{
        Object_Getter_Function::GetObject(ananame, anakey) };
      if (object != nullptr) {
        m_analyses.back()->AddObject(object);
      }
    }
    msg_Tracking()<<"   }\n";
  }
  msg_Info()<<"}"<<std::endl;
  if (success) {
    m_write=true;
  }
  return true;
}

void Analysis_Handler::DoAnalysis(const ATOOLS::Blob_List *bloblist,
				  const double weight)
{
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) (*ait)->DoAnalysis(bloblist,weight); 
}

void Analysis_Handler::CleanUp()
{ 
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) (*ait)->ClearAllData(); 
}

bool Analysis_Handler::WriteOut()
{
  if (rpa->gen.BatchMode()&1) m_write=false;
  if (!m_write) return true;
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
  if (OutputPath()[OutputPath().length()-1]=='/') {
    if (!MakeDir(OutputPath())) {
      msg_Error()<<"Analysis_Handler::Finish(..): "
		 <<"Cannot create directory '"<<OutputPath()
		 <<"'."<<std::endl; 
    }
  }
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) {
    (*ait)->FinishAnalysis(OutputPath(),0);
    (*ait)->RestoreAnalysis();
  }
  return true;
}

bool Analysis_Handler::Finish()
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()==0)
#endif
  if (OutputPath()[OutputPath().length()-1]=='/') {
    if (!MakeDir(OutputPath())) {
      msg_Error()<<"Analysis_Handler::Finish(..): "
		 <<"Cannot create directory '"<<OutputPath()
		 <<"'."<<std::endl; 
    }
  }
  msg_Info()<<"Analysis_Handler::Finish(..): {\n";
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) {
    msg_Info()<<"   Writing to '"<<OutputPath()<<(*ait)->OutputPath()
	      <<"'."<<std::endl; 
    (*ait)->FinishAnalysis(OutputPath(),0);
    (*ait)->RestoreAnalysis();
  }
  msg_Info()<<"}"<<std::endl;
  return true;
}

bool Analysis_Handler::Run(ATOOLS::Blob_List *const bl)
{
  DoAnalysis(bl,1.0);
  return true;
}

DECLARE_GETTER(Analysis_Handler,"Internal",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,Analysis_Handler>::
operator()(const Analysis_Arguments &args) const
{
  Analysis_Handler *analysis(new ANALYSIS::Analysis_Handler());
  analysis->SetOutputPath(args.m_outpath);
  return analysis;
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,
		    Analysis_Handler>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"internal analysis interface";
}

