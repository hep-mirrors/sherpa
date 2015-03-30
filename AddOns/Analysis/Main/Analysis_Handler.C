#include "AddOns/Analysis/Main/Analysis_Handler.H"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Math/Variable.H"
#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Org/My_MPI.H"

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

size_t Analysis_Handler::s_maxanalyses=100;

Analysis_Handler::Analysis_Handler():
  Analysis_Interface("Internal"),
  m_weighted(0), m_write(false)
{
  if(s_kftable.find(kf_bjet)==s_kftable.end()) // if not initialized yet
  s_kftable[kf_bjet] = new
    Particle_Info(kf_bjet,0.,0.,0,1, 2,1,1,1,0,"bj","bj",1,1);
}

Analysis_Handler::~Analysis_Handler()
{
  Clean();
  exh->RemoveTesterObject(this);
}

void Analysis_Handler::Clean()
{
  while (m_analyses.size()>0) {
    delete m_analyses.back();
    m_analyses.pop_back();
  }
}

Argument_Matrix
Analysis_Handler::FindArguments(const Argument_Matrix &strings,
				size_t &starty,size_t &startx)
{
  size_t j=0, open=0;
  Argument_Matrix result;
  if (strings[starty].size()>startx) j=startx;
  for (size_t i=starty;i<strings.size();++i) {
    result.push_back(std::vector<std::string>(strings[i].size()-j));
    for (size_t k=0;j<strings[i].size();++j,++k) {
      result.back()[k]=strings[i][j];
      size_t opos=result.back()[k].find("{");
      if (opos!=std::string::npos) {
	++open;
	if (open==1) {
	  if (opos+1<result.back()[k].size()) 
	    result.back()[k]=result.back()[k].substr(opos+1);
	  else 
	    result.back()[k]="";
	  if (result.back()[k].length()==0) {
	    --k;
	    result.back().pop_back();
	    if (result.back().size()==0) result.pop_back();
	    continue;
	  }
	}
      }
      if (open>0) {
	size_t cpos=result.back()[k].find("}");
	if (cpos!=std::string::npos) --open;
	if (open==0) {
	  result.back()[k]=result.back()[k].substr(0,cpos);
	  result.back().resize(k+1);
	  if (k==0 && result.back()[0].length()==0) 
	    result.resize(result.size()-1);
	  return result;
	}
      }
    }
    if (open==0) break;
    j=0;
  }  
  return result;
}

void Analysis_Handler::ShowSyntax(const int i)
{
  ATOOLS::Variable_Base<double>::ShowVariables(i);
  ATOOLS::Particle_Qualifier_Base::ShowQualifiers(i);
  if (!msg_LevelIsInfo() || i==0) return;
  msg_Out()<<"Analysis_Handler::ShowSyntax(): {\n\n"
	   <<"   ..     -  mandatory variable\n"
	   <<"   ..|..  -  mandatory selection\n"
	   <<"   [..]   -  optional variable\n"
	   <<"   -> ..  -  depends on\n\n"
	   <<"   list   -  particle list specifier, \n\n"
	   <<"   BEGIN_ANALYSIS {\n\n"
	   <<"   LEVEL      [ME]|[Shower]|[Hadron]\n\n"
	   <<"   PATH_PIECE path\n\n"
	   <<"   // observable listing\n\n";
  Observable_Getter_Function::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n   // detector/trigger & tools listing\n\n";
  Object_Getter_Function::PrintGetterInfo(msg->Out(),15);
  msg_Out()<<"\n   } END_ANALYSIS\n\n"
	   <<"}"<<std::endl;
}

bool Analysis_Handler::Init()
{
  msg_Info()<<"Analysis_Handler::ReadIn(): {\n";
  bool success=false;
  std::vector<std::string> helpsv;
  std::vector<std::vector<std::string> > helpsvv;
  Data_Reader reader(" ",";","//");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.SetInputPath(InputPath());
  std::string infile(InputFile());
  if (infile.find('|')!=std::string::npos)
    infile=infile.substr(0,infile.find('|'));
  reader.SetInputFile(infile);
  reader.AddComment("#");
  reader.SetFileBegin("BEGIN_ANALYSIS");
  reader.SetFileEnd("END_ANALYSIS");
  reader.AddFileBegin("BEGIN_ANALYSIS{");
  reader.AddFileEnd("}END_ANALYSIS");
  for (size_t i=0;i<s_maxanalyses;++i) {
    reader.SetOccurrence(i);
    reader.RescanInFile();
    if (!reader.VectorFromFile(helpsv,"LEVEL")) break;
    bool split=false;
    int mode=ANALYSIS::fill_all|
      ANALYSIS::splitt_jetseeds|ANALYSIS::split_sh;
    for (size_t j=0;j<helpsv.size();++j) {
      if (split) mode=mode|ANALYSIS::splitt_phase;
      else split=true;
      if (helpsv[j].find("MENLO")!=std::string::npos) 
	mode=mode|ANALYSIS::do_menlo;
      else if (helpsv[j].find("ME")!=std::string::npos) 
	mode=mode|ANALYSIS::do_me;
      else if (helpsv[j].find("MI")!=std::string::npos) 
	mode=mode|ANALYSIS::do_mi;
      else if (helpsv[j].find("Shower")!=std::string::npos) 
	mode=mode|ANALYSIS::do_shower;
      else if (helpsv[j].find("Hadron")!=std::string::npos) 
	mode=mode|ANALYSIS::do_hadron;
      else {
	msg_Error()<<"Analysis_Handler::ReadIn(): "
		   <<"Invalid analysis mode '"<<helpsv[j]
		   <<"'"<<std::endl;
	continue;
      }
    }
    success=true;
    std::string outpath;
    if (!reader.ReadFromFile(outpath,"PATH_PIECE")) outpath="";
    msg_Info()<<"   new Primitive_Analysis(\""<<outpath<<"\") -> "<<helpsv[0];
    for (size_t j=1;j<helpsv.size();++j) msg_Info()<<","<<helpsv[j];
    msg_Info()<<"\n";
    msg_Tracking()<<"   new Primitive_Analysis(..) {\n";
    mode=mode|m_weighted;
    m_analyses.push_back(new Primitive_Analysis(this,ToString(i),mode));
    m_analyses.back()->SetOutputPath(outpath);
    std::string maxjettag;
    if (!reader.ReadFromFile(maxjettag,"NMAX_JETS")) maxjettag="";
    m_analyses.back()->SetMaxJetTag(maxjettag);
    int splitjetconts;
    if (!reader.ReadFromFile(splitjetconts,"JETCONTS")) splitjetconts=1;
    m_analyses.back()->SetSplitJetConts(splitjetconts);
    reader.MatrixFromFile(helpsvv,"");
    Argument_Matrix arguments(helpsvv);
    for (size_t k=0;k<helpsvv.size();++k) {
      if (arguments[k].size()>0) {
	if (arguments[k][0]=="{" || arguments[k][0]=="}") continue;
      }
      size_t col=1;
      Argument_Matrix mat=FindArguments(arguments,k,col);
      ANALYSIS::Primitive_Observable_Base *observable = 
	Observable_Getter_Function::GetObject
	(arguments[k][0],mat(m_analyses.back()));
      if (observable!=NULL) {
	m_analyses.back()->AddObject(observable);
	if (msg_LevelIsTracking()) {
	  msg_Out()<<"      new Primitive_Observable_Base(\""
		   <<arguments[k][0]<<"\",";
	  for (size_t i=0;i<mat.size();++i) {
	    msg_Out()<<"{"<<(mat[i].size()>0?mat[i][0]:"");
	    for (size_t j=1;j<mat[i].size();++j) 
	      msg_Out()<<","<<mat[i][j];
	    msg_Out()<<"}";
	  }
	  msg_Out()<<")\n";
	}
      }
      ANALYSIS::Analysis_Object *object = 
	Object_Getter_Function::GetObject
	(arguments[k][0],mat(m_analyses.back()));
      if (object!=NULL) {
	m_analyses.back()->AddObject(object);
	if (msg_LevelIsTracking()) {
	  msg_Out()<<"      new Analysis_Object(\""
		   <<arguments[k][0]<<"\",";
	  for (size_t i=0;i<mat.size();++i) {
	    msg_Out()<<"{"<<(mat[i].size()>0?mat[i][0]:"");
	    for (size_t j=1;j<mat[i].size();++j) 
	      msg_Out()<<","<<mat[i][j];
	    msg_Out()<<"}";
	  }
	  msg_Out()<<")\n";
	}
      }
    }
    msg_Tracking()<<"   }\n";
  }
  msg_Info()<<"}"<<std::endl;
  if (success) {
    m_write=true;
    exh->AddTesterObject(this);
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

bool Analysis_Handler::ApproveTerminate()
{
  if (rpa->gen.BatchMode()&1) m_write=false;
  return true;
}

bool Analysis_Handler::WriteOut()
{
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
    (*ait)->FinishAnalysis(OutputPath());
    (*ait)->RestoreAnalysis();
  }
  return true;
}

bool Analysis_Handler::Finish()
{
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
    (*ait)->FinishAnalysis(OutputPath());
    (*ait)->RestoreAnalysis();
  }
  msg_Info()<<"}"<<std::endl;
  if (m_analyses.size()) {
    exh->RemoveTesterObject(this);
  }
  return true;
}

bool Analysis_Handler::Run(ATOOLS::Blob_List *const bl)
{
  Blob *sp(bl->FindFirst(btp::Signal_Process));
  // if no signal process (i.e. hadrons execs etc.), assume weight to be one
  if (!sp) {
    DoAnalysis(bl,1.);
    return true;
  }
  Blob_Data_Base *xs((*sp)["Weight"]);
  if (xs==NULL) THROW(fatal_error,"No weight information");
  double wgt(xs->Get<double>());
  DoAnalysis(bl,wgt);
  return true;
}

DECLARE_GETTER(Analysis_Handler,"Internal",
	       Analysis_Interface,Analysis_Arguments);

Analysis_Interface *ATOOLS::Getter
<Analysis_Interface,Analysis_Arguments,Analysis_Handler>::
operator()(const Analysis_Arguments &args) const
{
  Analysis_Handler *analysis(new ANALYSIS::Analysis_Handler());
  analysis->SetInputPath(args.m_inpath);
  analysis->SetInputFile(args.m_infile);
  analysis->SetOutputPath(args.m_outpath);
  return analysis;
}

void ATOOLS::Getter<Analysis_Interface,Analysis_Arguments,
		    Analysis_Handler>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"internal analysis interface";
}

