#include "Analysis_Handler.H"

#include "Shell_Tools.H"
#include "Data_Reader.H"
#include "Run_Parameter.H"

#ifdef PROFILE__all
#define PROFILE__Analysis_Handler
#endif
#ifdef PROFILE__Analysis_Handler
#include "prof.hh"
#else 
#define PROFILE_HERE
#define PROFILE_LOCAL(LOCALNAME)
#endif

#define Getter_Function ANALYSIS::Primitive_Observable_Base::Getter_Function
using namespace ANALYSIS;

size_t Analysis_Handler::s_maxanalyses=100;

Analysis_Handler::Analysis_Handler() {}

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

String_Matrix
Analysis_Handler::FindArguments(const String_Matrix &strings,
				const size_t starty,const size_t startx)
{
  size_t j=0, open=0;
  String_Matrix result;
  if (strings[starty].size()>startx) j=startx;
  for (size_t i=starty;i<strings.size();++i) {
    result.push_back(std::vector<std::string>(strings[i].size()-j));
    for (size_t k=0;j<strings[i].size();++j,++k) {
      result.back()[k]=strings[i][j];
      size_t opos=result.back()[k].find("{");
      if (opos!=std::string::npos) {
	++open;
	if (open==1) {
	  result.back()[k]=result.back()[k].substr(opos+1);
	  if (result.back()[k].length()==0) {
	    --k;
	    result.back().resize(result.back().size()-1);
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
	  if (k==0 && result.back()[0].length()==0) result.resize(result.size()-1);
	  return result;
	}
      }
    }
    if (open==0) break;
    j=0;
  }  
  return result;
}

void Analysis_Handler::ShowSyntax(const size_t i)
{
  if (!ATOOLS::msg.LevelIsInfo() || i==0) return;
  ATOOLS::msg.Out()<<"Analysis_Handler::ShowSyntax(): {\n\n"
		   <<"   ..     -  mandatory variable\n"
		   <<"   ..|..  -  mandatory selection\n"
		   <<"   [..]   -  optional variable\n"
		   <<"   -> ..  -  depends on\n\n"
                   <<"   list   -  particle list specifier, \n"
		   <<"             user defined lists can be created using e.g. Trigger \n"
                   <<"             predefined lists include:\n"
		   <<"              FinalState, Parton, ChargedParticle, Hadron, ChargedHadron,\n"
		   <<"              NeutralHadron, ChargePion, ChargedKaon, ProtonAntiproton \n\n"
		   <<"   BEGIN_ANALYSIS {\n\n"
		   <<"   LEVEL      [ME]|[MI]|[Shower]|[Hadron]\n\n"
		   <<"   PATH_PIECE path\n\n";
  Getter_Function::PrintGetterInfo(ATOOLS::msg.Out(),10);
  ATOOLS::msg.Out()<<"\n   } END_ANALYSIS\n\n"
		   <<"}"<<std::endl;
}

bool Analysis_Handler::ReadIn()
{
  msg_Info()<<"Analysis_Handler::ReadIn(): {\n";
  bool success=false;
  std::vector<std::string> helpsv;
  std::vector<std::vector<std::string> > helpsvv;
  ATOOLS::Data_Reader reader;
  reader.SetVectorType(reader.VHorizontal);
  reader.SetMatrixType(reader.MTransposed);
  reader.SetInputPath(InputPath());
  reader.SetInputFile(InputFile());
  reader.AddComment("!");
  reader.AddComment("%");
  reader.AddIgnore("+");
  reader.SetFileBegin("BEGIN_ANALYSIS");
  reader.SetFileEnd("END_ANALYSIS");
  for (size_t i=0;i<s_maxanalyses;++i) {
    reader.SetOccurrence(i);
    reader.RereadInFile();
    if (!reader.VectorFromFile(helpsv,"LEVEL")) break;
    bool split=false, trigger=false;
    int mode=ANALYSIS::fill_all|ANALYSIS::splitt_jetseeds;
    for (size_t j=0;j<helpsv.size();++j) {
      if (split) mode=mode|ANALYSIS::splitt_phase;
      else split=true;
      if (helpsv[j]=="ME") mode=mode|ANALYSIS::do_me;
      else if (helpsv[j]=="MI") mode=mode|ANALYSIS::do_mi;
      else if (helpsv[j]=="Shower") mode=mode|ANALYSIS::do_shower;
      else if (helpsv[j]=="Hadron") mode=mode|ANALYSIS::do_hadron;
      else {
	ATOOLS::msg.Error()<<"Analysis_Handler::ReadIn(): "
			   <<"Invalid analysis mode '"<<helpsv[j]
			   <<"'"<<std::endl;
	continue;
      }
      success=true;
      msg_Info()<<"   new Primitive_Analysis(\""<<helpsv[j]<<"\")\n";
      msg_Tracking()<<"   new Primitive_Analysis(\""<<helpsv[j]<<"\") {\n";
      m_analyses.push_back(new ANALYSIS::Primitive_Analysis(helpsv[j],mode));
      std::string outpath;
      if (!reader.ReadFromFile(outpath,"PATH_PIECE")) outpath="";
      m_analyses.back()->SetOutputPath(outpath);
      reader.MatrixFromFile(helpsvv,"");
      for (size_t k=0;k<helpsvv.size();++k) {
	String_Matrix mat=FindArguments(helpsvv,k,1);
	ANALYSIS::Primitive_Observable_Base *observable = 
	  Getter_Function::GetObject(helpsvv[k][0],mat);
	if (observable!=NULL) {
	  m_analyses.back()->AddObservable(observable);
	  if (helpsvv[k][0]=="Trigger") trigger=true;
	  if (ATOOLS::msg.LevelIsTracking()) {
	    ATOOLS::msg.Out()<<"      new Primitive_Observable_Base(\""
			     <<helpsvv[k][0]<<"\",";
	    for (size_t i=0;i<mat.size();++i) {
	      ATOOLS::msg.Out()<<"{"<<(mat[i].size()>0?mat[i][0]:"");
	      for (size_t j=1;j<mat[i].size();++j) 
		ATOOLS::msg.Out()<<","<<mat[i][j];
	      ATOOLS::msg.Out()<<"}";
	    }
	    ATOOLS::msg.Out()<<")\n";
	  }
	}
      }
      if (!trigger) {
	ANALYSIS::Primitive_Observable_Base *observable = 
	  Getter_Function::GetObject("Trigger",String_Matrix());
	m_analyses.back()->AddObservable(observable);
      }
      msg_Tracking()<<"   }\n";
    }
  }
  msg_Info()<<"}"<<std::endl;
  if (success) ATOOLS::Exception_Handler::AddTerminatorObject(this);
  return success;
}

void Analysis_Handler::DoAnalysis(const ATOOLS::Blob_List *bloblist,
				  const double weight)
{
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) (*ait)->DoAnalysis(bloblist,weight); 
}

void Analysis_Handler::Clear()
{ 
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) (*ait)->ClearAllData(); 
}

void Analysis_Handler::PrepareTerminate()
{
  Finish();
}

void Analysis_Handler::Finish(const std::string &path)
{
  if (OutputPath()[OutputPath().length()-1]=='/') {
    if (!ATOOLS::MakeDir(OutputPath(),448)) {
      ATOOLS::msg.Error()<<"Analysis_Handler::Finish(..): "
			 <<"Cannot create directory '"<<OutputPath()
			 <<"'."<<std::endl; 
      SetOutputPath(path);
      if (!ATOOLS::MakeDir(OutputPath(),448)) {
	ATOOLS::msg.Error()<<"Analysis_Handler::Finish(..): "
			   <<"Cannot create directory '"<<OutputPath()
			   <<"'."<<std::endl; 
	SetOutputPath(ATOOLS::rpa.gen.Variable("SHERPA_RUN_PATH"));
      }
    }
  }
  msg_Info()<<"Analysis_Handler::Finish(..): {\n";
  for (Analyses_Vector::const_iterator ait=m_analyses.begin();
       ait!=m_analyses.end();++ait) {
    msg_Info()<<"   Writing to '"<<OutputPath()<<(*ait)->OutputPath()<<"'."<<std::endl; 
    (*ait)->FinishAnalysis(OutputPath()); 
  }
  msg_Info()<<"}"<<std::endl;
  if (m_analyses.size()) ATOOLS::Exception_Handler::RemoveTerminatorObject(this);
}

