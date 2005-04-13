#include "Grid_Creator.H"

#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Data_Reader.H"
#include "Data_Writer.H"
#include "Blob.H"
#include "Run_Parameter.H"

#ifdef PROFILE__all
#define PROFILE__Grid_Creator
#endif
#ifdef PROFILE__Grid_Creator
#include "prof.hh"
#else
#define PROFILE_HERE
#endif

using namespace AMISIC;

Grid_Creator::Grid_Creator(Amisic_Histogram_Map *histograms,
			   EXTRAXS::XS_Group *const processes):
  p_histograms(histograms),
  p_processes(processes),
  m_xsextension("_xs.dat"),
  m_mcextension("MC"),
  m_datatag("[x,w,w2,max,n] = "),
  m_events(0)
{
  if (p_processes==NULL) {
    THROW(fatal_error,"Process handler is not initialized");
  }
  if (!CollectProcesses(p_processes)) {
    THROW(fatal_error,"Process handler does not own any process");
  }
}

Grid_Creator::~Grid_Creator()
{
}

bool Grid_Creator::CollectProcesses(EXTRAXS::XS_Base *const process)
{
  if (process->Size()==0) return false;
  if ((*process)[0]==process) {
    (*p_histograms)[process->Name()] = new Amisic_Histogram<double>(4);
    return true;
  }
  for (size_t i=0;i<process->Size();++i) {
    if (!CollectProcesses((*process)[i])) return false;
  }
  return true;
}

std::string Grid_Creator::MakeString(std::vector<std::string> input) const
{
  for (unsigned int i=1;i<input.size();++i) {
    input[0]+=std::string(" ")+input[i];
  }
  return input.size()>0 ? input[0] : "";
}

bool Grid_Creator::ReadInArguments(std::string tempifile,
				   std::string tempipath)
{
  if (tempipath!="") SetInputPath(tempipath);
  if (tempifile!="") SetInputFile(tempifile);
  if (!CheckInputFile()) return false;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  reader->SetInputFile(InputPath()+InputFile());
  reader->SetVectorType(reader->VHorizontal);
  std::vector<std::string> helps;
  if (!reader->VectorFromFile(helps,"X_VARIABLE")) m_gridxvariable="p_\\perp";
  else m_gridxvariable=MakeString(helps);
  if (!reader->VectorFromFile(helps,"Y_VARIABLE")) m_gridyvariable="";
  else m_gridyvariable=MakeString(helps);
  if (m_gridxmin==0.0) 
    m_gridxmin=sqrt(ATOOLS::Max(p_processes->ISR()->PDF(0)->Q2Min(),
				p_processes->ISR()->PDF(1)->Q2Min()));
  m_gridxmin=ATOOLS::Max(m_gridxmin,1.e-3);
  if (!reader->ReadFromFile(m_griddeltax,"GRID_DELTA_X")) 
    m_griddeltax=(log(m_gridxmax)-log(m_gridxmin))/250.;
  double helpd;
  if (!reader->ReadFromFile(helpd,"INITIAL_EVENTS")) helpd=0;
  m_initevents=(long unsigned)helpd;
  if (!reader->ReadFromFile(helpd,"MAX_EVENTS")) helpd=100000;
  m_maxevents=(long unsigned)helpd;
  if (!reader->ReadFromFile(m_binerror,"GRID_ERROR")) m_binerror=0.05;
  if (!reader->ReadFromFile(m_outputlevel,"GRID_CREATOR_OUTPUT")) 
    m_outputlevel=0;
  if (!reader->ReadFromFile(m_gridxscaling,"HISTO_X_SCALING")) 
    m_gridxscaling="Log_B_10";
  if (!reader->ReadFromFile(m_gridyscaling,"HISTO_Y_SCALING")) 
    m_gridyscaling="Id";
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    hit->second->XAxis()->SetVariable(m_gridxvariable);
    hit->second->YAxis()->SetVariable(m_gridyvariable);
    hit->second->XAxis()->SetScaling(m_gridxscaling);
    hit->second->YAxis()->SetScaling(m_gridyscaling);
    Amisic_Histogram_Type::Axis_Type *axis=hit->second->XAxis();
    if (!hit->second->Initialize(m_gridxmin,m_gridxmax,
				 abs((int)(((*axis)(m_gridxmax)-
					    (*axis)(m_gridxmin))/
					   m_griddeltax)))) 
      THROW(critical_error,"Cannot initialize histogram.");
  }
  delete reader;
  p_xaxis=p_histograms->begin()->second->XAxis();
  p_yaxis=p_histograms->begin()->second->YAxis();
  p_variable=p_xaxis->Variable();
  return true;
}

void Grid_Creator::Clear()
{
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    hit->second->Initialize(m_gridxmin,m_gridxmax,
			    abs((int)(((*p_xaxis)(m_gridxmax)-
				       (*p_xaxis)(m_gridxmin))/
				      m_griddeltax)));
  }    
}

bool Grid_Creator::ReadInGrid()
{
  PROFILE_HERE;
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    if (hit->second->ReadIn(OutputPath()+hit->first+m_xsextension,
			    m_datatag)) {
      if (hit->second->XMin()-m_gridxmin>m_gridxmin*1.0e-7 ||
	  m_gridxmax-hit->second->XMax()>m_gridxmax*1.0e-7 ||
	  hit->second->Entries()<m_initevents) {
	Clear();
	return false;
      }
    }
    else {
      Clear();
      return false;
    }
  }
  return true;
}

bool Grid_Creator::InitializeCalculation()
{
  int helpi=0;
  m_criterion=p_xaxis->Variable()->SelectorID();
  std::vector<ATOOLS::Flavour> flavours(1,(ATOOLS::kf::jet));
  p_processes->SelectorData()->
    SetData(m_criterion,flavours,helpi,m_gridxmin,m_gridxmax);
  p_processes->ResetSelector(p_processes->SelectorData());
  p_processes->Reset();
  p_processes->CalculateTotalXSec(OutputPath()+OutputFile()
				  +MCExtension(),true);
  return true;
}

bool Grid_Creator::UpdateHistogram(EXTRAXS::XS_Base *const process)
{
  if (process->Size()==0) return false;
  if ((*process)[0]==process) {
    process->Parent()->SetSelected(process);
    ATOOLS::Blob_Data_Base *xsdata=process->Parent()->SameWeightedEvent();
    PHASIC::Weight_Info info=xsdata->Get<PHASIC::Weight_Info>();
    delete xsdata;
    Amisic_Histogram_Type *histo=(*p_histograms)[process->Name()];
    const ATOOLS::Vec4D *p=process->Momenta();
    double value=(*p_variable)(&p[0]);
    for (size_t i=1;i<4;++i) value=ATOOLS::Max(value,(*p_variable)(&p[i]));
    histo->Add(value,info.weight);
    for (size_t i=1;i<info.ntrial;++i) histo->Add(value,0.0); 
    return true;
  }
  for (size_t i=0;i<process->Size();++i) {
    if (!UpdateHistogram((*process)[i])) return false;
  }
  return true;
}

bool Grid_Creator::CreateOptimizedGrid()
{
  msg_Info()<<"Grid_Creator::CreateOptimizedGrid(): "
	    <<"Optimizing grid for MI.\n";
  double starttime=ATOOLS::rpa.gen.Timer().UserTime();
  msg_Info()<<ATOOLS::tm::curoff;
  for (;m_events<m_maxevents;++m_events) {
    if (!UpdateHistogram(p_processes)) return false;
    if ((m_events%(m_maxevents/100))==0 && m_events>0) {
      double diff=ATOOLS::rpa.gen.Timer().UserTime()-starttime;
      msg_Info()<<"   "<<((100*m_events)/m_maxevents)<<" % ( "
		<<int(diff)<<" s elapsed / "
		<<int((m_maxevents-m_events)/(double)m_events*diff)
		<<" s left / "<<int(m_maxevents/(double)m_events*diff)
		<<" s total )   "<<ATOOLS::bm::cr<<std::flush;
    }
  }
  msg_Info()<<ATOOLS::tm::curon;
  return true;
}

bool Grid_Creator::CreateInitialGrid()
{
  msg_Info()<<"Grid_Creator::CreateInitialGrid(): "
	    <<"Initializing grid for MI.\n";
  double starttime=ATOOLS::rpa.gen.Timer().UserTime();
  msg_Info()<<ATOOLS::tm::curoff;
  for (;m_events<m_initevents;++m_events) {
    if (!UpdateHistogram(p_processes)) return false;
    if ((m_events%(m_initevents/100))==0 && m_events>0) {
      double diff=ATOOLS::rpa.gen.Timer().UserTime()-starttime;
      msg_Info()<<"   "<<((100*m_events)/m_initevents)<<" % ( "
		<<int(diff)<<" s elapsed / "
		<<int((m_initevents-m_events)/(double)m_events*diff)
		<<" s left / "<<int(m_initevents/(double)m_events*diff)
		<<" s total )   "<<ATOOLS::bm::cr<<std::flush;
    }
  }
  msg_Info()<<ATOOLS::tm::curon;
  return true;
}

bool Grid_Creator::WriteOutGrid(std::vector<std::string> addcomments) 
{
  PROFILE_HERE;
  bool success=true;
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    hit->second->Finish();
    if (!hit->second->WriteOut(OutputPath()+hit->first
			       +m_xsextension,m_datatag,
			       addcomments)) success=false;
  }
  return success;
}

bool Grid_Creator::CreateGrid()
{
  bool success=true;
  int formerlevel=ATOOLS::msg.Level();
  msg_Info()<<"Grid_Creator::CreateGrid(): "
	    <<"Calculating grid {"<<std::endl;
  ATOOLS::msg.SetLevel(m_outputlevel);
  if (!InitializeCalculation()) {
    ATOOLS::msg.Error()<<"Grid_Creator_Base::CreateGrid(..): "
		       <<"Initialization failed! Abort."<<std::endl;
    return false;
  }
  if (!CreateInitialGrid()) {
    ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		     <<"Initial grid creation failed."<<std::endl;
    success=false;
  }
  if (!CreateOptimizedGrid()) {
    ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		     <<"Sorry, grid cannot be optimized."<<std::endl;
    success=false;
  }
  if (!WriteOutGrid()) {
    ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		     <<"Sorry, grid cannot be written to '"
		     <<OutputFile()<<"'"<<std::endl;
    success=false;
  }
  ATOOLS::msg.SetLevel(formerlevel);
  msg_Info()<<"\n}"<<std::endl;
  return success;
}

