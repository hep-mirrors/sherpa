#include "Grid_Creator.H"

#include "Phase_Space_Handler.H"
#include "ISR_Handler.H"
#include "Data_Writer.H"
#include "Blob.H"
#include "Run_Parameter.H"

using namespace AMISIC;

Grid_Creator::Grid_Creator(Grid_Handler_Map *gridhandlers,
			   Grid_Handler_Map *maxhandlers,
			   EXTRAXS::XS_Group *const processes):
  p_gridhandlers(gridhandlers),
  p_maxhandlers(maxhandlers),
  p_processes(processes),
  m_xsextension("_xs.dat"),
  m_maxextension("_max.dat"),
  m_storemax(true)
{
  if (p_processes==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			    "Process handler is not initialized",
			    "Grid_Creator","Grid_Creator"));
  }
  if (!CollectProcesses(p_processes)) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			    "Process handler does not own any process",
			    "Grid_Creator","Grid_Creator"));
  }
}

Grid_Creator::~Grid_Creator()
{
  while (m_histograms.size()>0) {
    delete m_histograms.begin()->second;
    m_histograms.erase(m_histograms.begin());
  }
}

bool Grid_Creator::CollectProcesses(EXTRAXS::XS_Base *const process)
{
  if (process->Size()==0) return false;
  if ((*process)[0]==process) {
    m_histograms[process->Name()] = new Amisic_Histogram<double>();
    (*p_gridhandlers)[process->Name()] = new Grid_Handler_Type();
    if (m_storemax) 
      (*p_maxhandlers)[process->Name()] = new Grid_Handler_Type();
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
  if (!reader->VectorFromFile(helps,"X_VARIABLE")) 
    m_gridxvariable="p_\\perp";
  else m_gridxvariable=MakeString(helps);
  if (!reader->VectorFromFile(helps,"Y_VARIABLE")) 
    m_gridyvariable="\\frac{d\\sigma}{dp_\\perp}";
  else m_gridyvariable=MakeString(helps);
  if (!reader->ReadFromFile(m_gridxmin,"GRID_X_MIN")) 
    m_gridxmin=sqrt(ATOOLS::Max(p_processes->ISR()->PDF(0)->Q2Min(),
				p_processes->ISR()->PDF(1)->Q2Min()));
  if (!reader->ReadFromFile(m_gridxmax,"GRID_X_MAX")) 
    m_gridxmax=ATOOLS::rpa.gen.Ecms()/2.0;
  if (!reader->ReadFromFile(m_griddeltax,"GRID_DELTA_X")) 
    m_griddeltax=(log(m_gridxmax)-log(m_gridxmin))/1000.;
  double helpd;
  if (reader->ReadFromFile(helpd,"INITIAL_EVENTS")) 
    m_initevents=(long unsigned)helpd;
  else m_initevents=100000;
  if (reader->ReadFromFile(helpd,"MAX_EVENTS")) 
    m_maxevents=(long unsigned)helpd;
  else m_maxevents=1000000;
  if (!reader->ReadFromFile(m_binerror,"GRID_ERROR")) m_binerror=0.05;
  if (!reader->ReadFromFile(m_outputlevel,"GRID_CREATOR_OUTPUT")) {
    m_outputlevel=ATOOLS::msg.Level();
  }
  if (!reader->ReadFromFile(m_gridxscaling,
			    "HISTO_X_SCALING")) m_gridxscaling="Log_B_10";
  if (!reader->ReadFromFile(m_gridyscaling,
			    "HISTO_Y_SCALING")) m_gridyscaling="Id";
  for (Amisic_Histogram_Map::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();++hit) {
    hit->second->XAxis()->SetVariable(m_gridxvariable);
    hit->second->YAxis()->SetVariable(m_gridyvariable);
    hit->second->XAxis()->SetScaling(m_gridxscaling);
    hit->second->YAxis()->SetScaling(m_gridyscaling);
    Amisic_Histogram_Type::Axis_Type *axis=hit->second->XAxis();
    hit->second->Initialize(m_gridxmin,m_gridxmax,
			    abs((int)(((*axis)(m_gridxmax)-
				       (*axis)(m_gridxmin))/m_griddeltax)));
  }
  if (!reader->ReadFromFile(m_gridxscaling,
			    "GRID_X_SCALING")) m_gridxscaling="Log_B_10";
  if (!reader->ReadFromFile(m_gridyscaling,
			    "GRID_Y_SCALING")) m_gridyscaling="Log_B_10";
  for (Grid_Handler_Map::iterator git=p_gridhandlers->begin();
       git!=p_gridhandlers->end();++git) {
    git->second->Grid()->XAxis()->SetVariable(m_gridxvariable);
    git->second->Grid()->YAxis()->SetVariable(m_gridyvariable);
    git->second->Grid()->XAxis()->SetScaling(m_gridxscaling);
    git->second->Grid()->YAxis()->SetScaling(m_gridyscaling);
  }
  if (m_storemax) {
    for (Grid_Handler_Map::iterator mit=p_maxhandlers->begin();
	 mit!=p_maxhandlers->end();++mit) {
      mit->second->Grid()->XAxis()->SetVariable(m_gridxvariable);
      mit->second->Grid()->YAxis()->SetVariable(m_gridyvariable);
      mit->second->Grid()->XAxis()->SetScaling(m_gridxscaling);
      mit->second->Grid()->YAxis()->SetScaling(m_gridyscaling);
    }
  }
  delete reader;
  p_xaxis=m_histograms.begin()->second->XAxis();
  p_yaxis=m_histograms.begin()->second->YAxis();
  return true;
}

bool Grid_Creator::ReadInGrid()
{
  Grid_Handler_Map::iterator hit=p_gridhandlers->begin();
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader("=",";","#");
  reader->SetInputPath(OutputPath());
  reader->SetInputFile(hit->first+m_xsextension);
  double help;
  if (!reader->ReadFromFile(help,"events")) m_events=0;
  else m_events=(long unsigned int)help;
  delete reader;
  if (m_events<m_maxevents) return false;
  for (;hit!=p_gridhandlers->end();++hit) {
    if (!hit->second->ReadIn(ATOOLS::Type::TFStream,OutputPath()+
			     hit->first+m_xsextension)) return false;
  }
  for (Grid_Handler_Map::iterator mit=p_maxhandlers->begin();
       mit!=p_maxhandlers->end();++mit) {
    if (!mit->second->ReadIn(ATOOLS::Type::TFStream,OutputPath()+
			     mit->first+m_xsextension)) return false;
  }
  return true;
}

bool Grid_Creator::InitializeCalculation()
{
  m_criterion=ATOOLS::Variable::TypeToSelectorID(p_xaxis->Variable().Type());
  m_momdata=p_processes->SelectorData()->RemoveData(m_criterion);
  p_processes->SelectorData()->
    SetData(m_criterion,m_momdata.flavs,m_momdata.help,m_gridxmin,m_gridxmax);
  p_processes->ResetSelector(p_processes->SelectorData());
  p_processes->Reset();
  p_processes->CalculateTotalXSec(OutputPath()+OutputFile()+MCExtension());
  return true;
}

bool Grid_Creator::UpdateHistogram(EXTRAXS::XS_Base *const process,
				   const size_t trials)
{
  if (process->Size()==0) return false;
  if ((*process)[0]==process) {
    Amisic_Histogram_Type *histo=m_histograms[process->Name()];
    double weight=process->Differential(p_momenta)*
      process->PSHandler(false)->PSWeight();
    if (process->ISR()->On()==3) {
      weight+=process->Differential2()*process->PSHandler(false)->PSWeight();
    }
    const ATOOLS::Vec4D *p=process->Momenta();
    double value=p[0].PPerp();
    for (size_t i=1;i<4;++i) {
      value=ATOOLS::Max(value,p[i].PPerp());
    }
    histo->Add(value,weight);
    for (size_t i=1;i<trials;++i) histo->Add(value,0.0); 
    return true;
  }
  for (size_t i=0;i<process->Size();++i) {
    if (!UpdateHistogram((*process)[i],trials)) return false;
  }
  return true;
}

bool Grid_Creator::CreateOptimizedGrid()
{
  for (;m_events<m_maxevents;++m_events) {
    ATOOLS::Blob_Data_Base *xsdata=p_processes->WeightedEvent();
    PHASIC::Weight_Info info=xsdata->Get<PHASIC::Weight_Info>();
    p_momenta=p_processes->Selected()->Momenta();
    if (!UpdateHistogram(p_processes,info.ntrial)) return false;
    if (m_events%1000==0) 
      msg_Tracking()<<"Grid_Creator::CreateOptimizedGrid(): "
		    <<"Event "<<m_events<<std::endl;
    if ((m_events%(m_maxevents/100))==0) {
      ATOOLS::msg.Out()<<"\r   "<<(100*m_events)/m_maxevents<<" % ( "
		       <<ATOOLS::rpa.gen.Timer().UserTime()
		       <<" s )   "<<std::flush;
    }
  }
  return true;
}

bool Grid_Creator::CreateInitialGrid()
{
  ATOOLS::rpa.gen.Timer().Start();
  for (m_events=0;m_events<m_initevents;++m_events) {
    ATOOLS::Blob_Data_Base *xsdata=p_processes->WeightedEvent();
    PHASIC::Weight_Info info=xsdata->Get<PHASIC::Weight_Info>();
    p_momenta=p_processes->Selected()->Momenta();
    if (!UpdateHistogram(p_processes,info.ntrial)) return false;
    if (m_events%1000==0) 
      msg_Tracking()<<"Grid_Creator::CreateInitialGrid(): "
		    <<"Event "<<m_events<<std::endl;
    if ((m_events%(m_maxevents/100))==0) 
      ATOOLS::msg.Out()<<"\r   "<<(100*m_events)/m_maxevents<<" % ( "
		       <<ATOOLS::rpa.gen.Timer().UserTime()
		       <<" s )   "<<std::flush;
  }
  return true;
}

bool Grid_Creator::ExportHistogram(const std::string &name) const
{
  Amisic_Histogram_Type *const histo=m_histograms.find(name)->second;
  std::vector<double> xdata(histo->NBins()-2);
  std::vector<double> ydata(histo->NBins()-2);
  std::vector<double> ymax(histo->NBins()-2);
  for (size_t i=1;i<histo->NBins()-1;++i) {
    xdata[i-1]=histo->BinXMean(i);
    ydata[i-1]=histo->BinContent(i);
    ymax[i-1]=histo->BinMax(i);
  }
  Grid_Function_Type *const grid=p_gridhandlers->find(name)->second->Grid();
  grid->SetMonotony(grid->None);
  if (!grid->Import(&xdata,&ydata)) return false;
  if (m_storemax) {
    Grid_Function_Type *const max=p_maxhandlers->find(name)->second->Grid();
    max->SetMonotony(max->None);
    if (!max->Import(&xdata,&ymax)) return false;
  }
  return true;
}

bool Grid_Creator::TranslateGrid() 
{
  bool success=true;
  for (Amisic_Histogram_Map::iterator hit=m_histograms.begin();
       hit!=m_histograms.end();++hit) {
    hit->second->Finish();
    if (!ExportHistogram(hit->first)) success=false;
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
  if (!TranslateGrid()) {
    ATOOLS::msg.Out()<<"Grid_Creator_Base::CreateGrid(..): "
		     <<"Sorry, grid cannot be translated."<<std::endl;
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

bool Grid_Creator::WriteOutGrid(std::vector<std::string> addcomments,
				std::string tempopath)
{
  if (tempopath!="") SetOutputPath(tempopath);
  bool success=true;
  addcomments.push_back("--------------------");
  addcomments.push_back("--------------------");
  for (Grid_Handler_Map::iterator git=p_gridhandlers->begin();
       git!=p_gridhandlers->end();++git) {
    SetOutputFile(git->first+m_xsextension);
    if (m_storemax) {
      addcomments.back()=std::string("max file : ")+
	git->first+m_maxextension;
    }
    if (!WriteSingleGrid(git->second,addcomments)) success=false;
    if (m_storemax) {
      SetOutputFile(git->first+m_maxextension);
      addcomments.back()=std::string("xs file : ")+
	git->first+m_xsextension;
      success=success&&WriteSingleGrid((*p_maxhandlers)[git->first],
				       addcomments);
    }
  }
  SetOutputFile("");
  return success;
}

bool Grid_Creator::WriteSingleGrid(Grid_Handler_Type *grid,
				   std::vector<std::string> addcomments)
{
  if (!CheckOutputFile()) return false;
  std::vector<std::string> comments;
  std::string xvar, yvar;
  xvar=grid->Grid()->XAxis()->Variable().Name();
  yvar=grid->Grid()->YAxis()->Variable().Name();
  comments.push_back(std::string("x : ")+xvar);
  comments.push_back(std::string("y : ")+yvar);
  comments.push_back("--------------------");
  comments.push_back(std::string("x scale : ")+
		     grid->Grid()->XAxis()->Scaling()->Name());
  comments.push_back(std::string("y scale : ")+
		     grid->Grid()->YAxis()->Scaling()->Name());
  comments.push_back("--------------------");
  comments.push_back("boundary conditions ");
  comments.push_back("--------------------");
  comments.push_back(std::string("x_{min} = ")+ATOOLS::ToString(m_gridxmin));
  comments.push_back(std::string("x_{max} = ")+ATOOLS::ToString(m_gridxmax));
  comments.push_back("--------------------");
  comments.push_back(std::string("\\Delta x_{max} = ")+
		     ATOOLS::ToString(m_griddeltax));
  comments.push_back("--------------------");
  comments.push_back(std::string("events = ")+
		     ATOOLS::ToString(m_events));
  if (addcomments.size()!=0) {
    comments.push_back("--------------------");
    for (unsigned int i=0;i<addcomments.size();++i) 
      comments.push_back(addcomments[i]);
  }
  comments.push_back("--------------------");
  comments.push_back("  Data Set follows  ");
  comments.push_back("--------------------");
  msg_Debugging()<<"Grid_Creator::WriteSingleGrid(..): Writing grid to '"
		 <<OutputFile()<<"'"<<std::endl;
  return grid->WriteOut(ATOOLS::Type::TFStream,
			OutputPath()+OutputFile(),comments);
}
