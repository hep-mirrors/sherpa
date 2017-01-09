#include "AMISIC++/Model/Grid_Creator.H"

#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/Weight_Info.H"

using namespace AMISIC;
using namespace ATOOLS;

Grid_Creator::Grid_Creator(Amisic_Histogram_Map *histograms,
                           const std::vector<EXTRAXS::Process_Group*>& processes):
  p_histograms(histograms),
  p_processes(processes),
  m_xsextension("_xs.dat"),
  m_datatag("[x,w,w2,max,n]"),
  m_events(0)
{
  if (p_processes.empty()) {
    THROW(fatal_error,"Process handler is not initialized");
  }
  for (size_t i=0; i<p_processes.size(); ++i) {
    if (!CollectProcesses(p_processes[i])) {
      THROW(fatal_error,"Process handler does not own any process");
    }
  }
}

Grid_Creator::~Grid_Creator()
{
}

bool Grid_Creator::CollectProcesses(PHASIC::Process_Base *const process)
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
  if (InputFile()=="") return false;
  ATOOLS::Default_Reader reader;
  reader.SetInputPath(InputPath());
  reader.SetInputFile(InputFile());
  std::vector<std::string> helps;
  if (!reader.ReadVector(helps,"X_VARIABLE")) m_gridxvariable="PT";
  else m_gridxvariable=MakeString(helps);
  if (!reader.ReadVector(helps,"Y_VARIABLE")) m_gridyvariable="";
  else m_gridyvariable=MakeString(helps);
  if (m_gridxmin==0.0) {
    m_gridxmin=std::numeric_limits<double>::max();
    for (size_t i=0; i<p_processes.size(); ++i) {
      m_gridxmin=std::min(m_gridxmin, sqrt
        (ATOOLS::Max(p_processes[i]->Integrator()->ISR()->PDF(0)->Q2Min(),
		   p_processes[i]->Integrator()->ISR()->PDF(1)->Q2Min())));
    }
  }

  m_gridxmin=ATOOLS::Max(m_gridxmin,1.e-3);
  m_griddeltax = reader.Get("GRID_DELTA_X", (log(m_gridxmax)-log(m_gridxmin))/250.);
  m_initevents = reader.Get("INITIAL_EVENTS", 0);
  m_maxevents = reader.Get("MAX_EVENTS", 100000);
  m_binerror = reader.Get("GRID_ERROR", 0.05);
  m_outputlevel = reader.Get("GRID_CREATOR_OUTPUT", 2);
  m_gridxscaling = reader.Get<std::string>("HISTO_X_SCALING", "Log_B_10");
  m_gridyscaling = reader.Get<std::string>("HISTO_Y_SCALING", "Id");

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
  My_In_File::OpenDB(OutputPath());
  double sum=0.0;
  msg_Info()<<METHOD<<"(): Reading grid ";
  msg_Tracking()<<"{\n";
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    if (hit->second->ReadIn(OutputPath()+hit->first+m_xsextension,
			    m_datatag)) {
      if (hit->second->XMin()-m_gridxmin>m_gridxmin*1.0e-7 ||
	  m_gridxmax-hit->second->XMax()>m_gridxmax*1.0e-7 ||
	  hit->second->Entries()<m_initevents) {
	Clear();
	My_In_File::CloseDB(OutputPath());
	return false;
      }
      double cur=hit->second->Norm()*rpa->Picobarn();
      msg_Info()<<'.'<<std::flush;
      msg_Tracking()<<"  '"<<hit->first<<"' -> "<<cur<<" pb\n";
      sum+=cur;
    }
    else {
      Clear();
      My_In_File::CloseDB(OutputPath());
      return false;
    }
  }
  msg_Info()<<" done."<<std::endl;
  msg_Tracking()<<"} -> sum = "<<sum<<" pb\n";
  My_In_File::CloseDB(OutputPath());
  return true;
}

bool Grid_Creator::InitializeCalculation()
{
  m_criterion=p_xaxis->Variable()->SelectorID();
  std::vector<std::string> bounds(3);
  bounds[0]=ToString(kf_jet);
  bounds[1]=ToString(m_gridxmin);
  bounds[2]=ToString(m_gridxmax);
  Data_Reader read;
  for (size_t i=0; i<p_processes.size(); ++i) {
    PHASIC::Selector_Key skey(p_processes[i]->Integrator(),&read);
    skey.SetData(m_criterion,bounds);
    p_processes[i]->SetSelector(skey);
    p_processes[i]->Integrator()->PSHandler()->InitCuts();
  }
  return true;
}

bool Grid_Creator::UpdateHistogram(PHASIC::Process_Base *const process)
{
  if (process->IsGroup()) {
    for (size_t i=0;i<process->Size();++i)
      if (!UpdateHistogram((*process)[i])) return false;
    return true;
  }
  Amisic_Histogram_Type *histo=(*p_histograms)[process->Name()];
  ATOOLS::Weight_Info *info=process->OneEvent(0);
  const ATOOLS::Vec4D_Vector &p=process->Integrator()->Momenta();
  double value=(*p_variable)(&p[0]);
  for (size_t i=1;i<4;++i) value=ATOOLS::Max(value,(*p_variable)(&p[i]));
  if (info==NULL) {
    histo->Add(value,0.0); 
    return true;
  }
  histo->Add(value,info->m_weight);
  delete info;
  return true;
}

bool Grid_Creator::CreateGridInternal()
{
  msg_Info()<<METHOD<<": Initializing grid for MI.\n";
  double starttime=ATOOLS::rpa->gen.Timer().UserTime();
  msg_Info()<<ATOOLS::tm::curoff;
  for (;m_events<m_maxevents;++m_events) {
    for (size_t i=0; i<p_processes.size(); ++i) {
      if (!UpdateHistogram(p_processes[i])) return false;
    }
    if ((m_events%(m_maxevents/100))==0 && m_events>0) {
      double diff=ATOOLS::rpa->gen.Timer().UserTime()-starttime;
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

bool Grid_Creator::WriteOutGrid(std::vector<std::string> addcomments,
				const std::string &path) 
{
#ifdef USING__MPI
  if (MPI::COMM_WORLD.Get_rank()) return true;
#endif
  bool success=true;
  for (Amisic_Histogram_Map::iterator hit=p_histograms->begin();
       hit!=p_histograms->end();++hit) {
    if (!hit->second->WriteOut((path!=""?path:OutputPath())
			       +hit->first+m_xsextension,m_datatag,
			       addcomments)) success=false;
  }
  return success;
}

bool Grid_Creator::CreateGrid()
{
  bool success=true;
  int formerlevel=msg->Level();
  My_Out_File::OpenDB(OutputPath());
  msg_Info()<<"Grid_Creator::CreateGrid(): "
	    <<"Calculating grid {"<<std::endl;
  msg->SetLevel(m_outputlevel);
  for (size_t i=0; i<p_processes.size(); ++i) {
    p_processes[i]->CalculateTotalXSec(OutputPath()+OutputFile(),true);
  }
  if (!CreateGridInternal()) {
    msg_Out()<<"Grid_Creator_Base::CreateGrid(..): "
	     <<"Grid creation failed."<<std::endl;
    success=false;
  }
  My_Out_File::ExecDB(OutputPath(),"begin");
  if (!WriteOutGrid()) {
    msg_Out()<<"Grid_Creator_Base::CreateGrid(..): "
		     <<"Sorry, grid cannot be written to '"
		     <<OutputFile()<<"'"<<std::endl;
    success=false;
  }
  My_Out_File::ExecDB(OutputPath(),"commit");
  My_Out_File::CloseDB(OutputPath());
  msg->SetLevel(formerlevel);
  msg_Info()<<"\n}"<<std::endl;
  return success;
}

