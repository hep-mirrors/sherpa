#include "SHERPA/Single_Events/Output_Phase.H"

#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Output_Phase::Output_Phase(Output_Vector *const outputs,Event_Handler *const h):
  Event_Phase_Handler(""), 
  p_outputs(outputs), m_wit(std::numeric_limits<size_t>::max())
{
  p_eventhandler=h;
  m_type=eph::Analysis;
  for (Output_Vector::iterator it=p_outputs->begin(); it!=p_outputs->end(); ++it) {
    (*it)->SetEventHandler(p_eventhandler);
    (*it)->Header();
    m_name+=(*it)->Name()+"+";
  }
  if (m_name.length()>0) m_name.erase(m_name.length()-1);
  Data_Reader read(" ",";","!","=");
  read.SetAllowUnits(true);
  double wit;
  if (read.ReadFromFile(wit,"FILE_SIZE")) {
    if (wit<1.0) {
      if (wit*rpa->gen.NumberOfEvents()>1.0)
        m_wit=(size_t)(wit*rpa->gen.NumberOfEvents());
    }
    else m_wit=(size_t)(wit);
    msg_Info()<<METHOD<<"(): Set output interval "<<m_wit<<" events.\n";
  }
}

Return_Value::code Output_Phase::Treat(Blob_List *bloblist,double &weight) 
{
  if (!bloblist->empty())
    for (Output_Vector::iterator it=p_outputs->begin(); it!=p_outputs->end(); ++it) {
      (*it)->SetXS(p_eventhandler->TotalXS(),
		   p_eventhandler->TotalErr());
      (*it)->Output(bloblist,weight);
    }
  if (rpa->gen.NumberOfGeneratedEvents()>0 &&
      rpa->gen.NumberOfGeneratedEvents()%m_wit==0 &&
      rpa->gen.NumberOfGeneratedEvents()<rpa->gen.NumberOfEvents()) 
    for (Output_Vector::iterator it=p_outputs->begin(); 
	 it!=p_outputs->end(); ++it)
      (*it)->ChangeFile();
  return Return_Value::Nothing;
}

void Output_Phase::CleanUp(const size_t & mode) 
{
}

void Output_Phase::Finish(const std::string &path)
{
  for (Output_Vector::iterator it=p_outputs->begin(); 
       it!=p_outputs->end(); ++it)
    (*it)->Footer();
}
