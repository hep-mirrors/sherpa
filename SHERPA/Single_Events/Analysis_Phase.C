#include "SHERPA/Single_Events/Analysis_Phase.H"

#include "SHERPA/Tools/Analysis_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Analysis_Phase::Analysis_Phase(Analysis_Vector *const analyses):
  Event_Phase_Handler(""), 
  p_analyses(analyses),
  m_wit(std::numeric_limits<size_t>::max())
{
  m_type=eph::Analysis;
  for (Analysis_Vector::iterator it=p_analyses->begin(); it!=p_analyses->end(); ++it) {
    m_name+=(*it)->Name()+"+";
    m_inits[*it]=false;
  }
  if (m_name.length()>0) m_name.pop_back();

  Settings& s = Settings::GetMainSettings();
  const double wit{ s["ANALYSIS_WRITEOUT_INTERVAL"]
      .SetDefault(static_cast<double>(std::numeric_limits<size_t>::max()))
      .Get<double>() };
  if (wit<1.0) {
    if (wit*rpa->gen.NumberOfEvents()>1.0)
      m_wit=(size_t)(wit*rpa->gen.NumberOfEvents());
  } else {
    m_wit=(size_t)(wit);
  }
}

Return_Value::code Analysis_Phase::Treat(Blob_List* bloblist)
{
  if (!bloblist->empty())
    for (Analysis_Vector::iterator it=p_analyses->begin(); it!=p_analyses->end(); ++it) {
      if (!m_inits[*it]) m_inits[*it]=(*it)->Init();
      if (!(*it)->Run(bloblist)) return Return_Value::New_Event;
    }
  if (rpa->gen.NumberOfGeneratedEvents()>0 &&
      rpa->gen.NumberOfGeneratedEvents()%m_wit==0 &&
      rpa->gen.NumberOfGeneratedEvents()<rpa->gen.NumberOfEvents()) 
    for (Analysis_Vector::iterator it=p_analyses->begin(); 
	 it!=p_analyses->end(); ++it)
      (*it)->WriteOut();
  return Return_Value::Nothing;
}

void Analysis_Phase::CleanUp(const size_t & mode) 
{
  for (Analysis_Vector::iterator it=p_analyses->begin(); 
       it!=p_analyses->end(); ++it)
    (*it)->CleanUp();
}

void Analysis_Phase::Finish(const std::string &path)
{
  for (Analysis_Vector::iterator it=p_analyses->begin(); 
       it!=p_analyses->end(); ++it)
    (* it)->Finish();
}
