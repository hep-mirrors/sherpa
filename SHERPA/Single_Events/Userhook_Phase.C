#include "SHERPA/Single_Events/Userhook_Phase.H"

#include "SHERPA/Single_Events/Event_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include <limits>

using namespace SHERPA;
using namespace ATOOLS;

Userhook_Phase::Userhook_Phase(Sherpa* sherpa):
  Event_Phase_Handler("")
{
  m_type=eph::Userhook;
  InitializeHooks(sherpa);

  if (m_userhooks.empty()) {
    m_name="None";
  }
  else {
    for (Userhook_Vector::iterator it=m_userhooks.begin(); it!=m_userhooks.end(); ++it) {
      m_name+=(*it)->Name()+"+";
    }
    if (m_name.length()>0) m_name.pop_back();
  }
}

Userhook_Phase::~Userhook_Phase()
{
  while (m_userhooks.size()>0) {
    delete m_userhooks.back();
    m_userhooks.pop_back();
  }
}

void Userhook_Phase::InitializeHooks(Sherpa* sherpa)
{
  const auto userhooks = Settings::GetMainSettings()["USERHOOKS"]
    .SetDefault("None")
    .UseNoneReplacements()
    .GetVector<std::string>();
  for (size_t i=0; i<userhooks.size(); ++i) {
    if (userhooks[i]=="None") continue;
    Userhook_Base* userhook=Userhook_Base::Getter_Function::GetObject
      (userhooks[i],Userhook_Arguments(sherpa));
    if (userhook==NULL) {
      if (!s_loader->LoadLibrary("SherpaUserhook"+userhooks[i]))
	THROW(missing_module,"Cannot load userhook library SherpaUserhook"+userhooks[i]);
      userhook=Userhook_Base::Getter_Function::GetObject
        (userhooks[i],Userhook_Arguments(sherpa));
    }
    if (userhook==NULL) THROW(fatal_error,"Cannot initialize "+userhooks[i]+" userhook");
    m_userhooks.push_back(userhook);
  }

}

Return_Value::code Userhook_Phase::Treat(Blob_List* bloblist)
{
  int success(0), newevent(0);
  for (Userhook_Vector::iterator it=m_userhooks.begin(); it!=m_userhooks.end(); ++it) {
    switch ((*it)->Run(bloblist)) {
    case Return_Value::Nothing :
      break;
    case Return_Value::Success :
      success++;
      break;
    case Return_Value::New_Event :
      newevent++;
      break;
    case Return_Value::Error :
      return Return_Value::Error;
    default:
      THROW(fatal_error,"Invalid return value");
    }
  }
  
  if (newevent) return Return_Value::New_Event;
  else if (success) return Return_Value::Success;
  else return Return_Value::Nothing;
}

void Userhook_Phase::Finish(const std::string &path)
{
  for (Userhook_Vector::iterator it=m_userhooks.begin(); it!=m_userhooks.end(); ++it) {
    (*it)->Finish();
  }
}

void Userhook_Phase::CleanUp(const size_t & mode) 
{
  for (Userhook_Vector::iterator it=m_userhooks.begin(); it!=m_userhooks.end(); ++it) {
    (*it)->CleanUp();
  }
}
