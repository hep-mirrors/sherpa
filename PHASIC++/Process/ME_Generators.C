#include "PHASIC++/Process/ME_Generators.H"

#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace ATOOLS;
using namespace PHASIC;

ME_Generators::ME_Generators(const std::string &path,
			     const std::string &file):
  m_path(path), m_file(file)
{
  Default_Reader reader;
  reader.SetInputPath(m_path);
  reader.SetInputFile(m_file);
  std::vector<std::string> megens;
  if (!reader.ReadStringVectorNormalisingNoneLikeValues(megens,"ME_SIGNAL_GENERATOR")) {
    megens.push_back("Comix");
    megens.push_back("Amegic");
    megens.push_back("Internal");
  }
  for (size_t i(0);i<megens.size();++i) {
    if (megens[i]=="None") continue;
    push_back(ME_Generator_Getter::GetObject(megens[i],ME_Generator_Key()));
    if (back()==NULL) {
      msg_Info()<<METHOD<<"(): Try loading '"<<megens[i]
		<<"' from 'libSherpa"<<megens[i]<<"'."<<std::endl;
      if (s_loader->LoadLibrary("Sherpa"+megens[i]))
	back()=ME_Generator_Getter::GetObject(megens[i],ME_Generator_Key());
    }
    if (back()==NULL) {
      msg_Error()<<METHOD<<"(): ME generator '"<<megens[i]
                 <<"' not found. Ignoring it."<<std::endl;
      pop_back();
    }
  }
  for (size_t i(0);i<size();++i) {
    rpa->gen.SetVariable(at(i)->Name(),ToString(at(i)));
  }
}

ME_Generators::~ME_Generators()
{
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    delete *mit;
  }
}

bool ME_Generators::LoadGenerator(const std::string &name)
{
  for (size_t i(0);i<size();++i)
    if (at(i)->Name()==name) return true;
  push_back(ME_Generator_Getter::GetObject(name,ME_Generator_Key()));
  if (back()==NULL) {
    msg_Info()<<METHOD<<"(): Try loading '"<<name
	      <<"' from 'libSherpa"<<name<<"'."<<std::endl;
    if (s_loader->LoadLibrary("Sherpa"+name))
      back()=ME_Generator_Getter::GetObject(name,ME_Generator_Key());
  }
  if (back()==NULL) {
    msg_Error()<<METHOD<<"(): ME generator '"<<name
	       <<"' not found. Ignoring it."<<std::endl;
    pop_back();
    return false;
  }
  if (!back()->Initialize(m_path,m_file,p_model,p_beam,p_isr)) return false;
  back()->SetGenerators(this);
  return true;
}

bool ME_Generators::InitializeGenerators(MODEL::Model_Base *model,
                                         BEAM::Beam_Spectra_Handler *beam,
                                         PDF::ISR_Handler *isr)
{
  p_isr=isr;
  p_beam=beam;
  p_model=model;
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if (!(*mit)->Initialize(m_path,m_file,model,beam,isr)) return false;
    (*mit)->SetGenerators(this);
  }
  return true;
}

int ME_Generators::PerformTests()
{ 
  int result(1);
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    int ret((*mit)->PerformTests());
    if (ret==0) return 0;
    else if (ret==-1)
      result = -1;
  }
  return result;
}

bool ME_Generators::NewLibraries()
{
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if ((*mit)->NewLibraries()) return true;
  }
  return false;
}

Process_Base* ME_Generators::InitializeProcess(const Process_Info &pi, bool add)
{
  DEBUG_FUNC(&pi);
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if (pi.m_megenerator!="" && (*mit)->Name()!=pi.m_megenerator) continue;
    DEBUG_INFO("trying "<<(*mit)->Name());
    Process_Base *proc((*mit)->InitializeProcess(pi,add));
    if (proc) {
      DEBUG_INFO("found "<<proc->Name());
      return proc;
    }
  }
  DEBUG_INFO("couldn't initialize process.");
  return NULL;
}

