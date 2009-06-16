#include "PHASIC++/Process/ME_Generators.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Process/ME_Generator_Base.H"

using namespace ATOOLS;
using namespace PHASIC;

ME_Generators::ME_Generators(const std::string &path,
                                       const std::string &file) :
  std::vector<ME_Generator_Base*>(), m_path(path), m_file(file)
{
  Data_Reader read(" ",";","!","=");
  read.SetInputPath(m_path);
  read.SetInputFile(m_file);
  std::vector<std::string> megens;
  if (!read.VectorFromFile(megens,"ME_SIGNAL_GENERATOR")) {
    megens.push_back("Internal");
    megens.push_back("Amegic");
    megens.push_back("Comix");
  }
  for (size_t i(0);i<megens.size();++i) {
    push_back(ME_Generator_Getter::GetObject(megens[i],ME_Generator_Key()));
    if (back()==NULL) {
      msg_Error()<<METHOD<<"(): ME generator '"<<megens[i]
                 <<"' not found. Ignoring it."<<std::endl;
      pop_back();
    }
  }
  for (size_t i(0);i<size();++i) {
    rpa.gen.SetVariable(at(i)->Name(),ToString
                        ((long unsigned int)at(i)));
  }
}

ME_Generators::~ME_Generators()
{
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    delete *mit;
  }
}

bool ME_Generators::InitializeGenerators(MODEL::Model_Base *model,
                                         BEAM::Beam_Spectra_Handler *beam,
                                         PDF::ISR_Handler *isr)
{
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if (!(*mit)->Initialize(m_path,m_file,model,beam,isr)) return false;
  }
  return true;
}

bool ME_Generators::PerformTests()
{
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if (!(*mit)->PerformTests()) return false;
  }
  return true;
}

Process_Base* ME_Generators::InitializeProcess(const Process_Info &pi, bool add)
{
  DEBUG_FUNC(&pi);
  size_t nf(pi.m_fi.NExternal());
  for (ME_Generators::const_iterator mit=begin(); mit!=end(); ++mit) {
    if (pi.m_fi.NLOType()==nlo_type::loop && pi.m_loopgenerator!=(*mit)->Name())
      continue;
    DEBUG_INFO("trying "<<(*mit)->Name());
    if (nf>4 && mit!=--end()) continue;
    Process_Base *proc((*mit)->InitializeProcess(pi,add));
    if (proc) {
      DEBUG_INFO("found "<<proc->Name());
      return proc;
    }
  }
  DEBUG_INFO("couldn't initialize process.");
  return NULL;
}

