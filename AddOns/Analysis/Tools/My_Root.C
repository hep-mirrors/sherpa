#include "AddOns/Analysis/Tools/My_Root.H"

#ifdef USING__ROOT
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"
#include "TStyle.h"

#include <sys/stat.h>

using namespace MYROOT;

My_Root *MYROOT::myroot=NULL;

My_Root::My_Root():
  p_file(NULL)
{
  std::string path, file;
  std::string inputstring;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  reader->SetString(inputstring);
  if (!reader->ReadFromString(path,"ROOT_PATH")) path="./Analysis/";
  if (!reader->ReadFromString(file,"ROOT_FILE")) file="output.root";
  delete reader;
  SetOutputPath(path);
  SetOutputFile(file);
  int argcf=1;
  char **argvf = new char*[1];
  argvf[0] = "";
  p_root = new TApplication("MyRoot",&argcf,argvf);
  delete [] argvf;
} 

My_Root::~My_Root() 
{
  if (p_file!=NULL) {
    p_file->Write();
    delete p_file;
  }
  delete p_root;
}

TObject *My_Root::GetObject(const std::string &key)
{
  return m_objects.find(key)->second;
}

std::string My_Root::GetKey()
{
  return "Key_"+ATOOLS::ToString(m_objects.size());
}

void My_Root::PrepareTerminate()
{
  if (p_file!=NULL) p_file->Write();
}

void My_Root::InitFile()
{
  if (p_file) return;
  if (OutputPath()=="" && OutputFile()=="") return;
  ATOOLS::MakeDir(OutputPath());
  struct stat fst;
  if (stat((OutputPath()+OutputFile()).c_str(),&fst)!=-1 && 
      (fst.st_mode&S_IFMT)==S_IFREG) {
    remove((OutputPath()+OutputFile()).c_str());
  }
  p_file = new TFile((OutputPath()+OutputFile()).c_str(),"recreate");
}

bool My_Root::AddObject(TObject *const object,const std::string &key) 
{ 
  if (m_objects.find(key)==m_objects.end()) {
    m_objects.insert(String_Object_Map::value_type(key,object)); 
    return true;
  }
  return false;
}

#endif
