#ifdef ROOT_SUPPORT
#include "My_Root.H"

using namespace MYROOT;

My_Root::My_Root(const int argc,char **const argv):
  p_file(NULL),
  m_draw(false)
{
  std::string path, file;
#ifndef USING__My_Root_only
  std::string inputstring;
  for (int i=0;i<argc;++i) inputstring+=std::string(argv[i])+std::string("; ");
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  reader->SetString(inputstring);
  if (!reader->ReadFromString(path,"ROOT_PATH")) path="./";
  if (!reader->ReadFromString(file,"ROOT_FILE")) file="output.root";
  delete reader;
#else
  path="./";
  file="output.root";
#endif      
  SetOutputPath(path);
  SetOutputFile(file);
  int argcf=1;
  char **argvf = new char*[1];
  p_root = new TApplication("MyRoot",&argcf,argvf);
  if ((OutputPath()+OutputFile())!="") {
    if (!system((std::string("test -f ")+OutputPath()+OutputFile()).c_str())) {
      system((std::string("rm -f ")+OutputPath()+OutputFile()).c_str());
    }
    p_file = new TFile((OutputPath()+OutputFile()).c_str(),"recreate");
  }
} 

My_Root::~My_Root() 
{
  Draw();
  p_root->Run(kTRUE);
  if (p_file!=NULL) {
    p_file->Write();
    delete p_file;
  }
  delete p_root;
}

void My_Root::Draw() 
{
  if (!m_draw) return;
  for (String_Object_Map::const_iterator oit=m_objects.begin();
       oit!=m_objects.end();++oit) {
    new TCanvas(oit->first.c_str(),oit->first.c_str());
    oit->second->Draw();
  }
}

#ifndef USING__My_Root_only
void My_Root::PrepareTerminate()
{
  Draw();
  p_root->Run(kTRUE);
  if (p_file!=NULL) p_file->Write();
}
#endif

bool My_Root::AddObject(TObject *const object,const std::string &key) 
{ 
  if (m_objects.find(key)==m_objects.end()) {
    m_objects.insert(String_Object_Map::value_type(key,object)); 
    return true;
  }
  return false;
}

#endif
