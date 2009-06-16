#include "ATOOLS/Org/My_Root.H"

#ifdef USING__ROOT
#include "ATOOLS/Org/Shell_Tools.H"
#include "TStyle.h"

#include <sys/stat.h>

using namespace MYROOT;

My_Root *MYROOT::myroot=NULL;

My_Root::My_Root(const int argc,char **const argv):
  p_file(NULL),
  m_drawmode(0)
{
  std::string path, file;
#ifndef USING__My_Root_only
  std::string inputstring;
  for (int i=0;i<argc;++i) inputstring+=std::string(argv[i])+std::string("; ");
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  reader->SetString(inputstring);
  if (!reader->ReadFromString(path,"ROOT_PATH")) path="./Analysis/";
  if (!reader->ReadFromString(file,"ROOT_FILE")) file="output.root";
  if (!reader->ReadFromString(m_drawmode,"DRAW_ROOT_RESULTS")) m_drawmode=0;
  delete reader;
#else
  path="./";
  file="output.root";
#endif      
  SetOutputPath(path);
  SetOutputFile(file);
  int argcf=1;
  char **argvf = new char*[1];
  argvf[0] = "";
  p_root = new TApplication("MyRoot",&argcf,argvf);
  if ((OutputPath()+OutputFile())!="") {
    ATOOLS::MakeDir(OutputPath());
    struct stat fst;
    if (stat((OutputPath()+OutputFile()).c_str(),&fst)!=-1 && 
	(fst.st_mode&S_IFMT)==S_IFREG) {
      remove((OutputPath()+OutputFile()).c_str());
    }
    p_file = new TFile((OutputPath()+OutputFile()).c_str(),"recreate");
  }
  delete [] argvf;
} 

My_Root::~My_Root() 
{
  Draw();
  if (p_file!=NULL) {
    p_file->Write();
    delete p_file;
  }
  delete p_root;
}

void My_Root::Draw() 
{
  if (m_drawmode==0) return;
  gStyle->SetPalette(1);
  for (String_Object_Map::const_iterator oit=m_objects.begin();
       oit!=m_objects.end();++oit) {
    new TCanvas((oit->first+"_c").c_str(),(oit->first+"_c").c_str());
    if (m_drawmode&2) gPad->SetLogx();
    if (m_drawmode&4) gPad->SetLogy();
    if (m_drawmode&8) gPad->SetLogz();
    oit->second->Draw(m_drawoption.c_str());
  }
  if (m_objects.size()>0) p_root->Run(kTRUE);
}

#ifndef USING__My_Root_only
void My_Root::PrepareTerminate()
{
  Draw();
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
