#include "Amisic.H"

int main(int argc, const char *argv[]) 
{
  std::string arguments;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  for (int i=0;i<argc;++i) arguments+=std::string(argv[i])+std::string(" ");
  reader->SetString(arguments);
  std::string inputpath, inputfile;
  if (!reader->ReadFromString(inputpath,"-p")) { 
    ATOOLS::msg.Out()<<"main: No input path specified. Using default value './'"<<std::endl;
    inputpath=std::string("./");
  }
  if (!reader->ReadFromString(inputfile,"-f")) { 
    ATOOLS::msg.Out()<<"main: No input file specified. Using default value 'MI.dat'"<<std::endl;
    inputfile=std::string("MI.dat");
  }
  std::vector<std::string> options;
  bool creategrid=false;
  if (reader->VectorFromString(options,"-",arguments,reader->VHorizontal)) { 
    for (std::vector<std::string>::iterator it=options.begin();it!=options.end();++it) {
      if (*it==std::string("c")) creategrid=true;
    }
  }
  delete reader;

  int result=0;

  AMISIC::Amisic *test = new AMISIC::Amisic();
  result=1*(1-test->Initialize(inputpath,inputfile,creategrid));
  delete test;
  return result;
}
