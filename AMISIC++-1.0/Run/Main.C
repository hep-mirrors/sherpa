#include "Amisic.H"

int main(int argc, const char *argv[]) 
{
  std::string arguments;
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  for (int i=0;i<argc;++i) arguments+=std::string(argv[i])+std::string(" ");
  reader->SetString(arguments);
  std::string inputpath, inputfile;
  int testevents;
  if (!reader->ReadFromString(inputpath,"-p")) { 
    ATOOLS::msg.Out()<<"main: No input path specified. Using default value './'"<<std::endl;
    inputpath=std::string("./");
  }
  if (!reader->ReadFromString(inputfile,"-f")) { 
    ATOOLS::msg.Out()<<"main: No input file specified. Using default value 'MI.dat'"<<std::endl;
    inputfile=std::string("MI.dat");
  }
  if (!reader->ReadFromString(testevents,"-e")) { 
    testevents=0;
  }
  std::vector<std::string> options;
  bool creategrid=false, debug=false;
  if (reader->VectorFromString(options,ATOOLS::nullstring,arguments,reader->VHorizontal)) { 
    for (std::vector<std::string>::iterator it=options.begin();it!=options.end();++it) {
      if (*it==std::string("-c")) creategrid=true;
      if (*it==std::string("-d")) debug=true;
    }
  }
  delete reader;

  int result=0;

  AMISIC::Amisic *test = new AMISIC::Amisic();
  result=1*(1-test->Initialize(inputpath,inputfile,creategrid));
  if (result!=1) {
    for (int i=0;i<testevents;++i) {
      ATOOLS::Blob_List *blobs=test->CreateProcesses();
      if (debug) { 
	std::cout<<"***** new try *****"<<std::endl;
	for (ATOOLS::Blob_Iterator bit=blobs->begin();bit!=blobs->end();std::cout<<*bit++<<std::endl);
      }
      delete blobs;
    }
  }
  delete test;
  return result;
}
