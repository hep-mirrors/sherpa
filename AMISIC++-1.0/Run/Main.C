#include "Amisic.H"
#include "Particle.H"
#include "TApplication.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

int main(int argc, const char *argv[]) 
{
  std::string arguments;
  ATOOLS::msg.SetLevel(15);
  ATOOLS::Data_Reader *reader = new ATOOLS::Data_Reader();
  for (int i=0;i<argc;++i) arguments+=std::string(argv[i])+std::string(" ");
  reader->SetString(arguments);
  std::string inputpath, inputfile, outputfile;
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
  bool debug=false, analyse=false, display=false;
  if (reader->VectorFromString(options,ATOOLS::nullstring,arguments,reader->VHorizontal)) { 
    for (std::vector<std::string>::iterator it=options.begin();it!=options.end();++it) {
      if (*it==std::string("-d")) debug=true;
      if (*it==std::string("-a")) analyse=true;
      if (*it==std::string("-D")) display=true;
    }
  }
  TApplication *myroot=NULL;
  TFile *myfile=NULL;
  TH1D *multiplicity=NULL, *check=NULL;
  if (analyse) {
    if (!reader->ReadFromString(outputfile,"-o")) { 
      ATOOLS::msg.Out()<<"main: No output file specified. Using default value 'output.root'"<<std::endl;
      outputfile=std::string("output.root");
    }
    int argcf=1;
    char **argvf = new char*[1];
    myroot = new TApplication("MyRoot",&argcf,argvf);
    delete [] argvf;
    myfile = new TFile(outputfile.c_str(),"recreate",outputfile.c_str());
    multiplicity = new TH1D("multiplicity","multiplicity",120,0.0,120.0);
  }
  delete reader;
  
  int result=0;

  AMISIC::Amisic *test = new AMISIC::Amisic();
  test->SetInputPath(inputpath);
  test->SetInputFile(inputfile);
  result=1-test->Initialize();
  double mfw[6], msqrfw[6], mfwbw[6];
  for (unsigned int i=0;i<6;++i) mfw[i]=msqrfw[i]=mfwbw[i]=0.0;
  if (result!=1) {
    ATOOLS::Blob_List *blobs = new ATOOLS::Blob_List();
    for (int i=0;i<testevents;++i) {
      if (i%(testevents/10)==0) std::cout<<" Event "<<i<<std::endl;
      test->CleanUp();
      test->GenerateEvent(blobs);
      if (analyse) {
 	multiplicity->Fill((double)blobs->size(),1.0);
      }
      while (blobs->size()>0) {
	delete *blobs->begin();
	blobs->erase(blobs->begin());
      }
    }
    delete blobs;
  }
  delete test;
  if (analyse) {
    myfile->Write();
    if (display) {
      TCanvas *mcanvas = new TCanvas("mresults","results");
      mcanvas->cd();
      multiplicity->Draw();
      myroot->Run(kTRUE);
    }
    delete myfile;
    delete myroot;
  }
  return result;
}
