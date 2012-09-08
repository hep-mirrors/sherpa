#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <termios.h>
#include <unistd.h>
#include <stdlib.h>
#include <vector>

using namespace ATOOLS;
using namespace std;

void PrintInfo()
{
  cout<<"usage: Combine_Analysis [options] OUTPUT_DIRECTORY INPUT_DIRECTORY1 INPUT_DIRECTORY2 ..."<<endl;
  cout<<endl;
  cout<<" options: -m=xxx  specify combination mode:"<<endl;
  cout<<"             cmb  (default) combines histograms as if they were generated in one run"<<endl; 
  cout<<"             opt  combines histogram bins, weighted according to their statistical errors"<<endl;
  cout<<"             add  adds histograms"<<endl;
  cout<<"             min  determines bin by bin minima"<<endl;
  cout<<"             max  determines bin by bin maxima"<<endl;
  cout<<endl;
  cout<<"          -check  check statistical compatibility between data sets"<<endl;
  cout<<endl;
  cout<<"          -filter=<filter rule> only combine histogram files that pass the <filter rule>"<<endl;
  cout<<endl;
  cout<<"          -noheader  remove histogram header in the output;"<<endl;
  cout<<"                     further combination is not possible without the header"<<endl;
  cout<<endl;
  cout<<"          -rescale=<factor> rescale output by <factor>"<<endl;
}

#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#ifdef USING__MPI
#include "mpi.h"
#endif

int main(int argc,char **argv)
{
#ifdef USING__MPI
  MPI::Init(argc,argv);
#endif
  ATOOLS::exh = new Exception_Handler();
  ATOOLS::msg = new Message();
  if (argc<3){
    PrintInfo();
    return 0;
  }
  int i=1;
  int mode=0,noheader=0,check=0;
  double factor=1.;
  string filter("");
  while (argv[i][0]=='-') {
    string argstr=argv[i];
    if (argstr==string("-noheader")) noheader=1;
    else if (argstr==string("-check")) check=1;
    else if (argstr==string("-m=add")) mode=1;
    else if (argstr==string("-m=opt")) mode=2;
    else if (argstr==string("-m=min")) mode=3;
    else if (argstr==string("-m=max")) mode=4;
    else if (argstr.find("-rescale=")!=string::npos) {
      mode=5;
      argstr=argstr.substr(9);
      factor=ToType<double>(argstr);
      cout<<"rescaling by "<<factor<<endl;
    }
    else if (argstr.find("-filter=")!=string::npos) {
      argstr=argstr.substr(8);
      if (argstr[0]=='"') argstr=argstr.substr(1);
      if (argstr[argstr.length()-1]=='"') argstr=argstr.substr(0,argstr.length()-1);
      filter=argstr+" ";
    }
    else if (argstr!=string("-m=cmb")) {
      cout<<"unrecognized option "<<argstr<<endl;
      return 0;
    }
    i++;
    if (mode==1) check=0;
  }
  if (argc-i<2){
    PrintInfo();
    return 0;
  }

  string output=argv[i];
  output+="/";
  i++;
  vector<string> inlist;
  for (;i<argc;++i) {
    inlist.push_back(string(argv[i])+"/");
  }
  MakeDir(output);
  string flname=output+"fl.tmp";
  string tmp="ls "+filter+inlist[0]+" > "+flname;
  system(tmp.c_str());
  vector<string> filelist;
  std::string buf;
  ifstream from(flname.c_str());
  while (from) {
    getline(from,buf);
    if (buf.find(".dat")!=string::npos)
      if (FileExists(inlist[0]+buf)) filelist.push_back(buf);
  }
  from.close();
  system(("rm "+flname).c_str());

  if (check) {
    double** csmatrix=new double*[inlist.size()];
    int** nummatrix=new int*[inlist.size()];
    string** namematrix=new string*[inlist.size()];
    for (size_t j=0;j<inlist.size();j++){
      csmatrix[j]=new double[inlist.size()];
      nummatrix[j]=new int[inlist.size()];
      namematrix[j]=new string[inlist.size()];
      for (size_t k=0;k<inlist.size();k++) {
	csmatrix[j][k]=1.;
	nummatrix[j][k]=0;
      }
    }
    for (i=0;i<filelist.size();i++) {
      vector<Histogram*> hvec;
      for (size_t j=0;j<inlist.size();j++) {
	hvec.push_back(new Histogram (inlist[j]+filelist[i]));
      }
      for (size_t j=0;j<inlist.size();j++) {
	for (size_t k=j+1;k<inlist.size();k++) {
	  double avgdev(0.), stddev(0.);
	  if (hvec[j]->CheckStatistics(*(hvec[k]),avgdev,stddev)>1) {
	    if (stddev>csmatrix[j][k]) {
	      csmatrix[j][k]=stddev;
	      namematrix[j][k]=filelist[i];
	    }
	    if (stddev<csmatrix[k][j]) {
	      csmatrix[k][j]=stddev;
	      namematrix[k][j]=filelist[i];
	    }
	    if (stddev>1.5) nummatrix[j][k]++;
	    if (stddev<0.5) nummatrix[k][j]++;
	  }
	}
      }
      for (size_t j=0;j<hvec.size();j++) delete hvec[j];
      hvec.clear();
    }
    bool report=0;
    cout<<"Statistics check:"<<endl;
    for (size_t j=0;j<inlist.size();j++) {
      for (size_t k=j+1;k<inlist.size();k++) {
	if (csmatrix[j][k]>1.5) {
	  cout<<" Data sets "<<inlist[j]<<" and "<<inlist[k]<<" incompatible in "
	      <<nummatrix[j][k]<<" histograms"<<endl;
	  cout<<"   Maximal deviation="<<csmatrix[j][k]<<" sigma/bin in "<<namematrix[j][k]<<endl; 
	  report=1;
	}
      }
    }
    for (size_t j=0;j<inlist.size();j++) {
      for (size_t k=j+1;k<inlist.size();k++) {
	if (csmatrix[k][j]==0.) {
	  cout<<" Data sets "<<inlist[j]<<" and "<<inlist[k]<<" identical in "
	      <<nummatrix[k][j]<<" histograms"<<endl;
	  report=1;
	}
	else if (csmatrix[k][j]<0.5) {
	  cout<<" Data sets "<<inlist[j]<<" and "<<inlist[k]<<" are likely to be correlated in "
	      <<nummatrix[k][j]<<" histograms"<<endl;
	  cout<<"   Minimal deviation="<<csmatrix[k][j]<<" sigma/bin in "<<namematrix[k][j]<<endl; 
	  report=1;
	}
      }
    }
    if (report==0) cout<<" ...passed"<<endl;cout<<endl;
    for (size_t j=0;j<inlist.size();j++) delete[] csmatrix[j];
    for (size_t j=0;j<inlist.size();j++) delete[] nummatrix[j];
    for (size_t j=0;j<inlist.size();j++) delete[] namematrix[j];
    delete[] csmatrix;
    delete[] nummatrix;
    delete[] namematrix;
  }

  int sc=0;
  for (i=0;i<filelist.size();i++) {
    Histogram h0(inlist[0]+filelist[i]);
    if (h0.Fills()>0) {
      bool valid=1;
      if (mode==0) h0.Restore();
      for (int j=1;j<inlist.size();j++) {
	Histogram histo(inlist[j]+filelist[i]);
	if (histo.Fills()>0) {
	  if (mode==0) histo.Restore();
	  switch (mode) {
	  case 2: 
	    h0.Addopt(histo);
	    break;
	  case 3:
	    h0.BinMin(histo);
	    break;
	  case 4:
	    h0.BinMax(histo);
	    break;
	  default: h0+=histo;
	  }
	}
	else valid=0;
      }
      if (valid) {
	if (mode==0) h0.Finalize();
	if (mode==5) h0.Scale(factor);
	if (noheader) h0.SetFills(-1);
	h0.Output(output+filelist[i]);
	sc++;
      }
    }
  }
  if (sc>0) cout<<"successfully combined "<<inlist.size()
		<<" directories containing "<<sc<<" histograms"<<endl; 
  delete ATOOLS::msg;
  delete ATOOLS::exh;
#ifdef USING__MPI
  MPI::Finalize();
#endif
  return 0;
}
 
